/*
This code simulates a quantum particle on a circonference.
It uses the metropolis algorithm to try an update on the euclidean time lattice of the position.
{x} = {x_0 , ..., x_{N_t}}

Every point in the lattice is a different (euclidean) time.

The idea is to fix the temperature and gradually increase the N_t.
STEP in the code and in the theory indicates the spacing of the sites and is defined as follows:

STEP = beta * h^2 / (4 pi m R^2) * (1 / N_t)

The values used in the code are:
DELTA = 0.5
ETA = 2 / N_t

basically in main we call the sweep function N_SWEEPS times and that function calls
the update function N_SAMPLES times. This is equal to say that we call the site update
N_SAMPLES * N_SWEEPS times. The use of sweeps is fictitious and with the sole scope 
of acquisition of data.  
*/


/*======================================
 Libraries 
========================================*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pcg32.h>
#include <time.h>
#include <stdbool.h>
#include <errno.h>


/*======================================
Global variables definitions 
========================================*/


#include "parameters.h"										/*From this file we get the N_t lattice dimention input, python script modifies the value
																						in the file. We also get the BETA and DELTA values ecc.*/
#define MAX_RAND 4294967296.0
#define WINDOW_SIZE 8


double lattice[N_t];
long long  proposed_updates;
long long accepted_updates;


/*======================================
	Functions definition
========================================*/


void initialize_lattice_cold_start(){
	for (int i = 0; i < N_t  ;  i++){
		lattice[i]  = 0.0;
	}
}


double distance(double x, double y){
	double d = x - y;
	if (d > 0.5) d -= 1.0;
	if (d < -0.5) d += 1.0;
	
	return d;	
}


//updating a single site
void site_update(pcg32_random_t* rng_pointer){
	int i = pcg32_boundedrand_r(rng_pointer, N_t);
	
	int next_site = i + 1;
  if (next_site >= N_t) next_site = 0;
    
  int previous_site = i - 1;
  if (previous_site < 0) previous_site = N_t - 1;
	
	double r = pcg32_random_r(rng_pointer) / MAX_RAND;
	double tentative_i = lattice[i] + (1.0 - 2.0 * r) * DELTA;
	
	if (tentative_i > 1.0) tentative_i -= 1.0;
	if (tentative_i < 0.0) tentative_i += 1.0;
	
	double d_new_next = distance(tentative_i, lattice[next_site]);
	double d_new_prev = distance(lattice[previous_site], tentative_i);
	
	double d_old_next = distance(lattice[i], lattice[next_site]);
	double d_old_prev = distance(lattice[previous_site], lattice[i]);
	
	double delta_S = PREFACTOR * (d_new_next * d_new_next + d_new_prev * d_new_prev - d_old_next * d_old_next - d_old_prev * d_old_prev);

	proposed_updates++;
	if (delta_S <= 0.0 || (pcg32_random_r(rng_pointer) / MAX_RAND) < exp(- delta_S) ){
		lattice[i] = tentative_i;
		accepted_updates++;
	}
}


//a fixed number of update proposals
void lattice_sweep(pcg32_random_t* rng_pointer){
	for (int i = 0; i < N_SAMPLES; i++){
		site_update(rng_pointer);
	}
}


//a fixed number of sweeps to throw away
void termalization(pcg32_random_t *rng_pointer){
	for(int i = 0; i < N_THERMAL; i++){
		lattice_sweep(rng_pointer);
	}
}


//number of times particle goes around the ring, although on lattice it's still a number in Z
double winding_number(){
	double Q = 0.0;
	
	for (int i = 0; i < N_t; i++){
		int next_site = (i + 1) % N_t;
		Q += distance(lattice[i], lattice[next_site]);
	}
	return Q;
}


// function for evaluation of <Q^2>, ensamble avarage of square winding number
double mean_winding_number2(const char* data){
	FILE* file = fopen(data, "rb");

	double sum = 0.0;
	double q2 = 0.0;
	int counter = 0;
	while (fread(&q2, sizeof(double), 1, file) == 1){
		sum += q2;
		counter++;
	}
	fclose(file);

	return sum / (double)counter;
}


// function for error of <Q^2>, jakeknife method: removing one block at a time 
double error_winding_number2(const char* data, int block_size){

	int blocks_number = N_SWEEPS / block_size;
	long long  total_valid_points = (long long) blocks_number * block_size;
	 
	double* block_sums = (double*)malloc(blocks_number * sizeof(double));
	
	double total_sum = 0.0;
	double q2 = 0.0;
	
	FILE* file = fopen(data, "rb");
  if (!file){
  	perror("Error opening measurements file");
  	exit(1);
  }

  //optimized computation of total sum and block sums in only one file reading
  //the file reader pointer just goes on from where it was the last for cycle, this way we divide in blocks
  for (int i = 0; i < blocks_number; i++){
  	double current_block_sum = 0;
  	for (int j = 0; j < block_size; j++){
  		if (fread(&q2, sizeof(double), 1, file) != 1){
  			perror("Read error in total sum!");
  			exit(1);
  		}
  		current_block_sum += q2;
  	}
  	block_sums[i] = current_block_sum;
  	total_sum += current_block_sum;
  }
  fclose(file);

	//Calculate Jackknife average
	double jack_avarage = 0.0;
  for (int i = 0; i < blocks_number; i++) {
  	// (Total Sum - This Block) / (All valid points - block_size)
   	double jack_sample = (total_sum - block_sums[i]) / (double)(total_valid_points - block_size);
  	jack_avarage += jack_sample;
  }
  jack_avarage /= blocks_number;

	//Calculate Variance
	double variance = 0.0;
	for (int i = 0; i < blocks_number; i++) {
  	double jack_sample = (total_sum - block_sums[i]) / (double)(total_valid_points - block_size);
  	double diff = jack_sample - jack_avarage;
  	variance += diff * diff;
  }

	//Final error calculation
  double error = sqrt(((double)(blocks_number - 1) / blocks_number) * variance);  
  free(block_sums);
  return error;
}


//evaluating the variance of the errors in the given window, if it's under threshold the block_size is correct!
double window_variance(double *data, int n, double *mean_out){
    double mean = 0.0;
    for(int i = 0; i < n; i++){
    	mean += data[i];
    }
    mean /= n;
    
    double variance = 0.0;
    for(int i = 0; i < n; i++){
    	double diff = data[i] - mean;
    	variance += diff * diff;
    }
    variance /= n;
    
    if (mean_out != NULL) {
    	*mean_out = mean;
    }
    
    return variance;
}


/*======================================
		Main 
========================================*/

int main(){

	//random number generator
	pcg32_random_t rng = {0}; 
	uint64_t seed1 = (uint64_t)time(NULL);
	uint64_t seed2 = (uint64_t)clock();
	pcg32_srandom_r(&rng, seed1, seed2);
	
	initialize_lattice_cold_start();
	termalization(&rng);
	
	//fixed temporary file for Q² values (overwritten each run, deleted at end)
	const char *measurements_Q2 = "data1/measurements.bin";
	FILE *f_measurements = fopen(measurements_Q2, "wb");

	//mesuring <Q^2> by evolving the lattice
	for(int sweep = 0; sweep < N_SWEEPS; sweep++){
		lattice_sweep(&rng);
		double q = winding_number();
		double q2 = q * q;
		fwrite(&q2, sizeof(double), 1, f_measurements);
	}
	fclose(f_measurements);
	
	printf("Acceptance rate: %.2f%%\n", 100.0 * accepted_updates / (double)proposed_updates);
	

	/*======================================
			mesuring error of <Q^2> 
	========================================*/
	double current_error = 0.0;
	int optimal_block_size = N_SWEEPS/100; //initialize not zero in case is not found a value and the code uses the starting one
	double recent_errors[WINDOW_SIZE];
	bool plateau_found = false;
	double cv_threshold = 0.001;
	
	
	//file with  first column: "block_size" ; second column: "current_error"
	char plateau_name[64];
	sprintf(plateau_name, "data1/plateau_%.4f.txt", (float)BETA);

	//cycles until the plateau is found or the threshold value is too high, every time the threshold increases
	while(!plateau_found && cv_threshold <= 0.02){
		int window_index = 0;

		FILE *f_plateau = fopen(plateau_name, "w");		//opening file in cycle so that it resets and in putput we only get the one where plateau is found
			
		for (int block_size = 10; block_size <= N_SWEEPS/100; block_size = (int)(block_size * 1.1)){
			
			current_error = error_winding_number2(measurements_Q2, block_size);
				
			fprintf(f_plateau, "%d %.10f\n", block_size, current_error);

			recent_errors[window_index % WINDOW_SIZE] = current_error;
			window_index++;

			if(window_index >= WINDOW_SIZE && block_size > MINIMAL_BLOCK_SIZE){

				double mean;
				double variance = window_variance(recent_errors, WINDOW_SIZE, &mean);
				double coefficient_variation = sqrt(variance) / mean;

				if (coefficient_variation < cv_threshold){
					optimal_block_size = block_size;
					plateau_found = true;
					break;
				}
			}
		}
		
		fclose(f_plateau);
		
		if (!plateau_found){
			cv_threshold += 0.001;
		}
	}

	if (optimal_block_size == N_SWEEPS/100) {
		printf("Warning: No plateau found for N_t=%d. Using fallback block size.\n", N_t);
	}else{
		printf("Plateau found at  block size: %d,  cv_threshold: %.4f \n", optimal_block_size, cv_threshold);
	}
	
	current_error = error_winding_number2(measurements_Q2, optimal_block_size);
	
	//writing in file 1:"N_t"  2: "<Q^2>"  3:"<Q^2> error"  
	FILE *f_results = fopen("data1/results.txt", "a");
	fprintf(f_results, "%d %.5f  %.8f %.8f\n", N_t, BETA,  mean_winding_number2(measurements_Q2), current_error);

	fclose(f_results);
	remove(measurements_Q2);
	
	return 0;
}
