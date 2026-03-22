/*
2-D Ising model with cluster update: find the critical exponent values
Lattice sizes = 40, 60, 80, 100, 120, 140, 160, 180
10^6 single cluster updates for every beta, 40 beta for every lattice size
N_SWEEPS = tbd
N_THERMAL = tbd

observables: (all avarage observables)
energy densitity = e 
absolute magnetization
squared magnetization (susceptibility)
[susceptibility] - [squared absolute magnetization] (subtracted susceptibility)
Binder cumulant



run the code:
we will use a python script to run the code
*/



/*======================================
						 Libraries 
========================================*/


#include <stdio.h>
#include <math.h>
#include <pcg32.h>
#include <time.h>
#include <stdbool.h>
#include <stdlib.h>



/*======================================
			Global variables definition
========================================*/


#include "parameters.h"

#define MAX_RAND 4294967296.0
//number of times we sweep without acquisition
#define N_THERMAL 10000
//probability of adding the site to the cluster
#define P_ADD (1.0 - exp(-2.0 * BETA))

int lattice[L][L];
int stack_x[L * L];
int stack_y[L * L];




/*========================================
Functions to select and update the cluster
========================================*/


//does the spin flip? true or false return
bool spin_flip_function(int x, int y, pcg32_random_t* rng){
	if( (pcg32_random_r(rng)/MAX_RAND) < P_ADD){
		lattice[x][y] *= -1;
		return true;			
	}
	
	return false;
}


//generate and flip the cluster
void cluster_update(pcg32_random_t* rng){
	
	int stack_index = 0;
	int seed_x =  pcg32_boundedrand_r(rng, L);
	int seed_y =  pcg32_boundedrand_r(rng, L);
	int target_spin = lattice[seed_x][seed_y];

	lattice[seed_x][seed_y] *= -1;

	stack_x[stack_index] = seed_x;
	stack_y[stack_index] = seed_y;
	stack_index++;

	while(stack_index > 0){
		stack_index--;
		int coord_x = stack_x[stack_index];
		int coord_y = stack_y[stack_index];
		//left, right, up, down
		int neighbor_x[4] = { (coord_x - 1 + L)%L , (coord_x + 1)%L , coord_x , coord_x };
		int neighbor_y[4] = { coord_y , coord_y , (coord_y + 1)%L , (coord_y - 1 + L)%L };

		for(int i = 0; i < 4; i++){
			int x = neighbor_x[i];
			int y = neighbor_y[i];
			if(lattice[x][y] == target_spin){
				 if(spin_flip_function(x, y, rng)){
				 	stack_x[stack_index] = x;	
				 	stack_y[stack_index] = y;
				 	stack_index++;
				 }
			}
		}
			
	}
	
}



/*======================================
		Functions for prime observables
========================================*/


//magnetization of a given state of the lattice {s}: m
double magnetization(){
	int sum = 0;
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			sum += lattice[i][j];
		}	
	}
	return (double) sum  / (L*L);
}


//energy density of a given state of the lattice {s}: e
double energy_density(){
	int sum = 0;
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			sum += lattice[i][j] * lattice[(i+1)%L][j];
			sum += lattice[i][j] * lattice[i][(j+1)%L];
		}	
	}
	return (double) -sum / (L*L);
}



/*======================================
		General purpose functions
========================================*/


//random start for every lattice size
void initialize_lattice(pcg32_random_t* rng){
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			uint32_t r = pcg32_boundedrand_r(rng, 2);
			lattice[i][j] =  (r * 2) - 1;
		}	
	}	
}


//if we call it with NULL in the last two entries we are not measuring, just sweeping
void lattice_sweep(pcg32_random_t* rng, double* m, double* e){
	for(int i = 0; i < N_CLUSTER ; i++){
		cluster_update(rng);
	}
	if(m != NULL) *m = magnetization();
	if(e != NULL) *e = energy_density();
}


void thermalization(pcg32_random_t* rng){
	for(int i = 0; i < N_THERMAL; i++){
		lattice_sweep(rng, NULL, NULL);
	}
}



/*======================================
						observables
========================================*/


// <e> 
double energy(double* e_samples, int n){
	double sum = 0;
	for(int i = 0; i < n;  i++) sum += e_samples[i];
	return (double)sum/n;
}


// < e^2 > - < e >^2
double energy_variance(double* e_samples, int n){
	double sum_e = 0, sum_e2 = 0;
	for(int i = 0; i < n; i++){
		sum_e  += e_samples[i];
		sum_e2 += e_samples[i]*e_samples[i];
	}
	return (sum_e2/n) - (sum_e/n)*(sum_e/n);
}


// < |m| > 
double absolute_magnetization(double* m_samples, int n){
	double sum = 0;
	for(int i = 0; i < n;  i++) sum += fabs(m_samples[i]);
	return (double)sum/n;
}


// < m^2 >
double magnetization2(double* m_samples, int n){
	double sum = 0;
	for(int i = 0; i < n; i++ ) sum += m_samples[i] * m_samples[i];
	return sum/n;
}


// < m^2 > - < |m| >^2 
double absolute_magnetization_variance(double* m_samples, int n){
	double mag2_avg = magnetization2(m_samples, n);
	double abs_mag_avg = absolute_magnetization(m_samples, n);
	return mag2_avg - (abs_mag_avg * abs_mag_avg);
}


// U = < m^4 > / < m^2 >^2 
double binder(double* m_samples, int n){
	double sum_m2 = 0, sum_m4 = 0;
	for(int i = 0; i < n; i++){
		double m2 = m_samples[i]*m_samples[i];
		sum_m2 += m2;
		sum_m4 += m2 * m2; 
	}
	double avg_m2 = (double)sum_m2/n;
	double avg_m4 = (double)sum_m4/n;
	return (double) avg_m4 / (avg_m2 * avg_m2);
}



/*======================================
		Error finding general functions 
========================================*/


//gives standard deviation given the data and the size of the array passed
double standard_deviation(double* data, int n){
	double sum = 0, sum2 = 0;
	 for (int i = 0; i < n; i++){
	 	sum += data[i];
    sum2 += data[i] * data[i];
  }
  double avg = sum / n;
  double var = (sum2 / n) - (avg * avg);
	return sqrt(var > 0 ? var : 0);
}


//given a block-size and the data gives back the error evaluated with blocking
double blocking_error(double* data, int block_size, int n){
	int blocks_number = n/block_size;
	double* block_avarages = malloc(blocks_number*sizeof(double));
	
	for(int i=0; i < blocks_number; i++){
		double sum = 0;
		
		for(int j=0; j<block_size; j++){
			sum += data[(block_size * i) + j];
		}
		block_avarages[i] =  sum / (double)block_size;
	}

	double error = standard_deviation(block_avarages, blocks_number)/(double)(sqrt(blocks_number));
	free(block_avarages);
	return error;
}


//gives back the error at the block-size that realizes the plateau
double plateau_finder(double* data, int n){
	double previous_error = -1.0;
  double current_error  =  0.0;

  for(int block_size = 5; block_size <= n/40; block_size += 5){
  	current_error = blocking_error(data, block_size, n);
    if(previous_error > 0 && fabs(current_error - previous_error) / previous_error < 0.05){
    	return current_error;
    }
    previous_error = current_error;
  }

	fprintf(stderr, "WARNING: plateau not found in blocking (n=%d), returning last value %f\n", n, current_error);
	return current_error;
}


//takes the block size and the data array and gives back the jack-knife error
double jackknife_error_binder(double* data, int n, int block_size){
	int blocks_number = n / block_size;
	int effective_n = blocks_number * block_size;
	double full_sum_m2 = 0, full_sum_m4 = 0;
	double* block_sum_m2 = malloc(blocks_number * sizeof(double));
	double* block_sum_m4 = malloc(blocks_number * sizeof(double));

	for (int i = 0; i < blocks_number; i++) {
  	double s2 = 0, s4 = 0;
  	for (int j = 0; j < block_size; j++) {
   		double m2 = data[i * block_size + j] * data[i * block_size + j];
   		s2 += m2;
   		s4 += m2 * m2;
  	}
  	block_sum_m2[i] = s2;
  	block_sum_m4[i] = s4;
  	full_sum_m2 += s2;
  	full_sum_m4 += s4;
	}
 
 double full_avg_m2 = full_sum_m2 / effective_n;
 double full_avg_m4 = full_sum_m4 / effective_n;
 double full_estimate = full_avg_m4 / (full_avg_m2 * full_avg_m2);
 
 double sum_diff_sq = 0;
 int reduced_n = effective_n - block_size;

 for (int i = 0; i < blocks_number; i++) {
  double res_m2 = (full_sum_m2 - block_sum_m2[i]) / reduced_n;
  double res_m4 = (full_sum_m4 - block_sum_m4[i]) / reduced_n;
  double jack_estimate = res_m4 / (res_m2 * res_m2);
  
  double diff = jack_estimate - full_estimate;
  sum_diff_sq += diff * diff;
 }

 free(block_sum_m2); free(block_sum_m4);
 return sqrt((double)(blocks_number - 1) / blocks_number * sum_diff_sq);
}



double jackknife_error_variance(double* data, int n, int block_size) {
 int blocks_number = n / block_size;
 int effective_n = blocks_number * block_size;
 double full_sum = 0, full_sum_sq = 0;
 double* b_sum = malloc(blocks_number * sizeof(double));
 double* b_sum_sq = malloc(blocks_number * sizeof(double));

 for (int i = 0; i < blocks_number; i++) {
  double s = 0, s2 = 0;
  for (int j = 0; j < block_size; j++) {
   s += data[i * block_size + j];
   s2 += data[i * block_size + j] * data[i * block_size + j];
  }
  b_sum[i] = s; b_sum_sq[i] = s2;
  full_sum += s; full_sum_sq += s2;
 }

 double full_est = (full_sum_sq / effective_n) - pow(full_sum / effective_n, 2);
 double sum_diff_sq = 0;
 int red_n = effective_n - block_size;

 for (int i = 0; i < blocks_number; i++) {
  double res_avg = (full_sum - b_sum[i]) / red_n;
  double res_avg_sq = (full_sum_sq - b_sum_sq[i]) / red_n;
  double jack_est = res_avg_sq - (res_avg * res_avg);
  
  double diff = jack_est - full_est;
  sum_diff_sq += diff * diff;
 }

 free(b_sum); free(b_sum_sq);
 return sqrt((double)(blocks_number - 1) / blocks_number * sum_diff_sq);
}

 



/*======================================
	Error functions of the observables
========================================*/


//errors with blocking, primary observables


double energy_error(double* e_samples, int n){
	return plateau_finder(e_samples, n);
}


double absolute_magnetization_error(double* m_samples, int n){
	double* abs_m = malloc(n * sizeof(double));
	for(int i = 0; i < n; i++) abs_m[i] = fabs(m_samples[i]);
	double error = plateau_finder(abs_m, n);
	free(abs_m);
	return error;
}


double magnetization2_error(double* m_samples, int n){
	double* m2 = malloc(n * sizeof(double));
	for(int i = 0; i < n; i++) m2[i] = m_samples[i] * m_samples[i];
	double error = plateau_finder(m2, n);
	free(m2);
	return error;
}


//errors with jack-knife, secondary observables


//error with plateau finding logic, the relative change must be less than 3%
double binder_error(double* m_samples, int n){
 double previous_step_error = -1.0;
 double current_step_error = 0.0;
 
 for(int block_size = 5; block_size <= n / 40; block_size += 5){
  current_step_error = jackknife_error_binder(m_samples, n, block_size);
  
  if(previous_step_error > 0){
   double relative_change = fabs(current_step_error - previous_step_error) / previous_step_error;
   if(relative_change < 0.05){
    return current_step_error;
   }
  }
  previous_step_error = current_step_error;
 }
 fprintf(stderr, "WARNING: plateau not found in binder jackknife (n=%d), returning %f\n", n, current_step_error);
 return current_step_error;
}


double absolute_magnetization_variance_error(double* m_samples, int n){
 double previous_step_error = -1.0;
 double current_step_error = 0.0;
 
 double* absolute_magnetization_values = malloc(n * sizeof(double));
 for (int i = 0; i < n; i++) {
  absolute_magnetization_values[i] = fabs(m_samples[i]);
 }

 for (int block_size = 5; block_size <= n / 40; block_size += 5) {
  current_step_error = jackknife_error_variance(absolute_magnetization_values, n, block_size);
  
  if (previous_step_error > 0) {
   double relative_change = fabs(current_step_error - previous_step_error) / previous_step_error;
   if (relative_change < 0.05) {
    free(absolute_magnetization_values);
    return current_step_error;
   }
  }
  previous_step_error = current_step_error;
 }
 
 fprintf(stderr, "WARNING: plateau not found in magnetization variance jackknife (n=%d), returning %f\n", n, current_step_error);
 free(absolute_magnetization_values);
 return current_step_error;
}


double energy_variance_error(double* e_samples, int n){
 double previous_step_error = -1.0;
 double current_step_error = 0.0;
 
 for (int block_size = 5; block_size <= n / 40; block_size += 5) {
  current_step_error = jackknife_error_variance(e_samples, n, block_size);
  
  if (previous_step_error > 0) {
   double relative_change = fabs(current_step_error - previous_step_error) / previous_step_error;
   if (relative_change < 0.05) {
    return current_step_error;
   }
  }
  previous_step_error = current_step_error;
 }
 fprintf(stderr, "WARNING: plateau not found in energy variance jackknife (n=%d), returning %f\n", n, current_step_error);
 return current_step_error;
}

/*======================================
									MAIN
========================================*/


int main(){

	//random number generator
	pcg32_random_t rng = {0};
	uint64_t seed1 = (uint64_t)time(NULL);
	uint64_t seed2 = (uint64_t)clock();
	pcg32_srandom_r(&rng, seed1, seed2);


	//pre-analysis functions
	initialize_lattice(&rng);
	thermalization(&rng);


	//allocate arrays for the raw sampling
	int n_samples = N_SWEEPS / N_SKIP;
	double* m_samples = malloc(n_samples * sizeof(double));
	double* e_samples = malloc(n_samples * sizeof(double));


	//sweeping and acquisition of e and m every N_SKIP sweeps
	for(int i=0; i<N_SWEEPS; i++){
		if(i % N_SKIP == 0){
			lattice_sweep(&rng, &m_samples[i/N_SKIP], &e_samples[i/N_SKIP]);
		}else{
			lattice_sweep(&rng, NULL, NULL);	
		}
	}

	
	//writing in to a file L; beta; every operator value with error
	FILE* f = fopen("results/results.txt", "a");
	fprintf(f, "%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
    L,
    BETA,
    energy(e_samples, n_samples),                           energy_error(e_samples, n_samples),
    energy_variance(e_samples, n_samples),                  energy_variance_error(e_samples, n_samples),
    absolute_magnetization(m_samples, n_samples),                absolute_magnetization_error(m_samples, n_samples),
    absolute_magnetization_variance(m_samples, n_samples),  absolute_magnetization_variance_error(m_samples, n_samples),
    magnetization2(m_samples, n_samples),                   magnetization2_error(m_samples, n_samples),
    binder(m_samples, n_samples),                           binder_error(m_samples, n_samples)
	);
	fclose(f);


	free(m_samples);
	free(e_samples);
	return 0;
}
