/*
This code simulates the free scalar field in the contest of quantum field theory.
To achive such thing is needed a lattice discretization both  in the (D-1) spatial directions and in the time direction, with some algebra and going back to the continuous
we get to the partition function that contains the euclidean action, in short we get the path integral formulation.

Now we discretize the scalar field again in spatial and time directions in the context of path integral formulation in order to to simulate it.

Continuous limit:
Once defined the wall-wall operators to perform the continuous limit we need N_\mu ---> +inf

Lattice of size : N_t x N_s
temperature : N_t x a 
every site : phi_n 
evaluation : = <O_1>,<O_2>,<O_3>  ===> find epsilon ===> integrating find pressure
update method : 20% heat-bath, 80% microcanonical
boundary conditions : periodic both on time and space
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


/*======================================
	Global variables definitions 
========================================*/


#include "parameters.h"

//parameters fixed
#define N_SWEEPS 50000000
#define N_TERMAL 500000
#define MAX_RAND 4294967296.0
#define HEATH_BATH_PROBABILITY 0.2

#define M_HAT_SQUARED (M_HAT * M_HAT)
#define VARIANCE_DENOMINATOR (M_HAT_SQUARED + 4.0)
#define INVERSE_VARIANCE_DENOMINATOR (1 / VARIANCE_DENOMINATOR)
#define SIGMA sqrt(1.0 / VARIANCE_DENOMINATOR)
#define INVERSE_LATTICE_VOLUME (1.0 / (N_T * N_S))


double lattice[N_T][N_S];



/*======================================
			Functions to update sites
	20% heat-bath / 80% microcanonical
========================================*/


//See notebook for theory behind the heat-bath update 
void heat_bath_update(int time_index, int space_index, pcg32_random_t* random_number){

	double sum_of_neighbors;
	double gaussian_mean;
	double random_gaussian_number;
	double new_field_value;

	int time_index_forward = (time_index + 1) % N_T;
	int time_index_backward = (time_index - 1 + N_T) % N_T;
	int space_index_forward = (space_index + 1) % N_S;
	int space_index_backward = (space_index - 1 + N_S) % N_S;

	sum_of_neighbors = lattice[time_index_forward][space_index] + lattice[time_index_backward][space_index] + lattice[time_index][space_index_forward] + lattice[time_index][space_index_backward];

	gaussian_mean = (sum_of_neighbors * INVERSE_VARIANCE_DENOMINATOR);

	//box-muller agorithm, change of variables drom two uniformely distributed random numbers to one with gaussian distribution
	double uniform_random1;
	do{
		uniform_random1 = pcg32_random_r(random_number) / MAX_RAND;
	}while (uniform_random1 == 0.0);
	double uniform_random2 = pcg32_random_r(random_number) / MAX_RAND;
	random_gaussian_number = sqrt(-2.0 * log(uniform_random1)) * cos(2.0 * M_PI * uniform_random2);

	new_field_value = gaussian_mean + SIGMA * random_gaussian_number;
	lattice[time_index][space_index] = new_field_value;
}


//See notebook for notes on the microcanonical update
void microcanonical_update(int time_index, int space_index){

	double sum_of_neighbors;
	double old_field_value;
	double new_field_value;

	int time_index_forward = (time_index + 1) % N_T;
	int time_index_backward = (time_index - 1 + N_T) % N_T;
	int space_index_forward = (space_index + 1) % N_S;
	int space_index_backward = (space_index - 1 + N_S) % N_S;

	sum_of_neighbors = lattice[time_index_forward][space_index] + lattice[time_index_backward][space_index] + lattice[time_index][space_index_forward] + lattice[time_index][space_index_backward];

	old_field_value = lattice[time_index][space_index];

	new_field_value = 2 * (sum_of_neighbors * INVERSE_VARIANCE_DENOMINATOR) - old_field_value;
	lattice[time_index][space_index] = new_field_value;
}


void site_update(pcg32_random_t* random_number){
	int random_time_index;
	int random_space_index;
	double update_method_selector;

	random_time_index = (int)((pcg32_random_r(random_number) / MAX_RAND) * N_T);
	random_space_index = (int)((pcg32_random_r(random_number) / MAX_RAND) * N_S);

	update_method_selector = pcg32_random_r(random_number) / MAX_RAND ;

	if (update_method_selector < HEATH_BATH_PROBABILITY){
		heat_bath_update(random_time_index, random_space_index, random_number);
	} else {
		microcanonical_update(random_time_index, random_space_index);
		}
}



/*======================================
			observables functions
========================================*/

double O1_operator(){
	double operator_sum = 0.0;
	double field_value;
	
	for (int time_index = 0; time_index < N_T; time_index++){
		for (int space_index = 0; space_index < N_S; space_index++){
			
			field_value = lattice[time_index][space_index];
			operator_sum += field_value * field_value;
		}
	}
	return (double) (M_HAT_SQUARED * operator_sum) * INVERSE_LATTICE_VOLUME;
}


double O2_operator(){
	double operator_sum = 0.0;
	double field_difference; 
	int space_index_forward;
	
	for (int time_index = 0; time_index < N_T; time_index++){
		for (int space_index = 0; space_index < N_S; space_index++){
			
			space_index_forward = (space_index + 1) % N_S;
			field_difference = lattice[time_index][space_index_forward] - lattice[time_index][space_index];
			operator_sum += field_difference * field_difference;
		}
	}
	return (double) operator_sum * INVERSE_LATTICE_VOLUME; 
}


double O3_operator(){
	double operator_sum = 0.0;
	double field_difference;
	int time_index_forward;

	for (int time_index = 0; time_index < N_T; time_index++){
		for (int space_index = 0; space_index < N_S; space_index++){
			
			time_index_forward = (time_index + 1) % N_T;
			field_difference = lattice[time_index_forward][space_index] - lattice[time_index][space_index];
			operator_sum += field_difference * field_difference;
		}
	}
	return (double) operator_sum * INVERSE_LATTICE_VOLUME;
}



/*======================================
				General purpose	functions
========================================*/


void initialize_lattice_cold_start(){
	for(int t = 0; t < N_T; t++){
		for(int x = 0; x < N_S; x++){
			lattice[t][x] = 0.0;
		}
	}
}


void initialize_lattice_hot_start(pcg32_random_t* random_number){
	for(int t = 0; t < N_T; t++){
		for(int x = 0; x < N_S; x++){
			//use box-muller to generate the gaussian lattice values
			double uniform_random1 = pcg32_random_r(random_number) / MAX_RAND;
			double uniform_random2 = pcg32_random_r(random_number) / MAX_RAND;
			double random_gaussian_number = sqrt(-2.0 * log(uniform_random1)) * cos(2.0 * M_PI * uniform_random2);
			lattice[t][x] = SIGMA * random_gaussian_number;
		}
	}
}


void lattice_sweep(pcg32_random_t* random_number){
	int total_updates = N_T * N_S;

	for (int update_number = 0; update_number < total_updates; update_number++){
		site_update(random_number);
	}
}


void termalization(pcg32_random_t* random_number){
	for (int index = 0; index < N_TERMAL; index++){
		lattice_sweep(random_number);
	}
	
}



/*======================================
			Analysis and error functions
========================================*/


double jack_knife_error (double* data, int total_data, int block_size){

	int number_of_blocks = total_data / block_size;
	double total_sum = 0.0;

	double* block_sums = (double*)malloc(number_of_blocks * sizeof(double));
	if (block_sums == NULL) {
		printf("Error: Memory allocation failed in jackknife_error\n");
		exit(1);
	} 
	
	for (int block_index = 0; block_index < number_of_blocks; block_index++){
	
		double block_sum = 0.0;
		int start_index = block_index * block_size;

		for (int i = 0; i < block_size; i++){
			block_sum += data[start_index + i];	
		}

		block_sums[block_index] = block_sum;
		total_sum += block_sum;
	}

	double effective_N = block_size * number_of_blocks;
	double jack_N = effective_N - block_size;

	double jack_mean_global = 0.0;

	for (int i = 0; i < number_of_blocks; i++){
	
		double jack_sample_mean = (total_sum - block_sums[i]) / jack_N;
		block_sums[i] = jack_sample_mean;
		jack_mean_global += jack_sample_mean;
	}

	jack_mean_global /= number_of_blocks;

	double variance = 0.0;
	for (int i = 0; i < number_of_blocks; i++){
	
		double difference = jack_mean_global - block_sums[i];
		variance += difference * difference; 
	}

	double error = sqrt(((double)(number_of_blocks - 1) / number_of_blocks) * variance);
	free (block_sums);
	
	return error;
}



//finding the block size that realizes the plateau
int optimal_block_size(double* data, int total_data){

	double previous_error = 0.0;
	int stable_count = 0;
	int optimal_block_size = 100;
	
	for(int block_size = 100; block_size <= total_data / 50; block_size = (int)(block_size * 1.2)){

		double current_error = jack_knife_error(data, total_data, block_size);
		
		if (previous_error > 0.0){
			double relative_change = fabs(current_error - previous_error) / previous_error;
			if (relative_change < 0.01){
				stable_count++;
				if (stable_count >= 7){
					optimal_block_size = block_size;
					break;
				}
			} else {
				stable_count = 0;
			}
		}
		
		previous_error = current_error;
	}

	return optimal_block_size;
}



//function to obtain and save mean values and errors of the main observables
void compute_statistics(double* O1_operator_data, double* O2_operator_data, double* O3_operator_data, int total_measurements){

	double mean_O1 = 0.0, mean_O2 = 0.0, mean_O3 = 0.0;
	
	for (int i = 0; i < total_measurements; i++){
		mean_O1 += O1_operator_data[i];
		mean_O2 += O2_operator_data[i];
		mean_O3 += O3_operator_data[i];
	}

	mean_O1 /= total_measurements;
	mean_O2 /= total_measurements;
	mean_O3 /= total_measurements;

	int block_size_O1 = optimal_block_size(O1_operator_data, total_measurements);
	int block_size_O2 = optimal_block_size(O2_operator_data, total_measurements);
	int block_size_O3 = optimal_block_size(O3_operator_data, total_measurements);

	double error_O1 = jack_knife_error(O1_operator_data, total_measurements, block_size_O1);
	double error_O2 = jack_knife_error(O2_operator_data, total_measurements, block_size_O2);
	double error_O3 = jack_knife_error(O3_operator_data, total_measurements, block_size_O3);

	FILE* file_O1 = fopen("data/data_O1.txt", "a");
	FILE* file_O2 = fopen("data/data_O2.txt", "a");
	FILE* file_O3 = fopen("data/data_O3.txt", "a");

	if(file_O1 != NULL){
		fprintf(file_O1, "%.7f, %.9f, %d, %d \n", mean_O1, error_O1, N_T, N_S);
		fclose(file_O1);
	}
	
	if(file_O2 != NULL){
		fprintf(file_O2, "%.7f, %.9f, %d, %d \n", mean_O2, error_O2, N_T, N_S);
		fclose(file_O2);
	}
	
	if(file_O3 != NULL){
		fprintf(file_O3, "%.7f, %.9f, %d, %d \n", mean_O3, error_O3, N_T, N_S);
		fclose(file_O3);
	}
	
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

	initialize_lattice_hot_start(&rng);

	termalization(&rng);
	
	// allocate memory for storing operator measurements
	double* operator_O1_measurements = (double*)malloc(N_SWEEPS * sizeof(double));
	double* operator_O2_measurements = (double*)malloc(N_SWEEPS * sizeof(double));
	double* operator_O3_measurements = (double*)malloc(N_SWEEPS * sizeof(double));

	if(operator_O1_measurements == NULL || operator_O2_measurements == NULL || operator_O3_measurements == NULL){
		printf("Error: failed to allocate memory for operator measurements\n");
		return 1;
	}

	for (int sweep_index = 0; sweep_index < N_SWEEPS; sweep_index++){
		lattice_sweep(&rng);

		operator_O1_measurements[sweep_index] = O1_operator();
		operator_O2_measurements[sweep_index] = O2_operator();
		operator_O3_measurements[sweep_index] = O3_operator();	
	}

	compute_statistics(operator_O1_measurements, operator_O2_measurements, operator_O3_measurements, N_SWEEPS);
	
	free(operator_O1_measurements);
	free(operator_O2_measurements);
	free(operator_O3_measurements);

	return 0;
}

