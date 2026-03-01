#include <stdio.h>
#include <math.h>
#include <pcg32.h>
#include <time.h>

#define L 160
#define N_SAMPLES 60000
#define N_TERMAL 5000

int lattice[L][L];
float m[N_SAMPLES];
int stack_i[L * L];
int stack_j[L * L];
long running_sum =0;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



void initialize_lattice(pcg32_random_t* rng_pointer){
	running_sum = 0; 	//resetting running sum
	for (int i = 0;  i < L; i++){
		for (int j = 0; j < L; j++){
			uint32_t r = pcg32_boundedrand_r(rng_pointer, 2);   //random integer from 0 to 1
	  	lattice[i][j] = (r * 2) - 1; 	
	  	running_sum += lattice[i][j];
	 	}
	}
}



int wolff_step(pcg32_random_t* rng_pointer, double beta){
	float P_add = 1.0 - exp(-2.0 * beta);
	
	int stack_pointer = 0;
	int cluster_size = 0;
	
	int seed_i =  pcg32_boundedrand_r(rng_pointer, L);
	int seed_j =  pcg32_boundedrand_r(rng_pointer, L);
	
	int target_spin = lattice[seed_i][seed_j];
	
	//update running sum
	running_sum -= (lattice[seed_i][seed_j] * 2);
	//flip spin immediately
	lattice[seed_i][seed_j] *= -1; 
	
	//save coordinates to the stack
	stack_i[stack_pointer] = seed_i;
	stack_j[stack_pointer] = seed_j;
	
	//how many entries are filled in the stack = stack_pointer
	stack_pointer++;
	cluster_size++;
	
	while(stack_pointer > 0){
		stack_pointer--;
		int c_i = stack_i[stack_pointer];
		int c_j = stack_j[stack_pointer];
		
		//coordinates of the first neighbors
		// change 'ci' (row) while keeping 'cj' (column/horizontal) steady
		int ni[4] = {(c_i == 0) ? L - 1 : c_i - 1, (c_i == L - 1) ? 0 : c_i + 1, c_i, c_i};
		// keep 'ci' steady while changing 'cj'
		int nj[4] = {c_j, c_j, (c_j == 0) ? L - 1 : c_j - 1, (c_j == L - 1) ? 0 : c_j + 1};	
		
		for (int k = 0; k < 4; k++){
			int i_neigh = ni[k];
			int j_neigh = nj[k];
			
			if(lattice[i_neigh][j_neigh] == target_spin){
				if((pcg32_random_r(rng_pointer)/4294967296.0) < P_add){
				
					//update running sum_
					running_sum -= (lattice[i_neigh][j_neigh] * 2);
					
					//immediate flip
					lattice[i_neigh][j_neigh] *= -1;
						
					//push to the stack_i
					stack_i[stack_pointer] = i_neigh;
					stack_j[stack_pointer] = j_neigh;
					stack_pointer++;
					cluster_size++; 
				}
			}
		}
	}
	return cluster_size;	
}



void lattice_sweep(pcg32_random_t* rng_pointer, double beta){
	int total_flips = 0;
  while (total_flips < (L * L)){
  	total_flips += wolff_step(rng_pointer, beta);
  }
}



void termalization(pcg32_random_t* rng_pointer, double beta){
	for(int i = 0; i < N_TERMAL; i++){
		lattice_sweep(rng_pointer, beta);
	}
}



float magnetization_calculator(){
	return (float)running_sum / (L * L);
}



//<|m|> ensamble avarage of the abs_magnetization
float avarage_abs_magnetization(){
	float abs_msum = 0;
	for (int i = 0; i < N_SAMPLES; i++){
		abs_msum += fabs(m[i]); 
	}
	return (float) abs_msum / N_SAMPLES;
}


//error of <|m|> with blocking
float error_function_avg_abs_magnetization(){
	int blocks_number = 80;
	int block_size = N_SAMPLES / blocks_number;
	double block_avarage[blocks_number];
	double total_sum = 0;
	double total_sum_sq = 0;
	
	for (int i = 0; i < blocks_number ; i++){
		double block_sum = 0;
		for (int j = 0; j < block_size; j++){
			block_sum += fabs(m[i * block_size + j]);
		}
		
		block_avarage[i] = block_sum / block_size;
		total_sum += block_avarage[i]; 
	}
	
	double grand_avarage = total_sum / blocks_number; 
	for (int i = 0; i < blocks_number; i++){
		double difference = block_avarage[i] - grand_avarage;
		total_sum_sq += difference * difference;
	}

	double variance_of_blocks = total_sum_sq / (blocks_number - 1);
	
	//sigma = sqrt(Var(blocks) / number of blocks)
	return (float) sqrt(variance_of_blocks / blocks_number);
}



// susceptibility Chi = L^2 * <m^2>   
float susceptibility(){
	float m2_sum = 0;
	for(int i = 0; i < N_SAMPLES; i++){
		m2_sum += (double) m[i] * m[i];
	}
	m2_sum = (float)(m2_sum * (L * L));
	return (float)(m2_sum / N_SAMPLES);
}



// Subtracted Susceptibility: L^2 * (<m^2> - <|m|>^2)
float subtracted_susceptibility(){
	double m2_sum = 0;
	double abs_m_sum = 0;
    
		for(int i = 0; i < N_SAMPLES; i++){
    	double val = (double)m[i];
     	m2_sum += val * val;
      abs_m_sum += fabs(val);
    }
    
  double avg_m2 = m2_sum / N_SAMPLES;
	double avg_abs_m = abs_m_sum / N_SAMPLES;

	float chi_sub = (float)((L * L) * (avg_m2 - (avg_abs_m * avg_abs_m)));
    
  return chi_sub;
}


//jackknife for error evaluation of subtracted susceptibility
float error_function_subtracted_susceptibility(){
	int blocks_number = 80;
	int block_size = N_SAMPLES / blocks_number;
	double block_m2[blocks_number], block_abs_m[blocks_number], jack_chip[blocks_number];
	double global_m2 = 0, global_abs = 0;
	
	for (int i = 0; i < blocks_number; i++){
		double b_m2 = 0, b_abs = 0;
		for (int j = 0; j < block_size; j++){
			double val = (double) m[i*block_size + j];
			b_m2 += val * val;
			b_abs += fabs(val);
		}
		block_m2[i] = b_m2 / block_size;
		block_abs_m[i] = b_abs / block_size;
		global_m2 += block_m2[i];
		global_abs += block_abs_m[i];
	}
	
	//resampling (leaving one out)
	for (int i = 0; i < blocks_number; i++){
		double jack_m2 = (global_m2 - block_m2[i])/(blocks_number - 1);
		double jack_abs = (global_abs - block_abs_m[i])/(blocks_number - 1);
		jack_chip[i] = (double)(L * L) * (jack_m2 - (jack_abs * jack_abs));
	}
	
	//variance of jackknife samples
	double avarage_jack = 0;
	for (int i = 0; i < blocks_number; i++) avarage_jack += jack_chip[i];
	avarage_jack /= blocks_number;
	
	double sum_sq_diff = 0;
  	for (int i = 0; i < blocks_number; i++) {
    	double diff = jack_chip[i] - avarage_jack;
    	sum_sq_diff += diff * diff;
    }
    
	//Jackknife error formula: sqrt( [(K-1)/K] * sum of squared differences )
  return (float)sqrt(((double)(blocks_number - 1) / blocks_number) * sum_sq_diff);
}



//Binder cumulant: U = <m^4> / <m^2>^2
float binder_cumulant(){
	double m4_sum = 0;
	double m2_sum = 0;
	for (int i = 0; i < N_SAMPLES; i++){
		double val = m[i];
		m4_sum += (double) val * val * val * val;
		m2_sum += (double) val * val;
	}
	double avg_m2 = m2_sum / N_SAMPLES;
	double avg_m4 = m4_sum / N_SAMPLES;
	
	return (float) (avg_m4)/(avg_m2 * avg_m2); 
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


int main(){
	
	clock_t start = clock();
	
	pcg32_random_t rng = {0}; //defining the rng variable

	//seeds
	uint64_t seed1 = (uint64_t)time(NULL);
  uint64_t seed2 = (uint64_t)clock();
    

	//seeding function for the number generator
	pcg32_srandom_r(&rng, seed1, seed2);
	
	//initialize temperature
	double beta = 0.433;
	
	FILE *f = fopen("magbeta/magbeta160.txt", "w");
	FILE *g = fopen("chibeta/chibeta160.txt", "w");
	FILE *h = fopen("chipbeta/chipbeta160.txt", "w");
	FILE *l = fopen("ubeta/ubeta160.txt", "w");
	
	
	//cycle to take the <|m|> and increase beta, fixed L
	for (int i = 0; i < 40; i++){
	
		initialize_lattice(&rng);
		termalization(&rng, beta);
		
		for (int j = 0; j < N_SAMPLES; j++){
			lattice_sweep(&rng, beta);
			m[j] = magnetization_calculator();
		}
		
		fprintf(f, " %.5f , %.5f, %.5f \n", beta, avarage_abs_magnetization(), error_function_avg_abs_magnetization());
		fprintf(g, " %.5f , %.5f \n", beta, susceptibility());
		fprintf(h, " %.5f , %.5f , %.5f \n", beta, subtracted_susceptibility(), error_function_subtracted_susceptibility());
		fprintf(l, " %.5f , %.5f \n", beta, binder_cumulant());
		
		beta += 0.00025;
	}
	
	fclose(f);
	fclose(g);
	fclose(h);
	fclose(l);
	
	//CPU time
	clock_t end = clock();
	double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf(" Time: %.3f seconds \n", cpu_time);
	
	return 0;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















