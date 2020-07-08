//////////////////////////////////////////// Spin glass phase transition ////////////////////////////////////////

/* 

Phase space of the SK model.

Please cite the following paper when you use this code.

[Ezaki T, Fonseca dos Reis E, Watanabe T, Sakaki M, Masuda N. Closer to critical resting-state neural dynamics in individuals with higher fluid intelligence. Commun Biol 3:1 (2020).](https://www.nature.com/articles/s42003-020-0774-y)

Please do not distribute without contacting the authors above.

If you find a bug in this code, please contact the authors.

Author: Elohim Fonseca dos Reis

Date: July 8, 2020


Parameters:
----------
   This program maps the SK model phase space using the following parameters that have to be set by the user before compiling:
	- total number of spins (N);
	- thermal average dimension (tdim);
	- number of interaction configurations (conf_num);
	- number of equilibration sweeps (thermal);
	- mean interaction array (mu);
	- array of the standard deviation of the interactions (sd). 
   The mean interaction array is defined by the minimum value (mu_min), maximum value (mu_max) and the step (mu_step). The size of the mean interaction array
   is given by 1+(mu_max-mu_min)/mu_step.
   The array of the standard deviation of the interactions is define by the minimum value (sd_min), the maximum value (sd_max) and the step (sd_step).
   The size of the array of the standard deviation of the interactions is given by 1+(max_sd-min_sd)/sd_step.

Outputs:
-------
   The program computes: 
	- spin glass susceptibility (Xsg);
	- uniform susceptibility (Xuni);
	- spin glass order parameter (q);
	- magnetization (m);
	- specific heat (c).
   The outputs of the program are .txt files. Each file is the result of the simulation of one pair of (mu, sd), and it contains one line with seven values in this 
   order: mu, sd, Xsg, Xuni, q, m, c.
   The program has a built-in naive parallelization for the mu and sd arrays, i.e., after compiling the code, if the user runs the same program on 10 different instances,
   for example, the total time will decrease 10 times. So the user might be interested in running an array of jobs, each job running the same program for a given set of parameters.
   
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "mt19937ar.h"

// Parameters
#define N 264          // Total number of spins
#define tdim 10000     // Thermal average dimension size
#define conf_num 1000  // Number of interaction configurations
#define thermal 10000  // Equilibration sweeps

#define mu_min -0.002   // Minimum Average Interaction
#define mu_max 0.01     // Maximum Average Interaction
#define mu_step 0.0005  // Average interaction step
#define sd_min 0        // Minimum interaction standard deviation
#define sd_max 0.15     // Maximum interaction standard deviation
#define sd_step 0.0075  // Interaction sd step


// Matrix Serializations
#define SSArr(SS,i,j) SS[ tdim*i + j ] // Spin series array
#define JArr(J,i,j) J[ N*i + j ] // Interactions array
#define CovArr(var,i,j) var[N*i + j] // Covariance array

// Global Variables
int *s; // Spins
int *SS; // Spin series array (N x T)
double *J; // Interactions (Jij) drawn from a gaussian distribution

// Functions
double randnum(); // Real random number [0,1)
void init_randnum(); // Seed of the random number generator
void sweep(); // Sweep of the system with N metropolis steps
double cov(int *A, int k, int l); // Covariance cov[A(k),A(l)]
void elapsed_time();


int main()
{

///////////////////////////// Summary /////////////////////////////
	printf("\n\nSummary:\n\n");
	printf("Number of spins: %d\n", N);
	printf("Thermal average dimension: %d\n",tdim);
	printf("Number of configurations: %d\n",conf_num);
	printf("Equilibration sweeps: %d\n",thermal);
	printf("Number of parameter values: %d mean(J),  %d std(J)", (int)round(((mu_max - mu_min) / mu_step + 1)), (int)round(((sd_max - sd_min) / sd_step + 1)));
	printf("\n\n");

///////////////////////////// Model Parameters /////////////////////////////

	/* Average interaction array (main function variable) */
	int mu_size = (int)round( (mu_max - mu_min) / mu_step + 1); //printf("\n\nmu size: %d\n\n", mu_size);
	double *mu=malloc(mu_size*sizeof(double));
	if(!mu){
		printf("Out of memory.1\n");
		exit(1);
	}
	for(int i=0; i<mu_size; i++){
		mu[i] = (mu_min + i*mu_step); //printf("m[%d]=%f\n\n",i,mu[i]);
	}

	/* Interaction standard deviation array (main function variable) */
	int sd_size = (int)round( (sd_max - sd_min) / sd_step + 1); //printf("\n\nsd size: %d\n\n", sd_size);
	double *sd=malloc(sd_size*sizeof(double));
	if(!sd){
		printf("Out of memory.2\n");
		exit(1);
	}
	for(int i=0; i<sd_size; i++){
		sd[i] = (sd_min + i*sd_step); //printf("sd[%d]=%f\n\n",i,sd[i]);
	}




//////////////////////////////// Variables ///////////////////////////////


	/* Spin array (global) */
	s=malloc(N*sizeof(int));
	if(!s){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Spin series array (global) */
	SS=malloc(N*tdim*sizeof(int));
	if(!SS){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Interaction matrix (global) */
	J=malloc(N*N*sizeof(double));
	if(!J){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Covariance matrix configuration average sum */
	double *SumC=calloc(N*N,sizeof(double));
	if(!SumC){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Covariance matrix configuration average product sum */
	double *ProdC=calloc(N*N,sizeof(double));
	if(!ProdC){
		printf("Out of memory.\n");
		exit(1);
	}



	/* Magnetization */
	double *m=malloc(mu_size*sd_size*sizeof(double));
	if(!m){
		printf("Out of memory.\n");
		exit(1);
	}

	/* SG order parameter */
	double *q=malloc(mu_size*sd_size*sizeof(double));
	if(!q){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Uniform Susceptibiliy */
	double *Xuni=calloc(mu_size*sd_size,sizeof(double));
	if(!Xuni){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Spin Glass Susceptibiliy */
	double *Xsg=calloc(mu_size*sd_size,sizeof(double));
	if(!Xsg){
		printf("Out of memory.\n");
		exit(1);
	}

	/* Specific Heat */
	double *c=calloc(mu_size*sd_size,sizeof(double));
	if(!c){
		printf("Out of memory.\n");
		exit(1);
	}




/////////////////////////////// Program ///////////////////////////////

	/* Initialize mt rng */
	init_randnum();

	/* Interaction array variables */
	const gsl_rng_type * Q;
	gsl_rng * r;
	Q = gsl_rng_mt19937;
	r = gsl_rng_alloc (Q);

	//////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// AVERAGE LOOP (A) /////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	for(int A=0; A<mu_size; A++){

	//////////////////////////////////////////////////////////////////////////////////
	///////////////////////// STANDARD DEVIATION LOOP (B) ////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	for(int B=0; B<sd_size; B++){


		// Output file
		FILE *file;
		char *fname=malloc(100*sizeof(char));
		snprintf(fname,100,"file_%d_%d.txt",A,B);
		
		// Check if file exists
		if( access( fname, F_OK ) != -1 ) {
		    // file exists
		    continue;
		}


		// file doesn't exist: perform calculations
		
		// Save blank file
		if((file = fopen(fname,"w"))==NULL){
			printf("Cannot open file.\n");
		        exit(1);
		}
		fclose(file);
		printf("\nfile_%d_%d.txt  mu=%f  sd=%f",A,B, mu[A], sd[B]);

		//~~~~~~~~~ Beginning of calculations

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Initial condition
		for(int i=0; i<N; i++) s[i]=1;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Reset convariance array 
		memset(SumC, 0, N*N*sizeof(*SumC));
		memset(ProdC, 0, N*N*sizeof(*ProdC));

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Reset specific heat variables
		double ProdSumE = 0;
		double SumProdE = 0;
	
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Configuration Average
		for(int D=0; D<conf_num; D++){
	
	
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;
			// Interaction array
			for(int i=0; i<N; i++){
			for(int j=i; j<N; j++){
				JArr(J,i,j) = mu[A] + gsl_ran_gaussian_ziggurat(r,sd[B]);
				JArr(J,j,i) = JArr(J,i,j);
				JArr(J,i,i) = 0;
			}}
	
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Thermalization
			for(int i=0; i<thermal; i++){
				sweep(); 
			}
				
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Sampling
			for (int j=0; j<tdim; j++){
				sweep();
				for (int i=0; i<N; i++){
					SSArr(SS,i,j)=s[i];
				}
			}
	
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Magnetization and SG order parameter
			double mag=0;
			double sgop=0;
			double aux=0;
			for(int i=0; i<N; i++){
				aux=0;
				for(int t=0; t<tdim; t++){
					aux += SSArr(SS,i,t);
				}
				mag += aux;
				sgop += aux*aux;
			}
			m[sd_size*A + B] += mag;
			q[sd_size*A + B] += sgop;
			
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Covariance Matrix: configuration average sum (SumC) for Xuni and product (ProdC) for Xsg
			double cov_ij;
			for(int i=0; i<N; i++){
			for(int j=i; j<N; j++){
				cov_ij = cov(SS,i,j);
				CovArr(SumC,i,j) += cov_ij;
				CovArr(SumC,j,i) = CovArr(SumC,i,j);
				CovArr(ProdC,i,j) += cov_ij*cov_ij;
				CovArr(ProdC,j,i) = CovArr(ProdC,i,j);
			}}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Specific heat:
			double sumE=0, sumE2=0, E=0;
			for(int t=1; t<tdim; t++){
				E=0;
				for(int i=0; i<N; i++){
				for(int j=i; j<N; j++){
					E -= JArr(J,i,j) * SSArr(SS,i,t) *SSArr(SS,j,t);
				}}
				sumE += E;
				sumE2 += E*E;
			}
			ProdSumE += sumE*sumE;
			SumProdE += sumE2;
				
			
	
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		} // End of Configuration Average
	
		m[sd_size*A+B] = fabs(m[sd_size*A+B]/(1.0*N*tdim*conf_num));
		q[sd_size*A+B] = q[sd_size*A+B]/(1.0*N*tdim*tdim*conf_num);
	
		int NN=N*N;
		for(int i=0; i<NN; i++){
			Xuni[sd_size*A+B] += SumC[i];
			Xsg[sd_size*A+B] += ProdC[i];
		}
		Xuni[sd_size*A+B] = Xuni[sd_size*A+B]/(1.0*N*conf_num);
		Xsg[sd_size*A+B] = Xsg[sd_size*A+B]/(1.0*N*conf_num);

		c[sd_size*A+B] = (SumProdE/(1.0*tdim) - ProdSumE/(1.0*tdim*tdim))/(1.0*N*conf_num);

		

		// Save file with results
		if((file = fopen(fname,"w"))==NULL){
			printf("Cannot open file.\n");
        		exit(1);
        	}

		fprintf(file,"%f\t%f\t%f\t%f\t%f\t%f\t%f", mu[A],
							   sd[B],
							   Xsg[A*sd_size+B],
							   Xuni[A*sd_size+B],
							   q[A*sd_size+B],
							   m[A*sd_size+B],
							   c[A*sd_size+B]);

        	fclose(file);
		free(fname); fname = NULL;





	}/////////////////////////////////////////////////////////////////////////////////
	///////////////////////// END OF STANDARD DEVIATION LOOP (B) /////////////////////
	//////////////////////////////////////////////////////////////////////////////////


	}/////////////////////////////////////////////////////////////////////////////////
	////////////////////////////// END OF  AVERAGE LOOP (A) //////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	




// !!!! FREE POINTERS !!!!
	free(mu); mu=NULL;
	free(sd); sd=NULL;
	free(s); s=NULL;
	free(SS); SS=NULL;
	free(J); J=NULL;
	free(m); m=NULL;
	free(q); q=NULL;
	free(SumC); SumC=NULL;
	free(ProdC); ProdC=NULL;
	free(Xuni); Xuni=NULL;
	free(Xsg); Xsg=NULL;
	free(c); c=NULL;
	gsl_rng_free (r);


//-------- END OF THE PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	elapsed_time();
	return 0;

}




// 	FUNCTIONS	//

/*********************************************************************************************
		                      Metropolis algorithm 
**********************************************************************************************/

void sweep()
{
	int k;
	double Heff, delta;

	for(int i=0; i<N; i++)
	{
		/* Choose a site at random */
		k = N*randnum();

		/* Calculate the effective field */
		Heff=0;
		for(int j=0; j<N; j++) Heff += JArr(J,k,j)*s[j];

		/* Calculate the change in energy */
		delta = s[k]*(Heff);
		
		/* Decide whether to flip the spin */
		if (delta <= 0){
			s[k] = -s[k];
		}
		else if ( randnum() < exp(-2*delta) ){
			s[k] = -s[k];
		}	
	}
}



/*********************************************************************************************
                      Covariance: function to calculate correlations 
**********************************************************************************************/

double cov(int *A, int k, int l)
{
	long long int P=0, Qk=0, Ql=0;
	double R=0;

	for(int t=0; t<tdim; t++){
		P += SSArr(A,k,t) * SSArr(A,l,t);
		Qk += SSArr(A,k,t);
		Ql += SSArr(A,l,t);
	}

	R = P/(1.0*tdim) - Qk * Ql/(1.0*tdim*tdim);

	return R;
}




/*********************************************************************************************
                 Mersenne Twister random number generator: real number in [0,1) 
**********************************************************************************************/

/*Seeding the MT random number generator*/
void init_randnum()
{
	//init_genrand(time(NULL));
	init_genrand(0);
}
/*Generates a real number "rand()" between [0,1)*/
double randnum()
{
	double num;
	num=genrand_real2();
	return num;
}

/*********************************************************************************************
                                       Execution time 
**********************************************************************************************/

void elapsed_time(void)
{
	printf("\n.\n.\n.\nElapsed time: %f min\n\n", clock()/(60.0*CLOCKS_PER_SEC));
}
