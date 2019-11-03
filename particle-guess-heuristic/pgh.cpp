#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#define SUBS 10000

void discretize ( double *, const int, const int, const double );
void hoodUpdate ( const int, double *, const int, const double *, const int );
void readPrior ( double *, const int, const int );
void priorUpdAndNorm ( const int, double *, double * );
void cumulUpdate ( const int, double *, const double * );
double randomGen ( const int );
int binarySearch ( double *, double, int, int, const double );
void randPriorPhaseRecord ( double );
double phaseToVoltage ( double );
void writePrior ( double *, const int );


int main ( int argc, const char *argv[] ) {
  
	int index, i, click, N, iter;
	double rndm = 0, delta = (2 * M_PI)/SUBS;
	double x, randPriorTension, *middle, *cumPrior, *likelihood, *prior;
	
	middle = ( double * ) malloc ( SUBS * sizeof(double ) );
	prior = ( double * ) malloc( SUBS * sizeof( double ) );
	cumPrior = ( double * ) calloc( SUBS, sizeof( double ) );
	likelihood = ( double * ) calloc( SUBS, sizeof( double ) );						// CONSTANTS AND VECTORS
    
    click = atoi ( *( argv + 1 ) );	
    
    N = atoi ( *( argv + 2 ) ); 
    
    iter = atoi ( *( argv + 3 ) );
    
    	discretize ( middle , SUBS, 0, delta ); 				//DISCRETIZING MIDDLE
 					
		hoodUpdate( SUBS, likelihood, click, middle, N ); 	//LIKELIHOOD UPDATING
 		
		readPrior ( prior, SUBS, N );      							//READING INITIAL PRIOR
 		
		priorUpdAndNorm ( SUBS, prior, likelihood );					//PRIOR UPDATING WITH BAYES
 		
			free(likelihood);
		 				
		cumulUpdate ( SUBS, cumPrior, prior ); 						//UPDATING CUMULATIVE PRIOR
		
		rndm = randomGen( iter );											//RANDOM NUMBER GENERATOR				
	
		index = binarySearch( cumPrior, rndm, 0, SUBS-1, delta );			//INDEX OF FEEDBACK PHASE RANDOM DISTRIBUTED FROM PRIOR WITH BINARY SEARCH	 
		
		randPriorPhaseRecord ( *( middle + index ) );
		
			free(cumPrior);
			
		randPriorTension = phaseToVoltage( *( middle + index ) );				//FEEDBACK PHASE CONVERTED IN VOLTAGE
		
		writePrior ( prior, SUBS );									//PRIOR UPDATING ON FILE
		
		free(middle);
		free(prior);
 		free(cumPrior);
 		free(likelihood);
 		
		printf("%f", randPriorTension);
 
	return 0;

};
    

void discretize ( double *vec, const int length, const int ex , const double del ) {
	
	int i;
	
		for ( i = 0; i < length; i++ ) {
			
             	*( vec + i ) = ex + (2*i+1)*del/2; 
             	
             		}; 	
             		
};


void hoodUpdate( const int subs, double *likehood, const int detect, const double *mid, const int photon ) {
	
	int i;
	double picked,x;
	FILE *fp;	

	if ( photon == 1 ) { 	
		
		picked = 0;
		
		for ( i = 0; i < subs; i++ ) {
			
            *( likehood + i ) = detect * pow (sin(( *( mid + i )-picked)/2), 2) + 
            
				(1 - detect) * pow (cos(( *( mid + i )-picked)/2), 2); 
           
		    };
		
		} else{
		    	
			 if ( ( fp = fopen ("randPrior.dat", "r") ) == NULL ) {
				
					   	printf ("Error when opening the file initprior.dat\n");
					   	
					  		exit( EXIT_FAILURE );};
	   						
		 						for (i = 0; i < 1; i++) {
		 							
									fscanf(fp , "%lf", &x);
									
										picked = x;};
										
					fclose (fp);
		
		for ( i = 0; i < subs; i++ ) {
			
            *( likehood + i ) = detect * pow (sin(( *( mid + i )-picked)/2), 2) + 
            
				(1 - detect) * pow (cos(( *( mid + i )-picked)/2), 2); 
           
		    };  
			
	};
			
	
		    
		    
};


void readPrior ( double *vec, const int length, const int photon ) {

	FILE *fp;
	int i;
	double x;
	
		if( photon==1 ) {	
				
			if ( ( fp = fopen ("initprior.dat", "r") ) == NULL ) {
				
					   	printf ("Error when opening the file initprior.dat\n");
					   	
					  		exit( EXIT_FAILURE );};
	   						
		 						for (i = 0; i < SUBS; i++) {
		 							
									fscanf(fp , "%lf", &x);
									
										*( vec + i ) = x;}; 	
 
			} else {
			
				if ( ( fp = fopen ("newprior.dat", "r") ) == NULL ) {
					
	   						printf ("Error when opening the file newprior.dat\n");
	   						
	   							exit( EXIT_FAILURE );};
	   
	   								for (i = 0; i < SUBS; i++) {
	   									
										fscanf(fp , "%lf", &x);
										
											*( vec + i ) = x;};
	   
	   };
	   
	fclose(fp);

};
 	 	
 	
void priorUpdAndNorm ( const int subs, double *vec, double *likehood ) {
	
	int i;
	double norm=0.0;
	
			for ( i = 0; i < subs; i++ ) {
				
            	*( vec + i ) *= *( likehood + i ); 
            	
            		};  
 	 
    			for ( i = 0; i < subs; i++ ) {
    				
            		norm += *( vec + i ); 
            		
                		};   
    
	    				for ( i = 0; i < subs; i++ ) {
	    					
            				*( vec + i ) /= norm; 
            				
             					};          

};	       

		        
void cumulUpdate ( const int subs, double *cumVec, const double *vec ) {
	
	int i,j;
	
			for ( i = 0; i < subs; i++ ) {
				
        		for ( j = 0; j <= i; j++ ) {
        			
					*( cumVec + i ) += *( vec + j );
					
				  			};
				   }; 
				      
};
  	
  	
double randomGen( const int k ){

		
								
	int i;
	double x,nmbr;
	FILE *fp;

							if ( ( fp = fopen ("rndms.dat", "r+") ) == NULL ) {
	
	   									printf ("Error when opening the file newprior.dat\n");
	   
	  										exit( EXIT_FAILURE );}
     
     								for (i = 0; i <= k; i++) {
     	
									fscanf(fp , "%lf\n",  x ); 
									
									nmbr = x;
									
								};

 														
 
	fclose(fp);
		
		
		return nmbr;
};	
	  	
	  	
int binarySearch ( double *vec, double rndm, int first, int last , const double del) {
    	
	int mid;
	  	
     	while ( first <= last ) {
     		
			mid = rint ( (first + last) / 2 ) ;
			 	
				if ( fabs( *( vec + mid ) - rndm )  <= del/2   ) {
					
						return mid;} 
						
			 		 else {
					 	
					 	if ( ( rndm - *( vec + mid) ) > 0) {
					 		
						 	first = mid+1;} else {
						 		
						 					last = mid-1;}
					 	
				};
					
		};
	
};


void randPriorPhaseRecord ( double randPriorPhase ) {

	FILE *fp;
	int i;

							if ( ( fp = fopen ("randPrior.dat", "w") ) == NULL ) {
	
	   									printf ("Error when opening the file newprior.dat\n");
	   
	  										exit( EXIT_FAILURE );}
     
     							for (i = 0; i < 1; i++) {
     	
									fprintf(fp , "%lf", randPriorPhase);

 														};
 
	fclose(fp);


};	  	
	  	
	  	
double phaseToVoltage ( double phase ){

FILE *fp;
int i, lengthVec = 7299;
double k,x,min = 0, tension, value = 2*M_PI;


					if ( ( fp = fopen ("tensionphases.dat", "r") ) == NULL ) {
						
 								printf ("Error when opening the file tensionphases.dat\n");
 								
 									exit( EXIT_FAILURE );};

							for (i = 0; i <= lengthVec; i++) {
									
 								fscanf(fp , "%lf\t%lf", &k, &x);
 									
								min = fabs(phase-x);
									
										if( min < value ){
												
											value = min;
												
											tension = k;}; 
												
											};
												
	fclose(fp);

		return tension;
};	
	
	
void writePrior( double *vec, const int length ){

	int i;
	FILE *fp;

							if ( ( fp = fopen ("newprior.dat", "w") ) == NULL ) {
	
	   									printf ("Error when opening the file newprior.dat\n");
	   
	  										exit( EXIT_FAILURE );}
     
     							for (i = 0; i < length; i++) {
     	
									fprintf(fp , "%lf", *(vec + i));

 														};
 
	fclose(fp);
     		
};
	
	
	
