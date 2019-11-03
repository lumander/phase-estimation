#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#define SUBS 10000

void discretize ( double *, const int, const int, const double );
void initialize ( double *, const int, const double );
void hoodUpdate ( const int, double *, const int, const double *, const double );
void readPrior ( double *, const int, const int );
void priorUpdAndNorm ( const int, double *, double * );
void cumulUpdate ( const int, double *, const double * );
double randomGen ();
int binarySearch ( double *, double, int, int );
double phaseToVoltage ( double );
void writePrior ( double *, const int );


int main (){
  
	int index, i, click, N;
	double rndm = 0, delta = (2 * M_PI)/SUBS;
	double x, randPrior;
	double middle[SUBS];
	double prior[SUBS];
	double cumPrior[SUBS];
	double likelihood[SUBS];						// CONSTANTS AND VECTORS
    
    	
   		discretize ( middle, SUBS, 0, delta ); 		// 0;2Pi DISCRETIZATION
			
                initialize ( likelihood, SUBS, 0.0 );   		//LIKELIHOOD VECTOR INITIALIZATION
		
		initialize ( cumPrior, SUBS, 0.0 );			//CUMULATIVE VECTOR INITIALIZATION
 					
		hoodUpdate( SUBS, likelihood, click, middle, randPrior ); 	//LIKELIHOOD UPDATING
 		
		readPrior ( prior, SUBS, N );      							//READING INITIAL PRIOR
 		
		priorUpdAndNorm ( SUBS, prior, likelihood );					//PRIOR UPDATING WITH BAYES
 		 				
		cumulUpdate ( SUBS, cumPrior, prior ); 						//UPDATING CUMULATIVE PRIOR
		
		rndm = randomGen();											//RANDOM NUMBER GENERATOR				
	
		index = binarySearch( cumPrior, rndm, 0, SUBS-1 );			//INDEX OF FEEDBACK PHASE RANDOM DISTRIBUTED FROM PRIOR WITH BINARY SEARCH	 
	
		randPrior = phaseToVoltage( *( middle + index ) );				//FEEDBACK PHASE CONVERTED IN VOLTAGE
		
		writePrior ( prior, SUBS );									//PRIOR UPDATING ON FILE
		
 
	return 0;

};
    

void discretize ( double *vec, const int length, const int ex , const double del ) {
	
	int i;
	
		for ( i = 0; i < length; i++ ) {
			
             	*( vec + i ) = ex + (2*i+1)*del/2; 
             	
             		}; 	
             		
};

	
void initialize ( double *vec, const int length, const double initValue ) {

	int i;
	
		for ( i = 0; i < length; i++ ) {
			
          		*( vec + i ) = initValue;
          		
				  };
				  
};


void hoodUpdate( const int subs, double *likehood, const int detect, const double *mid, const double picked ) {
	
	int i;
		
		for ( i = 0; i < subs; i++ ) {
			
            *( likehood + i ) = detect * pow (sin(( *( mid + i )-picked)/2), 2) + 
            
				(1 - detect) * pow (cos(( *( mid + i )-picked)/2), 2); 
           
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
  	
  	
double randomGen(){

	double nmbr;
	int seed;
		
		seed = time(0);	
		
			srand(seed);	
								
				nmbr = (double) rand () / RAND_MAX; 
		
		
		return nmbr;
};	
	  	
	  	
int binarySearch ( double *vec, double rndm, const int first, const int last ) {
    	
	int mid,low,high;  
    	low = first;
   		high = last;
	  	
     	while ( low <= high ) {
     		
			mid = rint ( (low + high) / 2 ) ;
			 	
				if ( *( vec + mid-1 ) <= rndm && rndm <= *( vec + mid+1 ) ) {
					
			 		if ( *( vec + mid ) <= rndm )
			 		
			 			{return mid = mid + 1;} 
			 			
						 else {
						 	
						 return mid;}
						 
					 } else {
					 	
					 	if ( ( rndm - *( vec + mid-1 ) ) > 0) {
					 		
						 	low = mid;} else {
						 		
						 					high = mid;}
					 	
				};
					
		};
	
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
	
	
	
