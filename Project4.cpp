#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <omp.h>

const float GRAIN_GROWS_PER_MONTH =        8.0;
const float ONE_DEER_EATS_PER_MONTH =      0.5;

const float AVG_PRECIP_PER_MONTH =         6.0;
const float AMP_PRECIP_PER_MONTH =         6.0;
const float RANDOM_PRECIP =                2.0;

const float AVG_TEMP =                     50.0;
const float AMP_TEMP =                     20.0;
const float RANDOM_TEMP =                  10.0;

const float MIDTEMP =                      40.0;
const float MIDPRECIP =                    10.0;

const int STARTYEAR =                      2014;
const int ENDYEAR =                        2020;

unsigned int seed1;
unsigned int seed2;

int NowYear;      // 2014-2019
int NowMonth;     // 0 - 11

float NowPrecip;  // inches of rain per month
float NowTemp;    // temperature this month
float NowNeuroToxin; // deadly neurotoxin
float NowHeight;  // grain height in inches
int NowNumDeer;   // current deer population
float GF;         // grain growth factor
float TrainingData[ ( ENDYEAR - STARTYEAR ) * 12 ][6]; // training data array
float Weights[5]; // basic perceptron weights

float Ranf( float low, float high, unsigned int* seed );
void Temperature( );
void Precipitation( );
float GrowthFactor( );
void GrainDeer( );
void Grain( );
void Watcher( );


int 
main( )
{
	
#ifndef _OPENMP
        fprintf( stderr, "OpenMP is not supported here -- sorry.\n" );
        return 1;
#endif

	unsigned int seed1 = omp_get_wtime( );
	unsigned int seed2 = omp_get_wtime( );

	NowNumDeer = 1;
	NowHeight = 1.;
	NowMonth = 0;
	NowYear = STARTYEAR;
	
	Temperature( );
	Precipitation( );

	omp_set_num_threads(4);
	
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			GrainDeer();
		}
		#pragma omp section
		{
			Grain();
		}
		#pragma omp section
		{
			Watcher();
		}
		#pragma omp section
		{
			GLADoS();
		}
		// implied barrier: all sections must complete before we get here
	}
	
	return 0;

}

float Ranf( float low, float high, unsigned int* seed )
{
	
	float r = (float)rand_r(seed);    // 0 - RAND_MAX
	return( low + r * ( high - low ) / (float)RAND_MAX );
	
}

void Temperature( )
{
	
	float ang = ( 30.*(float)NowMonth + 15. ) * ( M_PI / 180. );
	float temp = AVG_TEMP - AMP_TEMP * cos( ang );
	
	NowTemp = temp + Ranf( -RANDOM_TEMP, RANDOM_TEMP, &seed1 );
	
}

void Precipitation( )
{
	
	float ang = ( 30.*(float)NowMonth + 15. ) * ( M_PI / 180. );
	float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin( ang );
	
	NowPrecip = precip + Ranf( -RANDOM_PRECIP, RANDOM_PRECIP, &seed2 );
	
	if( NowPrecip < 0. )
		NowPrecip = 0.;
	
}

float GrowthFactor( )
{
	
	float TF = exp( -pow( ( NowTemp - MIDTEMP ) / 10., 2 ) );
	float PF = exp( -pow( ( NowPrecip - MIDPRECIP ) / 10., 2 ) );
	
	float GF = TF * PF;
	
	return GF;
	
}

void GrainDeer( )
{
	
	while( NowYear < ENDYEAR )
	{
	
		int TmpNumberDeer;
		
		if( NowHeight < (float) NowNumDeer )
			TmpNumberDeer = NowNumDeer - 1;
		else
			TmpNumberDeer = NowNumDeer + 1;
		
		#pragma omp barrier
		
		NowNumDeer = TmpNumberDeer;
		
		#pragma omp barrier
		
		#pragma omp barrier
	
	}
	
}

void Grain( )
{
	
	while( NowYear < ENDYEAR )
	{
	
		GF = GrowthFactor( );
		
		float TmpNowHeight = NowHeight + GF * GRAIN_GROWS_PER_MONTH;
		TmpNowHeight -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
		
		if( TmpNowHeight < 0 )
			TmpNowHeight = 0;
		
		#pragma omp barrier
		
		NowHeight = TmpNowHeight;
		
		#pragma omp barrier
		
		#pragma omp barrier
	
	}
	
}

void Watcher( )
{
	
	while( NowYear < ENDYEAR )
	{
	
		#pragma omp barrier
		
		#pragma omp barrier
		
		std::ofstream dataFile("GrainDeerData.txt", std::ofstream::out | std::ofstream::app);
		dataFile << NowMonth << ", " << NowYear << ", " << (5. / 9.)*( NowTemp - 32. ) << ", " << NowPrecip * 2.54 << ", " << NowHeight * 2.54 << ", " << NowNumDeer << "\n";
		
		TrainingIdx = ( NowYear - STARTYEAR ) * 12 + NowMonth;
		TrainingData[ TrainingIdx ][1] = NowTemp;
		TrainingData[ TrainingIdx ][2] = NowPrecip;
		TrainingData[ TrainingIdx ][3] = NowHeight;
		TrainingData[ TrainingIdx ][4] = NowNumDeer;
		
		if( NowYear >= ( ENDYEAR - STARTYEAR ) / 2 )
			TrainingData[ TrainingIdx ][5] = NowNeuroToxin;
		if( ( GC >= 1 ) && ( NowNumDeer < NowHeight ) )
			TrainingData[ TrainingIdx ][6] = 1;
		else
			TrainingData[ TrainingIdx ][6] = -1;
		
		NowMonth += 1;
		
		if( NowMonth == 12 )
		{
			
			NowMonth = 0;
			NowYear += 1;
			
		}
		
		Temperature( );
		Precipitation( );
		
		#pragma omp barrier
	
	}
	
}

void GLADoS( )
{
	
	while( NowYear < ENDYEAR )
	{
		
		if( NowYear >= ( ENDYEAR - STARTYEAR ) / 2 )
		{
			
			int TotalEpochs = 5;
			int N = sizeof(TrainingData)/sizeof(TrainingData[0]);
			int M = sizeof(TrainingData[0]) / sizeof(TrainingData[0][0]);
			
			for( j = 0; j < TotalEpochs; j++ )
			{
				
				for( i = 0; i < N; i++ )
				{
					float yi = 0;
					
					for( k = 0; k < M; k++ )
					{
						
						yi += TrainingData[i][k] * Weights[k];
						
						if( yi <= 0 )
						{
							
							for( ii = 0; ii < M; ii++ )
								Weights[ii] += yi * TrainingData[i][ii];
							
						}
						
					}
					
				}
				
			}
			
			#pragma omp barrier
			
			int M = sizeof(TrainingData[0]) / sizeof(TrainingData[0][0]);
			
			for( i = 0; i < M; i++ )
			{
				
				
				
			}
			
			#pragma omp barrier
			
			 
			
		}
		
	}
	
}
