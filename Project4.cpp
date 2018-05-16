#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <omp.h>

const float GRAIN_GROWS_PER_MONTH =        6.0;
const float ONE_DEER_EATS_PER_MONTH =      0.1;
const float GRAIN_THRESHOLD =              120.0;

const float AVG_PRECIP_PER_MONTH =         6.0;
const float AMP_PRECIP_PER_MONTH =         6.0;
const float RANDOM_PRECIP =                2.0;

const float AVG_TEMP =                     50.0;
const float AMP_TEMP =                     20.0;
const float RANDOM_TEMP =                  10.0;

const float MIDTEMP =                      40.0;
const float MIDPRECIP =                    10.0;

const int STARTYEAR =                      1970;
const int ENDYEAR =                        2020;

unsigned int seed1;
unsigned int seed2;

int NowYear;      // 2014-2019
int NowMonth;     // 0 - 11

float NowPrecip;  // inches of rain per month
float NowTemp;    // temperature this month
float NowHeight;  // grain height in inches
int NowNeurotoxin; // deadly neurotoxin
int NowNumDeer;   // current deer population
float GF;         // grain growth factor
float TrainingData[ ( ENDYEAR - STARTYEAR ) * 12 ][8] = {0.}; // training data array
int TruthData[ ( ENDYEAR - STARTYEAR ) * 12 ]; // training data truth labels
float Weights[8]; // basic perceptron weights

float Ranf( float low, float high, unsigned int* seed );
void Temperature( );
void Precipitation( );
float FirstDerivative( float NewHeight, float PrevHeight );
float SecondDerivative( float NewHeight, float PrevHeight, float PrevPrevHeight );
float GrowthFactor( );
void GrainDeer( );
void Grain( );
void Watcher( );
void GLADoS( );


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
	
	int N = sizeof(Weights) / sizeof( Weights[0] );
	for( int i = 0; i < N; i++ )
		printf( "%8.2lf \n", Weights[i] );
	
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

float FirstDerivative( float NewHeight, float PrevHeight )
{
	
	float dHeight = NewHeight - PrevHeight;
	
	return dHeight;
	
}

float SecondDerivative( float NewHeight, float PrevHeight, float PrevPrevHeight )
{
	
	float d2Height = PrevPrevHeight - 2 * PrevHeight + NewHeight;
	
	return d2Height;
	
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
			TmpNumberDeer = NowNumDeer + 2;
		
		#pragma omp barrier
		
		NowNumDeer = TmpNumberDeer;
		
		#pragma omp barrier
		
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
		
		#pragma omp barrier
	
	}
	
}

void Watcher( )
{
	
	while( NowYear < ENDYEAR )
	{
	
		#pragma omp barrier
		
		#pragma omp barrier
		
		#pragma omp barrier
		
		std::ofstream dataFile("GrainDeerData.txt", std::ofstream::out | std::ofstream::app );
		dataFile << NowMonth << ", " << NowYear << ", " << (5. / 9.)*( NowTemp - 32. ) << ", " << NowPrecip * 2.54 << ", " << NowHeight * 2.54 << ", " << NowNumDeer << ", " << NowNeurotoxin << "\n";
		
		int TrainingIdx = ( NowYear - STARTYEAR ) * 12 + NowMonth;
		TrainingData[ TrainingIdx ][0] = NowTemp;
		TrainingData[ TrainingIdx ][1] = NowPrecip;
		TrainingData[ TrainingIdx ][2] = NowHeight;
		TrainingData[ TrainingIdx ][3] = NowNumDeer;
		
		float dHeight = FirstDerivative( TrainingData[ TrainingIdx ][2], TrainingData[ TrainingIdx - 1 ][2] );
		float d2Height = SecondDerivative( TrainingData[ TrainingIdx ][2], TrainingData[ TrainingIdx - 1 ][2], TrainingData[ TrainingIdx - 2 ][2] );
		float dDeer = FirstDerivative( TrainingData[ TrainingIdx ][3], TrainingData[ TrainingIdx - 1 ][3] );
		
		TrainingData[ TrainingIdx ][4] = dHeight;
		TrainingData[ TrainingIdx ][5] = d2Height;
		TrainingData[ TrainingIdx ][6] = dDeer;
		TrainingData[ TrainingIdx ][7] = 1;
		
		if( ( NowHeight >= GRAIN_THRESHOLD ) )
			TruthData[ TrainingIdx ] = -1;
		else
			TruthData[ TrainingIdx ] = 1;
		
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
		
		if( NowYear > STARTYEAR + ( ENDYEAR - STARTYEAR ) / 10 )
		{
			
			int TotalEpochs = 10;
			int N = ( NowYear - STARTYEAR ) * 12 + NowMonth;
			int M = sizeof(TrainingData[0]) / sizeof(TrainingData[0][0]);
			
			for( int kk = 0; kk < M; kk++ )
				Weights[kk] = 0.;
			
			for( int j = 0; j < TotalEpochs; j++ )
			{
				
				for( int i = N - 144; i < N; i++ )
				{
					float dotProduct = 0;
					float truth = TruthData[i];
					
					for( int k = 0; k < M; k++ )
					{
						
						dotProduct += TrainingData[i][k] * Weights[k];
						
					}
					
					//printf( "Dot Product: %8.2f \n", dotProduct );
						
					if( truth * dotProduct <= 0 )
					{
						
						for( int ii = 0; ii < M; ii++ )
							Weights[ii] += truth * TrainingData[i][ii];
						
					}
					
				}
				
			}
			
		}
			
		#pragma omp barrier
		
		#pragma omp barrier
		
		if( NowYear > STARTYEAR + ( ENDYEAR - STARTYEAR ) / 10 )
		{
		
			float Decision = 0.;
			int N = ( NowYear - STARTYEAR ) * 12 + NowMonth - 1;
			int M = sizeof(Weights) / sizeof(Weights[0]);
			
			for( int i = 0; i < M; i++ )
			{
				Decision += (float)TrainingData[N][i] * Weights[i];
				
			}
			
			//printf( "Decision: %8.2lf \n", Decision );
			
			int K = 1000;
			
			if( Decision > 0 )
			{
				NowNeurotoxin = (int) K * abs( GRAIN_THRESHOLD - NowHeight ) * ONE_DEER_EATS_PER_MONTH;
				NowNumDeer -= NowNeurotoxin;
			}
			else
			{
				NowNeurotoxin = 0;
			}
			
			if( NowNumDeer < 0 )
				NowNumDeer = 0;
			
		}
		else
		{
			NowNeurotoxin = 0;
		}
		
		#pragma omp barrier
		
		#pragma omp barrier
		
	}
	
}
