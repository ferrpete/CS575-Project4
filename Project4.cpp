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

unsigned int seed1;
unsigned int seed2;

int NowYear;      // 2014-2019
int NowMonth;     // 0 - 11

float NowPrecip;  // inches of rain per month
float NowTemp;    // temperature this month
float NowHeight;  // grain height in inches
int NowNumDeer;   // current deer population

float Ranf( float low, float high, unsigned int* seed );
void Temperature( );
void Precipitation( );
float ConditionFactor( );
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
	NowYear = 2014;
	
	Temperature( );
	Precipitation( );

	omp_set_num_threads(3);
	
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

float ConditionFactor( )
{
	
	float TF = exp( -pow( ( NowTemp - MIDTEMP ) / 10., 2 ) );
	float PF = exp( -pow( ( NowPrecip - MIDPRECIP ) / 10., 2 ) );
	
	float CF = TF * PF;
	
	return CF;
	
}

void GrainDeer( )
{
	
	while( NowYear < 2020 )
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
	
	while( NowYear < 2020 )
	{
	
		float CF = ConditionFactor( );
		
		float TmpNowHeight = NowHeight + CF * GRAIN_GROWS_PER_MONTH;
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
	
	while( NowYear < 2020 )
	{
	
		#pragma omp barrier
		
		#pragma omp barrier
		
		std::ofstream dataFile("GrainDeerData.txt", std::ofstream::out | std::ofstream::app);
		dataFile << NowMonth << ", " << NowYear << ", " << NowTemp << ", " << NowPrecip << ", " << NowHeight << ", " << NowNumDeer << "\n";
		
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
