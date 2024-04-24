#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


// Parameter setting
#define SEQNUM 131072
#define SEQLENGTH 178
#define MOTIFLENGTH 16
#define NUMSEEDS 3
#define NUMITER 1


float sequences[SEQNUM*SEQLENGTH];
float bestMotifs[SEQNUM*MOTIFLENGTH];
float results[SEQNUM*MOTIFLENGTH];
float motifs[SEQNUM*MOTIFLENGTH];
float updatedMotif[MOTIFLENGTH];
float psaa[MOTIFLENGTH];


// Elapsed time checker
static inline double timeCheckerCPU( void ) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Function for reading benchmark file
void readBenchmarkData( char* filename ) {
	FILE *data = fopen(filename, "rb");
	if( data == NULL ) {
		printf( "File not found: %s\n", filename );
		exit(1);
	}

	for ( int i = 0; i < SEQNUM; i ++ ) {
		for ( int j = 0; j < SEQLENGTH; j ++ ) {
			fread(&sequences, sizeof(float), 1, data);
		}
	}

	fclose(data);
}

// Function for writing result values
void writeResultData( char* filename ) {
	FILE *data = fopen(filename, "wb");

	fwrite(bestMotifs, sizeof(float), SEQNUM*MOTIFLENGTH, data);

	fclose(data);
}

// Manhattan
float manhattan( float *core, float *target) {
	float total = 0.0;
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		total += abs(core[i] - target[i]);
	}
	
	return total;
}

// Get a score of Results
float scoreResults() {
	// Phase 1
	// Get average motif values
	float average[MOTIFLENGTH];
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		float total = 0.0;
		for ( int j = 0; j < SEQNUM; j ++ ) {
			total += results[i + MOTIFLENGTH*j];
		}
		average[i] = total / SEQNUM;
	}

	// Phase 2
	// Compare between each motif and the average motif values
	// Get the score via Manhattan Distance
	float scoreTotal = 0;
	float motif[MOTIFLENGTH];
	for ( int i = 0; i < SEQNUM; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			motif[j] = motifs[i*MOTIFLENGTH + j];
		}
		scoreTotal += manhattan(average, motif);

	}

	return (scoreTotal / SEQNUM);
}

// Get a score of motifs
float scoreMotifs() {
	// Phase 1
	// Get average motif values
	float average[MOTIFLENGTH];
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		float total = 0.0;
		for ( int j = 0; j < SEQNUM; j ++ ) {
			total += motifs[i + MOTIFLENGTH*j];
		}
		average[i] = total / SEQNUM;
	}

	// Phase 2
	// Compare between each motif and the average motif values
	// Get the score via Manhattan Distance
	float scoreTotal = 0;
	float motif[MOTIFLENGTH];
	for ( int i = 0; i < SEQNUM; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			motif[j] = motifs[i*MOTIFLENGTH + j];
		}
		scoreTotal += manhattan(average, motif);

	}

	return (scoreTotal / SEQNUM);
}

// Make position specific average array
void makePSAA( int seqIdx ) {
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		float total = 0.0;
		for ( int j = 0; j < SEQNUM; j ++ ) {
			if ( j != seqIdx ) {
				total += motifs[i + MOTIFLENGTH*j];
			}
		}
		psaa[i] = total / SEQNUM;
	}
}

// Profiling
void profile( float *sequence ) {
	int sizeDists = SEQLENGTH - MOTIFLENGTH + 1;

	// Phase 1
	// Get manhattan distance values
	float sumDists = 0.0;
	float dists[sizeDists];
	float subSeq[MOTIFLENGTH];
	for ( int i = 0; i < sizeDists; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			subSeq[j] = sequence[i+j];
		}
		float dist = manhattan(subSeq, psaa);
		dists[i] = dist;
		sumDists += dist;
	}

	// Phase 2
	// Randomly designate the updated motif
	float partialSum = 0.0;
	float randVal = (float)rand()/RAND_MAX;
	int position = 0;
	for ( int i = 0; i < sizeDists; i ++ ) {
		partialSum += dists[i];
		if ( (partialSum/sumDists) >= randVal ) {
			position = i;
			break;
		}
	}
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		updatedMotif[i] = sequence[position+i];
	}
}

// Gibbs Sampler
void gibbsSampler() {
	// Phase 1
	// Designate the motifs randomly first
	for ( int i = 0; i < SEQNUM; i ++ ) {
		int sp = rand() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			float f = sequences[SEQLENGTH*i + sp + j];
			motifs[MOTIFLENGTH*i + j] = f;
			results[MOTIFLENGTH*i + j] = f;
		}
	}
	float resultsScore = scoreMotifs();

	// Phase 2
	// Update the motifs over #sequence
	float sequence[SEQLENGTH];
	for ( int i = 0; i < NUMITER; i ++ ) {
		for ( int j = 0; j < SEQNUM; j ++ ) {
			// Pick a sequence sequencially
			for ( int k = 0; k < SEQLENGTH; k ++ ) {
				sequence[k] = sequences[SEQLENGTH*j + k];
			}
			
			// Make PSAA first
			makePSAA(j);
			
			// Go to profiling step then
			profile(sequence);
			
			// Get an updated motif of the picked sequence
			for ( int k = 0; k < MOTIFLENGTH; k ++ ) {
				motifs[MOTIFLENGTH*j + k] = updatedMotif[k];
			}

			// Compare the scores & Update motifs or not based on the score
			float currentScore = scoreMotifs();
			if ( currentScore < resultsScore ) {
				memcpy(results, motifs, sizeof(results));
				resultsScore = currentScore;
			}
		}
	}
	
}

// Gibbs Sampler Wrapper
void gibbsSamplerWrapper() {
	float bestScore = 0.0;
	
	for ( int i = 0; i < NUMSEEDS; i ++ ) {
		gibbsSampler();
		if ( i == 0 ) {
			float currentScore = scoreResults();
			bestScore = currentScore;
			memcpy(bestMotifs, results, sizeof(bestMotifs));
		} else {
			float currentScore = scoreResults();
			if ( currentScore < bestScore ) {
				bestScore = currentScore;
				memcpy(bestMotifs, results, sizeof(bestMotifs));
			}
		}
		printf( "NUMSEEDS: %d\n", i );
		fflush( stdout );
	}
}


// Main
int main() {
	srand(time(NULL));
	//--------------------------------------------------------------------------------------------
	// Read the sequence fasta format file first
	//--------------------------------------------------------------------------------------------
	char benchmark_filename[] = "../../../dataset/sleep.bin";
	readBenchmarkData(benchmark_filename);
	printf( "Reading the sequence file finished!\n" );
	fflush( stdout );
	//--------------------------------------------------------------------------------------------
	// Run GibbsSampler
	//--------------------------------------------------------------------------------------------
	printf( "Motif Finder Started!\n" );
	fflush( stdout );
	double processStart = timeCheckerCPU();
	gibbsSamplerWrapper();
	double processFinish = timeCheckerCPU();
	printf( "Motif Finder Finished!\n" );
	fflush( stdout );
	double processTime = processFinish - processStart;
	//--------------------------------------------------------------------------------------------
	// Print elapsed time
	//--------------------------------------------------------------------------------------------
	printf( "Elapsed Time: %.8f\n", processTime );
	fflush( stdout );
	//--------------------------------------------------------------------------------------------
	// Write best motifs to bin file
	//--------------------------------------------------------------------------------------------
	char result_filename[] = "../../../dataset/sleep_result.bin";
	writeResultData(result_filename);

	return 0;
}
