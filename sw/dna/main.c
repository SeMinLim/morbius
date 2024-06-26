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
#define SEQLENGTH 1000
#define MOTIFLENGTH 16
#define NUMSEEDS 3
#define NUMITER 1


int base[4*MOTIFLENGTH];
int resultsScore = 0;
//int positionQ[SEQNUM*SEQLENGTH*NUMSEEDS];
//int posIdx = 0;


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Get a base matrix
void getBase( char *motifs ) {
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		int a = 1, c = 1, g = 1, t = 1;
		for ( int j = 1; j < SEQNUM; j ++ ) {
			if ( motifs[MOTIFLENGTH*j+i] == 'A' ) a += 1;
			else if ( motifs[MOTIFLENGTH*j+i] == 'C' ) c += 1;
			else if ( motifs[MOTIFLENGTH*j+i] == 'G' ) g += 1;
			else if ( motifs[MOTIFLENGTH*j+i] == 'T' ) t += 1;
		}
		base[i] = a;
		base[MOTIFLENGTH*1 + i] = c;
		base[MOTIFLENGTH*2 + i] = g;
		base[MOTIFLENGTH*3 + i] = t;
	}
}

// Get a score of motifs
int score( char *motifs, char *updatedMotif ) {
	// Phase 1
	// Pick the most frequent letter at each motif position
	char pattern[MOTIFLENGTH];
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		if ( updatedMotif[i] == 'A' ) base[i] += 1;
		else if ( updatedMotif[i] == 'C' ) base[MOTIFLENGTH*1 + i] += 1;
		else if ( updatedMotif[i] == 'G' ) base[MOTIFLENGTH*2 + i] += 1;
		else if ( updatedMotif[i] == 'T' ) base[MOTIFLENGTH*3 + i] += 1;

		int a = base[i];
		int c = base[MOTIFLENGTH*1 + i];
		int g = base[MOTIFLENGTH*2 + i];
		int t = base[MOTIFLENGTH*3 + i];

		if ( a >= c && a >= g && a >= t ) pattern[i] = 'A';
		else if ( c >= g && c >= t ) pattern[i] = 'C';
		else if ( g >= t ) pattern[i] = 'G';
		else pattern[i] = 'T';
	}

	// Phase 2
	// Compare between each motif and the picked string 
	// Get the score via Hamming Distance 
	int score = 0;
	char motif[MOTIFLENGTH];
	for ( int i = 0; i < SEQNUM; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			motif[j] = motifs[i*MOTIFLENGTH+j];
		}
	
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			if ( motif[j] != pattern[j] ) score += 1;
		}
	}

	return score;
}

// Make position specific score matrix
void makePSSM( float *pssm, char *motifs, int seqIdx ) {
	if ( seqIdx > 0 ) {
		for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
			if ( motifs[MOTIFLENGTH*seqIdx + i] == 'A' ) base[i] -= 1;
			else if ( motifs[MOTIFLENGTH*seqIdx + i] == 'C' ) base[MOTIFLENGTH*1 + i] -= 1;
			else if ( motifs[MOTIFLENGTH*seqIdx + i] == 'G' ) base[MOTIFLENGTH*2 + i] -= 1;
			else if ( motifs[MOTIFLENGTH*seqIdx + i] == 'T' ) base[MOTIFLENGTH*3 + i] -= 1;
		}
	}

	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		int a = base[i];
		int c = base[MOTIFLENGTH*1 + i];
		int g = base[MOTIFLENGTH*2 + i];
		int t = base[MOTIFLENGTH*3 + i];

		pssm[i] = (float)a / (float)SEQNUM;
		pssm[MOTIFLENGTH*1+i] = (float)c / (float)SEQNUM;
		pssm[MOTIFLENGTH*2+i] = (float)g / (float)SEQNUM;
		pssm[MOTIFLENGTH*3+i] = (float)t / (float)SEQNUM;
	}
}

// Profiling
void profile( char *updatedMotif, char *sequence, float *pssm ) {
	int sizeProbs = SEQLENGTH - MOTIFLENGTH + 1;

	// Phase 1
	// Get probabilities
	float sumProbs = 0.0;
	float probs[sizeProbs];
	for ( int i = 0; i < sizeProbs; i ++ ) {
		float prob = 1.0;
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			if ( sequence[i+j] == 'A' ) prob *= pssm[j];
			else if ( sequence[i+j] == 'C' ) prob *= pssm[MOTIFLENGTH*1+j];
			else if ( sequence[i+j] == 'G' ) prob *= pssm[MOTIFLENGTH*2+j];
			else if ( sequence[i+j] == 'T' ) prob *= pssm[MOTIFLENGTH*3+j];
		}

		probs[i] = prob;
		sumProbs += prob;
	}

	// Phase 2
	// Randomly designate the updated motif
	float partialSum = 0.0;
	float randVal = (float)rand()/RAND_MAX;
	int position = 0;
	for ( int i = 0; i < sizeProbs; i ++ ) {
		partialSum += probs[i];
		if ( (partialSum/sumProbs) >= randVal ) {
			position = i;
			//positionQ[posIdx++] = i;		
			break;
		}
	}
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		updatedMotif[i] = sequence[position+i];
	}
}

// Gibbs Sampler
void gibbsSampler( char *results, char *sequences ) {
	// Phase 1
	// Designate the motifs randomly first
	char motifs[SEQNUM*MOTIFLENGTH];
	for ( int i = 0; i < SEQNUM; i ++ ) {
		int sp = rand() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			char c = sequences[SEQLENGTH*i + sp + j];
			motifs[MOTIFLENGTH*i + j] = c;
			results[MOTIFLENGTH*i + j] = c;
		}
	}
	// Get a base matrix
	getBase(motifs);

	// Phase 2
	// Update the motifs over 20 x #sequence
	char sequence[SEQLENGTH];
	char updatedMotif[MOTIFLENGTH];
	float pssm[4*MOTIFLENGTH];
	for ( int i = 0; i < NUMITER; i ++ ) {
		for ( int j = 0; j < SEQNUM; j ++ ) {
			// Pick a sequence sequencially
			for ( int k = 0; k < SEQLENGTH; k ++ ) {
				sequence[k] = sequences[SEQLENGTH*j + k];
			}
			
			// Make PSSM
			makePSSM(pssm, motifs, j);
			
			// Go to profiling step then
			profile(updatedMotif, sequence, pssm);
			
			// Get an updated motif of the picked sequence
			for ( int k = 0; k < MOTIFLENGTH; k ++ ) {
				motifs[MOTIFLENGTH*j + k] = updatedMotif[k];
			}
			
			// Compare the scores & Update motifs or not based on the score
			int currentScore = score(motifs, updatedMotif);
			if ( j == 0 ) {
				for ( int k = 0; k < SEQNUM*MOTIFLENGTH; k ++ ) {
					results[k] = motifs[k];
				}
				resultsScore = currentScore;
			} else {
				if ( currentScore < resultsScore ) {
					for ( int k = 0; k < SEQNUM*MOTIFLENGTH; k ++ ) {
						results[k] = motifs[k];
					}
					resultsScore = currentScore;
				}
			}
		}
	}
}

// Gibbs Sampler Wrapper
void gibbsSamplerWrapper( char *bestMotifs, char *sequences ) {
	int bestScore = 0;
	
	char results[SEQNUM*MOTIFLENGTH + 1];
	for ( int i = 0; i < NUMSEEDS; i ++ ) {
		gibbsSampler(results, sequences);
		if ( i == 0 ) {
			bestScore = resultsScore;
			for ( int j = 0; j < SEQNUM*MOTIFLENGTH; j ++ ) {
				bestMotifs[j] = results[j];
			}
		} else {
			if ( resultsScore < bestScore ) {
				bestScore = resultsScore;
				for ( int j = 0; j < SEQNUM*MOTIFLENGTH; j ++ ) {
					bestMotifs[j] = results[j];
				}
			}
		}
	}
}


// Main
int main() {
	srand(time(NULL));
	//--------------------------------------------------------------------------------------------
	// Read the sequence fasta format file first
	//--------------------------------------------------------------------------------------------
	int seqIdx = 0;
	char seqLine[1024];
	char *sequences = (char*)malloc(sizeof(char)*SEQNUM*SEQLENGTH);
	char *filename_sequences = "../../dataset/DATASET_DNA_3.fasta";
	FILE* f_data_sequences = fopen(filename_sequences, "r");
	if ( f_data_sequences == NULL ) {
		printf( "File not found: %s\n", filename_sequences );
		exit(1);
	}
	while ( fgets(seqLine, sizeof(seqLine), f_data_sequences) != NULL ) {
		if ( seqLine[0] != '>' ) {
			for ( int i = 0; i < SEQLENGTH; i ++ ) {
				sequences[seqIdx++] = seqLine[i];
			}
		}
	}
	fclose(f_data_sequences);
	printf( "Reading the sequence file finished!\n" );
	fflush( stdout );
	//--------------------------------------------------------------------------------------------
	// Run GibbsSampler
	//--------------------------------------------------------------------------------------------
	char *bestMotifs = (char*)malloc(sizeof(char)*SEQNUM*MOTIFLENGTH);
	printf( "Motif Finder Started!\n" );
	fflush( stdout );
	double processStart = timeCheckerCPU();
	gibbsSamplerWrapper(bestMotifs, sequences);
	double processFinish = timeCheckerCPU();
	printf( "Motif Finder Finished!\n" );
	fflush( stdout );
	double processTime = processFinish - processStart;
	//--------------------------------------------------------------------------------------------
	// Print elapsed time
	//--------------------------------------------------------------------------------------------
	printf( "Elapsed Time: %.8f\n", processTime );
	fflush( stdout );
	//int s = 0;
	//for ( int i = 0; i < posIdx; i ++ ) {
	//	s += positionQ[i];
	//}
	//printf("%f\n", (float)s / (float)posIdx );
	//--------------------------------------------------------------------------------------------
	// Print the results
	//--------------------------------------------------------------------------------------------
	char *filename_result = "DATASET_DNA_3_result.fasta";
	FILE *f_data_result = fopen(filename_result, "w");
	printf( "Storing Motifs Started!\n" );
	fflush( stdout );
	for ( int i = 0; i < SEQNUM; i ++ ) {
		fprintf(f_data_result, ">%d\n", i);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			fputc(bestMotifs[MOTIFLENGTH*i + j], f_data_result);
		}
		fputc('\n', f_data_result);
	}
	fclose(f_data_result);
	printf( "Storing Motifs Finished!\n" );
	fflush( stdout );

	// Free
	free(sequences);
	free(bestMotifs);

	return 0;
}
