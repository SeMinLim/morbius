#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


#define POLY_MASK_32 0xB4BCD35C
#define POLY_MASK_31 0x7A5BC2E3


// Parameter setting
#define SEQSNUM 26
#define SEQLENGTH 111
#define NUMSEEDS 10
#define NUMITER 20
#define MOTIFLENGTH 11


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Random number generator based on LFSR
uint32_t lfsr32 = 0xACE8F;
uint32_t lfsr31 = 0x23456789;
uint32_t shift_lfsr( uint32_t *lfsr, uint32_t polynomial_mask ) {
	uint32_t feedback = *lfsr & 1;
	*lfsr >>= 1;
	if ( feedback == 1 ) *lfsr ^= polynomial_mask;
	return *lfsr;
}
uint32_t rand_generator() {
	shift_lfsr(&lfsr32, POLY_MASK_32);
	uint32_t tmp_1 = (shift_lfsr(&lfsr32, POLY_MASK_32) ^ shift_lfsr(&lfsr31, POLY_MASK_31)) & 0xffff;
	uint32_t tmp_2 = tmp_1 << 31;
	uint32_t value = tmp_2 >> 31;
	return value;
}

// Get a score of motifs
int score( char *motifs ) {
	// Phase 1
	// Pick the most frequent letter at each motif position
	char pattern[MOTIFLENGTH];
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		int a = 0, c = 0, g = 0, t = 0;
		for ( int j = 0; j < SEQSNUM; j ++ ) {
			if ( motifs[MOTIFLENGTH*j+i] == 'A' ) a += 1;
			else if ( motifs[MOTIFLENGTH*j+i] == 'C' ) c += 1;
			else if ( motifs[MOTIFLENGTH*j+i] == 'G' ) g += 1;
			else if ( motifs[MOTIFLENGTH*j+i] == 'T' ) t += 1;
		}

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
	for ( int i = 0; i < SEQSNUM; i ++ ) {
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
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		int a = 1, c = 1, g = 1, t = 1;
		for ( int j = 0; j < SEQSNUM; j ++ ) {
			// Build PSSM based on the motifs except for a picked sequence's motif
			if ( j != seqIdx ) {
				if ( motifs[MOTIFLENGTH*j+i] == 'A' ) a += 1;
				else if ( motifs[MOTIFLENGTH*j+i] == 'C' ) c += 1;
				else if ( motifs[MOTIFLENGTH*j+i] == 'G' ) g += 1;
				else if ( motifs[MOTIFLENGTH*j+i] == 'T' ) t += 1;
			}
		}

		float total = (float)(a + c + g + t);
		pssm[i] = (float)a / total;
		pssm[MOTIFLENGTH*1+i] = (float)c / total;
		pssm[MOTIFLENGTH*2+i] = (float)g / total;
		pssm[MOTIFLENGTH*3+i] = (float)t / total;
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
	char motifs[SEQSNUM*MOTIFLENGTH];
	for ( int i = 0; i < SEQSNUM; i ++ ) {
		int sp = rand_generator() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			char c = sequences[SEQLENGTH*i + sp + j];
			motifs[MOTIFLENGTH*i + j] = c;
			results[MOTIFLENGTH*i + j] = c;
		}
	}
	int resultsScore = score(motifs);

	// Phase 2
	// Update the motifs over 20 x #sequence
	char sequence[SEQLENGTH];
	char updatedMotif[MOTIFLENGTH];
	float pssm[4*MOTIFLENGTH];
	for ( int i = 0; i < NUMITER; i ++ ) {
		for ( int j = 0; j < SEQSNUM; j ++ ) {
			// Pick a sequence sequencially
			for ( int k = 0; k < SEQLENGTH; k ++ ) {
				sequence[k] = sequences[SEQLENGTH*j + k];
			}
			
			// Make PSSM first
			makePSSM(pssm, motifs, j);
			
			// Go to profiling step then
			profile(updatedMotif, sequence, pssm);
			
			// Get an updated motif of the picked sequence
			for ( int k = 0; k < MOTIFLENGTH; k ++ ) {
				motifs[MOTIFLENGTH*j + k] = updatedMotif[k];
			}

			// Compare the scores & Update motifs or not based on the score
			int currentScore = score(motifs);
			if ( currentScore < resultsScore ) {
				for ( int k = 0; k < SEQSNUM*MOTIFLENGTH; k ++ ) {
					results[k] = motifs[k];
				}
				resultsScore = currentScore;
			}
		}
	}
}

// Gibbs Sampler Wrapper
void gibbsSamplerWrapper( char *bestMotifs, char *sequences ) {
	int bestScore = 0;
	
	char results[SEQSNUM*MOTIFLENGTH + 1];
	for ( int i = 0; i < NUMSEEDS; i ++ ) {
		gibbsSampler(results, sequences);
		if ( i == 0 ) {
			int currentScore = score(results);
			bestScore = currentScore;
			for ( int j = 0; j < SEQSNUM*MOTIFLENGTH; j ++ ) {
				bestMotifs[j] = results[j];
			}
		} else {
			int currentScore = score(results);
			if ( currentScore < bestScore ) {
				bestScore = currentScore;
				for ( int j = 0; j < SEQSNUM*MOTIFLENGTH; j ++ ) {
					bestMotifs[j] = results[j];
				}
			}
		}
	}
}


// Main
int main() {
	srand(time(NULL));

	// Read the sequence text file first
	int seqsIdx = 0;
	char line[1024];
	char *sequences = (char*)malloc(sizeof(char)*SEQSNUM*SEQLENGTH);
	char *filename_query = "../dataset/MM9GATA4.fasta";
	FILE* f_data_query = fopen(filename_query, "r");
	if ( f_data_query == NULL ) {
		printf( "File not found: %s\n", filename_query );
		exit(1);
	}
	while ( fgets(line, sizeof(line), f_data_query) != NULL ) {
		if ( line[0] != '>' ) {
			for ( int i = 0; i < 111; i ++ ) {
				sequences[seqsIdx++] = line[i];
			}
		}
	}
	fclose(f_data_query);
	printf( "Read the query sequence file done!\n" );
	fflush( stdout );

	// Read the answer text file then
	int ansIdx = 0;
	char *answers = (char*)malloc(sizeof(char)*SEQSNUM*MOTIFLENGTH);
	char *filename_answers = "../dataset/MM9GATA4SOLUTION.fasta";
	FILE* f_data_answers = fopen(filename_answers, "r");
	if ( f_data_answers == NULL ) {
		printf( "File not found: %s\n", filename_answers );
		exit(1);
	}
	while ( !feof(f_data_answers) ) {
		char c;
		fread(&c, sizeof(char), 1, f_data_answers);
		if ( c != '\n' ) {
			if ( isupper(c) ) answers[ansIdx++] = c;
		}
	}
	fclose(f_data_answers);
	printf( "Read the answer motif file done!\n" );
	fflush( stdout );

	// Run GibbsSampler
	char *bestMotifs = (char*)malloc(sizeof(char)*SEQSNUM*MOTIFLENGTH);
	printf( "Motif Finder Started!\n" );
	fflush( stdout );
	double processStart = timeCheckerCPU();
	gibbsSamplerWrapper(bestMotifs, sequences);
	double processFinish = timeCheckerCPU();
	printf( "Motif Finder Finished!\n" );
	fflush( stdout );
	double processTime = processFinish - processStart;

	// Get accuracy & Print elapsed time
	int match = 0;
	char bestMotif[MOTIFLENGTH+1];
	char answer[MOTIFLENGTH+1];
	printf( "Motif Finding via Gibbs Sampling\n" );
	fflush( stdout );
	for ( int i = 0; i < SEQSNUM; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			bestMotif[j] = bestMotifs[MOTIFLENGTH*i + j];
			answer[j] = answers[MOTIFLENGTH*i + j];
		}
		answer[MOTIFLENGTH] = '\0';
		bestMotif[MOTIFLENGTH] = '\0';
		
		if ( strncmp(bestMotif, answer, 11) == 0 ) match++;
	}
	float accuracy = ((float)match / (float)SEQSNUM) * 100;
	printf( "Accuracy: %f\n", accuracy );
	printf( "Elapsed Time: %.8f\n", processTime );
	fflush( stdout );

	// Print the results
	char *filename_result = "MM9GATA4_result.fasta";
	FILE *f_data_result = fopen(filename_result, "w");
	printf( "Storing Motifs Started!\n" );
	fflush( stdout );
	for ( int i = 0; i < SEQSNUM; i ++ ) {
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
	free(answers);
	free(bestMotifs);

	return 0;
}