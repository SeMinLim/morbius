#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


// Parameter setting
#define SEQSNUM 26
#define SEQLENGTH 111
#define NUMSEEDS 20
#define NUMITER 1000
#define MOTIFLENGTH 11


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Hamming distance
int hammingDistance( char *motif, char *pattern ) {
	int diffs = 0;
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		if ( motif[i] != pattern[i] ) diffs += 1;
	}

	return diffs;
}

// Get a score of motifs
int score( char *motifs ) {
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

	int score = 0;
	char motif[MOTIFLENGTH];
	for ( int i = 0; i < SEQSNUM; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			motif[j] = motifs[i*MOTIFLENGTH+j];
		}
		score += hammingDistance(motif, pattern);
	}

	return score;
}

// Make position specific score matrix
void makePSSM( float *pssm, char *subsetMotifs ) {
	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		int a = 1, c = 1, g = 1, t = 1;
		for ( int j = 0; j < SEQSNUM - 1; j ++ ) {
			if ( subsetMotifs[MOTIFLENGTH*j+i] == 'A' ) a += 1;
			else if ( subsetMotifs[MOTIFLENGTH*j+i] == 'C' ) c += 1;
			else if ( subsetMotifs[MOTIFLENGTH*j+i] == 'G' ) g += 1;
			else if ( subsetMotifs[MOTIFLENGTH*j+i] == 'T' ) t += 1;
		}

		float total = (float)(a + c + g + t);
		pssm[i] = (float)a / total;
		pssm[MOTIFLENGTH*1+i] = (float)c / total;
		pssm[MOTIFLENGTH*2+i] = (float)g / total;
		pssm[MOTIFLENGTH*3+i] = (float)t / total;
	}
}

// Profiling
void profile( char *updatedMotif, char *subsetMotifs, char *sequences, 
	      int seqIdx ) {
	float pssm[4*MOTIFLENGTH];
	makePSSM(pssm, subsetMotifs);

	int sizeProbs = SEQLENGTH - MOTIFLENGTH + 1;
	float sumProbs = 0.0;
	float probs[sizeProbs];
	for ( int i = 0; i < sizeProbs; i ++ ) {
		float prob = 1.0;
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			if ( sequences[SEQLENGTH*seqIdx+i+j] == 'A' ) prob *= pssm[j];
			else if ( sequences[SEQLENGTH*seqIdx+i+j] == 'C' ) prob *= pssm[MOTIFLENGTH*1+j];
			else if ( sequences[SEQLENGTH*seqIdx+i+j] == 'G' ) prob *= pssm[MOTIFLENGTH*2+j];
			else if ( sequences[SEQLENGTH*seqIdx+i+j] == 'T' ) prob *= pssm[MOTIFLENGTH*3+j];
		}

		probs[i] = prob;
		sumProbs += prob;
	}

	float partialSum = 0.0;
	float randVal = (float) rand()/RAND_MAX;
	int position = 0;
	for ( int i = 0; i < sizeProbs; i ++ ) {
		partialSum += probs[i];
		if ( (partialSum/sumProbs) >= randVal ) {
			position = i;
			break;
		}
	}

	for ( int i = 0; i < MOTIFLENGTH; i ++ ) {
		updatedMotif[i] = sequences[SEQLENGTH*seqIdx+position+i];
	}
}

// Gibbs Sampler
void gibbsSampler( char *results, char *sequences ) {
	char motifs[SEQSNUM*MOTIFLENGTH];
	for ( int i = 0; i < SEQSNUM; i ++ ) {
		int sp = rand() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			motifs[MOTIFLENGTH*i+j] = sequences[SEQLENGTH*i+sp+j];
		}
	}
	
	strcpy(results, motifs);
	int resultsScore = score(motifs);
	
	char subsetMotifs[(SEQSNUM-1)*MOTIFLENGTH];
	char updatedMotif[MOTIFLENGTH];
	for ( int i = 0; i < NUMITER; i ++ ) {
		int j = rand() % SEQSNUM;
		int idx = 0;
		for ( int k = 0; k < j; k ++ ) {
			for ( int l = 0; l < MOTIFLENGTH; l ++ ) {
				subsetMotifs[MOTIFLENGTH*idx+l] = motifs[MOTIFLENGTH*k+l];
			}
			idx++;
		}
		for ( int k = j + 1; k < SEQSNUM; k ++ ) {
			for ( int l = 0; l < MOTIFLENGTH; l ++ ) {
				subsetMotifs[MOTIFLENGTH*idx+l] = motifs[MOTIFLENGTH*k+l];
			}
			idx++;
		}
		
		profile(updatedMotif, subsetMotifs, sequences, j);
		for ( int m = 0; m < MOTIFLENGTH; m ++ ) {
			motifs[MOTIFLENGTH*j+m] = updatedMotif[m];
		}

		int currentScore = score(motifs);
		if ( currentScore < resultsScore ) {
			strcpy(results, motifs);
			resultsScore = currentScore;
		}
	}
}

// Gibbs Sampler Wrapper
void multipleSeedsGibbsSampling( char *bestMotifs, char *sequences ) {
	int bestScore = 0;
	
	char results[SEQSNUM*MOTIFLENGTH];
	for ( int i = 0; i < NUMSEEDS; i ++ ) {
		gibbsSampler(results, sequences);
		if ( i == 0 ) {
			int currentScore = score(results);
			bestScore = currentScore;
			strcpy(bestMotifs, results);
		} else {
			int currentScore = score(results);
			if ( currentScore < bestScore ) {
				bestScore = currentScore;
				strcpy(bestMotifs, results);
			}
		}
	}
}

// Main
int main() {
	srand(time(NULL));

	// Read the sequence text file first
	int seqsIdx = 0;
	char *sequences = (char*)malloc(sizeof(char)*SEQSNUM*SEQLENGTH);
	char *filename_query = "../dataset/mm9Gata4MotifCollection.txt";
	FILE* f_data_query = fopen(filename_query, "r");
	if ( f_data_query == NULL ) {
		printf( "File not found: %s\n", filename_query );
		exit(1);
	}
	while ( !feof(f_data_query) ) {
		char c;
		fread(&c, sizeof(char), 1, f_data_query);
		if ( c != '\n' ) sequences[seqsIdx++] = c;
	}
	fclose(f_data_query);
	printf( "Read the query sequence file done!\n" );

	// Read the answer text file then
	int ansIdx = 0;
	char *answers = (char*)malloc(sizeof(char)*SEQSNUM*MOTIFLENGTH);
	char *filename_answers = "../dataset/mm9Gata4Solutions.txt";
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

	// Run GibbsSampler
	char *bestMotifs = (char*)malloc(sizeof(char)*SEQSNUM*MOTIFLENGTH);
	double processStart = timeCheckerCPU();
	multipleSeedsGibbsSampling(bestMotifs, sequences);
	double processFinish = timeCheckerCPU();
	double processTime = processFinish - processStart;

	// Print the results
	char bestMotif[MOTIFLENGTH+1];
	char answer[MOTIFLENGTH+1];
	printf( "Motif Finding via Gibbs Sampling\n" );
	printf( "System Best Pick       Match/Unmatch       Answer\n" );
	for ( int i = 0; i < SEQSNUM; i ++ ) {
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			bestMotif[j] = bestMotifs[MOTIFLENGTH*i + j];
			answer[j] = answers[MOTIFLENGTH*i + j];
		}
		answer[MOTIFLENGTH] = '\0';
		bestMotif[MOTIFLENGTH] = '\0';
		
		if ( strncmp(bestMotif, answer, 11) == 0 ) {
			printf( "%s            Match               %s\n", bestMotif, answer );
		} else {
			printf( "%s            Unmatch             %s\n", bestMotif, answer );
		}
		
	}
	printf( "Elapsed Time: %.8f\n", processTime );

	// Free
	free(sequences);
	free(answers);
	free(bestMotifs);

	return 0;
}
