#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Hamming distance
int hammingDistance( char *motif, char *pattern, int motifLength ) {
	int diffs = 0;
	for ( int i = 0; i < motifLength; i ++ ) {
		if ( motif[i] != pattern[i] ) diffs += 1;
	}

	return diffs;
}

// Get a score of motifs
int score( char *motifs, int motifLength, int numSeqs ) {
	char *pattern = (char*)malloc(sizeof(char)*motifLength);
	for ( int i = 0; i < motifLength; i ++ ) {
		int a = 0, c = 0, g = 0, t = 0;
		for ( int j = 0; j < numSeqs; j ++ ) {
			if ( motifs[motifLength*j+i] == 'A' ) a += 1;
			else if ( motifs[motifLength*j+i] == 'C' ) c += 1;
			else if ( motifs[motifLength*j+i] == 'G' ) g += 1;
			else if ( motifs[motifLength*j+i] == 'T' ) t += 1;
		}

		if ( a >= c and a >= g and a >= t ) pattern[i] = 'A';
		else if ( c >= g and c >= t ) pattern[i] = 'C';
		else if ( g >= t ) pattern[i] = 'G';
		else pattern[i] = 'T';
	}

	int score = 0;
	char *motif = (char*)malloc(sizeof(char)*motifLength);
	for ( int i = 0; i < numSeqs; i ++ ) {
		for ( int j = 0; j < motifLength; j ++ ) {
			motif[j] = motifs[i*motifLength+j];
		}
		score += hammingDistance(motif, pattern, motifLength);
	}

	free(pattern);
	free(motif);

	return score;
}

// Make position specific score matrix
void makePSSM( float *pssm, char *subsetMotifs, int motifLength, int numSubsetMotifs ) {
	for ( int i = 0; i < motifLength; i ++ ) {
		int a = 1, c = 1, g = 1, t = 1;
		for ( int j = 0; j < numSubsetMotifs; j ++ ) {
			if ( subsetMotifs[motifLength*j+i] == 'A' ) a += 1;
			else if ( subsetMotifs[motifLength*j+i] == 'C' ) c += 1;
			else if ( subsetMotifs[motifLength*j+i] == 'G' ) g += 1;
			else if ( subsetMotifs[motifLength*j+i] == 'T' ) t += 1;
		}

		float total = (float)(a + c + g + t);
		pssm[i] = (float)a / total;
		pssm[motifLength*1+i] = (float)c / total;
		pssm[motifLength*2+i] = (float)g / total;
		pssm[motifLength*3+i] = (float)t / total;
	}
}

// Profiling
void profile( char *updatedMotif, char *subsetMotifs, char *sequence, int motifLength, int numSubsetMotifs ) {
	float *pssm = (float*)malloc(sizeof(float)*4*motifLength);
	makePSSM(pssm, subsetMotifs, motifLength, numSubsetMotifs);

	int sizeProbs = strlen(sequence) - motifLength + 1;
	float sumProbs = 0.0;
	float *probs = (float*)malloc(sizeof(float)*sizeProbs);
	for ( int i = 0; i < sizeProbs; i ++ ) {
		float prob = 1.0;
		for ( int j = 0; j < motifLength; j ++ ) {
			if ( sequence[i+j] == 'A' ) prob *= pssm[j];
			else if ( sequence[i+j] == 'C' ) prob *= pssm[motifLength*1+j];
			else if ( sequence[i+j] == 'G' ) prob *= pssm[motifLength*2+j];
			else if ( sequence[i+j] == 'T' ) prob *= pssm[motifLength*3+j];
		}

		probs[i] = prob;
		sumProbs += prob;
	}

	float partialSum = 0.0;
	float randVal = (float)((rand()%10) / 10);
	int position = 0;
	for ( int i = 0; i < sizeProbs; i ++ ) {
		partialSum += probs[i];
		if ( (partialSum/sumProbs) >= randVal ) {
			position = i;
			break;
		}
	}

	for ( int i = 0; i < motifLength; i ++ ) {
		updatedMotif[i] = sequence[position+i];
	}

	free(pssm);
	free(probs);
}

// Gibbs Sampler
void gibbsSampler( char *results, char *sequences, int numIter, int motifLength, int numSeqs, int lengthSeq ) {
	char *motifs = (char*)malloc(sizeof(char)*numSeqs*motifLength);
	for ( int i = 0; i < numSeqs; i ++ ) {
		int sp = rand() % (lengthSeq - motifLength + 1);
		for ( int j = 0; j < motifLength; j ++ ) {
			motifs[motifLength*i+j] = sequences[motifLength*i+sp+j];
		}
	}

	strcpy(results, motifs);
	int resultsScore = score(motifs, motifLength, numSeqs);

	char *subsetMotifs = (char*)malloc(sizeof(char)*(numSeqs-1)*motifLength);
	char *updatedMotif = (char*)malloc(sizeof(char)*motifLength);
	for ( int i = 0; i < numIter; i ++ ) {
		int j = rand() % numSeqs;
		int idx = 0;
		if ( j == 0 ) {
			for ( int k = 1; k < numSeqs; k ++ ) {
				for ( int l = 0; l < motifLength; l ++ ) {
					subsetMotifs[motifLength*idx+l] = motifs[motifLength*k+l];
				}
				idx++;
			}
		} else {
			for ( int k = 0; k < j; k ++ ) {
				for ( int l = 0; i < motifLength; l ++ ) {
					subsetMotifs[motifLength*idx+l] = motifs[motifLength*k+l];
				}
				idx++;
			}
			for ( int k = j + 1; k < numSeqs; k ++ ) {
				for ( int l = 0; l < motifLength; l ++ ) {
					subsetMotifs[motifLength*idx+l] = motifs[motifLength*k+l];
				}
				idx++;
			}
		}
		
		char *sequence = (char*)malloc(sizeof(char)*lengthSeq);
		for ( int i = 0; i < lengthSeq; i ++ ) {
			sequence[i] = sequences[lengthSeq*j + i];
		}
		profile(updatedMotif, subsetMotifs, sequence, motifLength, (numSeqs-1));
		for ( int i = 0; i < motifLength; i ++ ) {
			motifs[motifLength*j+i] = updatedMotif[i];
		}

		int currentScore = score(motifs, motifLength, numSeqs);
		if ( currentScore < resultsScore ) {
			strcpy(results, motifs);
			resultsScore = currentScore;
		}
	}

	free(motifs);
	free(subsetMotifs);
	free(updatedMotif);
}

// Gibbs Sampler Wrapper
void multipleSeedsGibbsSampling( char *bestMotifs, char *sequences, 
				 int numSeeds, int numIter, 
				 int motifLength, int numSeqs, int lengthSeq ) {
	int bestScore = 0;
	
	char *results = (char*)malloc(sizeof(char)*numSeqs*motifLength);
	for ( int i = 0; i < numSeeds; i ++ ) {
		gibbsSampler(results, sequences, numIter, motifLength, numSeqs, lengthSeq);
		if ( i == 0 ) {
			int currentScore = score(results, motifLength, numSeqs);
			bestScore = currentScore;
			strcpy(bestMotifs, results);
		} else {
			int currentScore = score(results, motifLength, numSeqs);
			if ( currentScore < bestScore ) {
				bestScore = currentScore;
				strcpy(bestMotifs, results);
			}
		}
	}

	free(results);
}

// Main
int main() {
	srand(time(NULL));

	// Initial parameter setting
	int numSeq = 26;
	int lengthSeq = 111;
	int numSeeds = 20;
	int numIter = 1000;
	int motifLength = 11;

	// Read the sequence text file first
	char *line = (char*)malloc(sizeof(char)*lengthSeq);
	char *sequences = (char*)malloc(sizeof(char)*numSeq*lengthSeq);
	char *filename_query = "../dataset/mm9Gata4MotifCollection.txt";
	FILE* f_data_query = fopen(filename_query, "r");
	if ( f_data_query == NULL ) {
		printf( "File not found: %s\n", filename_query );
		exit(1);
	}
	for ( int i = 0; i < numSeq; i ++ ) {
		fgets(line, lengthSeq, f_data_query);
		for ( int j = 0; j < lengthSeq; j ++ ) {
			sequences[lengthSeq*i+j] = line[j];
		}
	}
	fclose(f_data_query);
	printf( "Read the query sequence file done!\n" );

	// Read the answer text file then
	int cntAns = 0;
	char *answers = (char*)malloc(sizeof(char)*numSeq*motifLength);
	char *filename_answers = "../dataset/mm9Gata4Solutions.txt";
	FILE* f_data_answers = fopen(filename_answers, "r");
	if ( f_data_answers == NULL ) {
		printf( "File not found: %s\n", filename_answers );
		exit(1);
	}
	for ( int i = 0; i < numSeq; i ++ ) {
		fgets(line, lengthSeq, f_data_answers);
		for ( int j = 0; j < lengthSeq; j ++ ) {
			if ( isupper(line[j]) ) {
				answers[cntAns] = line[j];
				cntAns++;
			}
		}
	}
	fclose(f_data_answers);
	printf( "Read the answer motif file done!\n" );

	// Run GibbsSampler
	char *bestMotifs = (char*)malloc(sizeof(char)*numSeq*motifLength);
	double processStart = timeCheckerCPU();
	multipleSeedsGibbsSampling(bestMotifs, sequences, numSeeds, numIter, motifLength, numSeq, lengthSeq);
	double processFinish = timeCheckerCPU();
	double processTime = processFinish - processStart;

	// Print the results
	printf( "Motif Finding via Gibbs Sampling\n" );
	printf( "System Best Pick        Match/Unmatch       Answer\n" );
	char *bestMotif = (char*)malloc(sizeof(char)*motifLength);
	char *answer = (char*)malloc(sizeof(char)*motifLength);
	for ( int i = 0; i < numSeq; i ++ ) {
		for ( int j = 0 ; j < motifLength; j ++ ) {
			bestMotif[j] = bestMotifs[motifLength*i+j];
			answer[j] = answers[motifLength*i+j];
		}
		if ( strcmp(bestMotif, answer) == 0 ) {
			printf( "%s   Match   %s\n", bestMotif, answer );
		} else {
			printf( "%s   Match   %s\n", bestMotif, answer );
		}
	}
	printf( "Elapsed Time: %.8f\n", processTime );

	// Free
	free(line);
	free(sequences);
	free(answers);
	free(bestMotifs);
	free(answer);
	free(bestMotif);

	return 0;
}
