#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


// Parameter setting
#define SEQNUM 32768
#define SEQLENGTH 1000
#define MOTIFLENGTH 16
#define NUMSEEDS 3
#define NUMITER 1


// Elapsed time checker
static inline double timeCheckerCPU(void) {
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
        return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

// Get a score of motifs
int score( char *motifs ) {
	// Phase 1
	// Pick the most frequent letter at each motif position
	char pattern[MOTIFLENGTH];
	for ( int cnt_1 = 0; cnt_1 < MOTIFLENGTH; cnt_1 ++ ) {
		int m = 0, t = 0, n = 0, k = 0, s = 0;
		int r = 0, v = 0, a = 0, d = 0, e = 0;
		int g = 0, f = 0, l = 0, y = 0, c = 0;
		int w = 0, p = 0, h = 0, q = 0, i = 0;	
		for ( int cnt_2 = 0; cnt_2 < SEQNUM; cnt_2 ++ ) {
			if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'M' ) m += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'T' ) t += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'N' ) n += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'K' ) k += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'S' ) s += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'R' ) r += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'V' ) v += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'A' ) a += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'D' ) d += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'E' ) e += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'G' ) g += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'F' ) f += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'L' ) l += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'Y' ) y += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'C' ) c += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'W' ) w += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'P' ) p += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'H' ) h += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'Q' ) q += 1;
			else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'I' ) i += 1;
		}

		if ( m >= t && m >= n && m >= k && m >= s && m >= r &&
		     m >= v && m >= a && m >= d && m >= e && m >= g &&
		     m >= f && m >= l && m >= y && m >= c && m >= w &&
		     m >= p && m >= h && m >= q && m >= i ) {
			pattern[cnt_1] = 'M';
		} else if ( t >= n && t >= k && t >= s && t >= r && t >= v &&
			    t >= a && t >= d && t >= e && t >= g && t >= f &&
			    t >= l && t >= y && t >= c && t >= w && t >= p &&
			    t >= h && t >= q && t >= i ) {
			pattern[cnt_1] = 'T';
		} else if ( n >= k && n >= s && n >= r && n >= v && n >= a &&
			    n >= d && n >= e && n >= g && n >= f && n >= l &&
			    n >= y && n >= c && n >= w && n >= p && n >= h &&
			    n >= q && n >= i ) {
			pattern[cnt_1] = 'N';
		} else if ( k >= s && k >= r && k >= v && k >= a && k >= d &&
			    k >= e && k >= g && k >= f && k >= l && k >= y &&
			    k >= c && k >= w && k >= p && k >= h && k >= q &&
			    k >= i ) {
			pattern[cnt_1] = 'K';
		} else if ( s >= r && s >= v && s >= a && s >= d && s >= e &&
			    s >= g && s >= f && s >= l && s >= y && s >= c &&
			    s >= w && s >= p && s >= h && s >= q && s >= i ) {
			pattern[cnt_1] = 'S';
		} else if ( r >= v && r >= a && r >= d && r >= e && r >= g &&
			    r >= f && r >= l && r >= y && r >= c && r >= w &&
			    r >= p && r >= h && r >= q && r >= i ) {
			pattern[cnt_1] = 'R';
		} else if ( v >= a && v >= d && v >= e && v >= g && v >= f &&
			    v >= l && v >= y && v >= c && v >= w && v >= p &&
			    v >= h && v >= q && v >= i ) {
			pattern[cnt_1] = 'V';
		} else if ( a >= d && a >= e && a >= g && a >= f && a >= l &&
			    a >= y && a >= c && a >= w && a >= p && a >= h &&
			    a >= q && a >= i ) {
			pattern[cnt_1] = 'A';
		} else if ( d >= e && d >= g && d >= f && d >= l && d >= y &&
			    d >= c && d >= w && d >= p && d >= h && d >= q &&
			    d >= i ) {
			pattern[cnt_1] = 'D';
		} else if ( e >= g && e >= f && e >= l && e >= y && e >= c &&
			    e >= w && e >= p && e >= h && e >= q && e >= i ) {
			pattern[cnt_1] = 'E';
		} else if ( g >= f && g >= l && g >= y && g >= c && g >= w &&
			    g >= p && g >= h && g >= q && g >= i ) {
			pattern[cnt_1] = 'G';
		} else if ( f >= l && f >= y && f >= c && f >= w && f >= p &&
			    f >= h && f >= q && f >= i ) {
			pattern[cnt_1] = 'F';
		} else if ( l >= y && l >= c && l >= w && l >= p && l >= h &&
			    l >= q && l >= i ) {
			pattern[cnt_1] = 'L';
		} else if ( y >= c && y >= w && y >= p && y >= h && y >= q &&
			    y >= i ) {
			pattern[cnt_1] = 'Y';
		} else if ( c >= w && c >= p && c >= h && c >= q && c >= i ) {
			pattern[cnt_1] = 'C';
		} else if ( w >= p && w >= h && w >= q && w >= i ) {
			pattern[cnt_1] = 'W';
		} else if ( p >= h && p >= q && p >= i ) {
			pattern[cnt_1] = 'P';
		} else if ( h >= q && h >= i ) {
			pattern[cnt_1] = 'H';
		} else if ( g >= i ) {
			pattern[cnt_1] = 'G';
		} else {
			pattern[cnt_1] = 'I';
		}
	}

	// Phase 2
	// Compare between each motif and the picked string 
	// Get the score via Hamming Distance 
	int score = 0;
	char motif[MOTIFLENGTH];
	for ( int cnt_1 = 0; cnt_1 < SEQNUM; cnt_1 ++ ) {
		for ( int cnt_2 = 0; cnt_2 < MOTIFLENGTH; cnt_2 ++ ) {
			motif[cnt_2] = motifs[MOTIFLENGTH*cnt_1 + cnt_2];
		}
	
		for ( int cnt_2 = 0; cnt_2 < MOTIFLENGTH; cnt_2 ++ ) {
			if ( motif[cnt_2] != pattern[cnt_2] ) score += 1;
		}
	}

	return score;
}

// Make position specific score matrix
void makePSSM( float *pssm, char *motifs, int seqIdx ) {
	for ( int cnt_1 = 0; cnt_1 < MOTIFLENGTH; cnt_1 ++ ) {
		int m = 1, t = 1, n = 1, k = 1, s = 1;
		int r = 1, v = 1, a = 1, d = 1, e = 1;
		int g = 1, f = 1, l = 1, y = 1, c = 1;
		int w = 1, p = 1, h = 1, q = 1, i = 1;	
		for ( int cnt_2 = 0; cnt_2 < SEQNUM; cnt_2 ++ ) {
			// Build PSSM based on the motifs except for a picked sequence's motif
			if ( cnt_2 != seqIdx ) {
				if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'M' ) m += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'T' ) t += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'N' ) n += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'K' ) k += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'S' ) s += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'R' ) r += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'V' ) v += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'A' ) a += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'D' ) d += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'E' ) e += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'G' ) g += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'F' ) f += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'L' ) l += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'Y' ) y += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'C' ) c += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'W' ) w += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'P' ) p += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'H' ) h += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'Q' ) q += 1;
				else if ( motifs[MOTIFLENGTH*cnt_2 + cnt_1] == 'I' ) i += 1;
			}
		}

		float total = (float)(m + t + n + k + s + r + v + a + d + e +
				      g + f + l + y + c + w + p + h + q + i);
		pssm[cnt_1] = (float)m / total;
		pssm[MOTIFLENGTH*1 + cnt_1] = (float)t / total;
		pssm[MOTIFLENGTH*2 + cnt_1] = (float)n / total;
		pssm[MOTIFLENGTH*3 + cnt_1] = (float)k / total;
		pssm[MOTIFLENGTH*4 + cnt_1] = (float)s / total;
		pssm[MOTIFLENGTH*5 + cnt_1] = (float)r / total;
		pssm[MOTIFLENGTH*6 + cnt_1] = (float)v / total;
		pssm[MOTIFLENGTH*7 + cnt_1] = (float)a / total;
		pssm[MOTIFLENGTH*8 + cnt_1] = (float)d / total;
		pssm[MOTIFLENGTH*9 + cnt_1] = (float)e / total;
		pssm[MOTIFLENGTH*10 + cnt_1] = (float)g / total;
		pssm[MOTIFLENGTH*11 + cnt_1] = (float)f / total;
		pssm[MOTIFLENGTH*12 + cnt_1] = (float)l / total;
		pssm[MOTIFLENGTH*13 + cnt_1] = (float)y / total;
		pssm[MOTIFLENGTH*14 + cnt_1] = (float)c / total;
		pssm[MOTIFLENGTH*15 + cnt_1] = (float)w / total;
		pssm[MOTIFLENGTH*16 + cnt_1] = (float)p / total;
		pssm[MOTIFLENGTH*17 + cnt_1] = (float)h / total;
		pssm[MOTIFLENGTH*18 + cnt_1] = (float)q / total;
		pssm[MOTIFLENGTH*19 + cnt_1] = (float)i / total;
	}
}

// Profiling
void profile( char *updatedMotif, char *sequence, float *pssm ) {
	int sizeProbs = SEQLENGTH - MOTIFLENGTH + 1;

	// Phase 1
	// Get probabilities
	float sumProbs = 0.0;
	float probs[sizeProbs];
	for ( int cnt_1 = 0; cnt_1 < sizeProbs; cnt_1 ++ ) {
		float prob = 1.0;
		for ( int cnt_2 = 0; cnt_2 < MOTIFLENGTH; cnt_2 ++ ) {
			if ( sequence[cnt_1 + cnt_2] == 'M' ) prob *= pssm[cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'T' ) prob *= pssm[MOTIFLENGTH*1 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'N' ) prob *= pssm[MOTIFLENGTH*2 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'K' ) prob *= pssm[MOTIFLENGTH*3 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'S' ) prob *= pssm[MOTIFLENGTH*4 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'R' ) prob *= pssm[MOTIFLENGTH*5 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'V' ) prob *= pssm[MOTIFLENGTH*6 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'A' ) prob *= pssm[MOTIFLENGTH*7 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'D' ) prob *= pssm[MOTIFLENGTH*8 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'E' ) prob *= pssm[MOTIFLENGTH*9 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'G' ) prob *= pssm[MOTIFLENGTH*10 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'F' ) prob *= pssm[MOTIFLENGTH*11 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'L' ) prob *= pssm[MOTIFLENGTH*12 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'Y' ) prob *= pssm[MOTIFLENGTH*13 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'C' ) prob *= pssm[MOTIFLENGTH*14 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'W' ) prob *= pssm[MOTIFLENGTH*15 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'P' ) prob *= pssm[MOTIFLENGTH*16 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'H' ) prob *= pssm[MOTIFLENGTH*17 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'Q' ) prob *= pssm[MOTIFLENGTH*18 + cnt_2];
			else if ( sequence[cnt_1 + cnt_2] == 'I' ) prob *= pssm[MOTIFLENGTH*19 + cnt_2];
		}

		probs[cnt_1] = prob;
		sumProbs += prob;
	}

	// Phase 2
	// Randomly designate the updated motif
	float partialSum = 0.0;
	float randVal = (float)rand()/RAND_MAX;
	int position = 0;
	for ( int cnt = 0; cnt < sizeProbs; cnt ++ ) {
		partialSum += probs[cnt];
		if ( (partialSum/sumProbs) >= randVal ) {
			position = cnt;
			break;
		}
	}

	for ( int cnt = 0; cnt < MOTIFLENGTH; cnt ++ ) {
		updatedMotif[cnt] = sequence[position + cnt];
	}
}

// Gibbs Sampler
void gibbsSampler( char *results, char *sequences ) {
	// Phase 1
	// Designate the motifs randomly first
	char motifs[SEQNUM*MOTIFLENGTH];
	for ( int cnt_1 = 0; cnt_1 < SEQNUM; cnt_1 ++ ) {
		int sp = rand() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int cnt_2 = 0; cnt_2 < MOTIFLENGTH; cnt_2 ++ ) {
			char c = sequences[SEQLENGTH*cnt_1 + sp + cnt_2];
			motifs[MOTIFLENGTH*cnt_1 + cnt_2] = c;
			results[MOTIFLENGTH*cnt_1 + cnt_2] = c;
		}
	}
	int resultsScore = score(motifs);

	// Phase 2
	// Update the motifs over #sequence
	char sequence[SEQLENGTH];
	char updatedMotif[MOTIFLENGTH];
	float pssm[20*MOTIFLENGTH];
	for ( int cnt_1 = 0; cnt_1 < NUMITER; cnt_1 ++ ) {
		for ( int cnt_2 = 0; cnt_2 < SEQNUM; cnt_2 ++ ) {
			// Pick a sequence sequencially
			for ( int cnt_3 = 0; cnt_3 < SEQLENGTH; cnt_3 ++ ) {
				sequence[cnt_3] = sequences[SEQLENGTH*cnt_2 + cnt_3];
			}
			
			// Make PSSM first
			makePSSM(pssm, motifs, cnt_2);
			
			// Go to profiling step then
			profile(updatedMotif, sequence, pssm);
			
			// Get an updated motif of the picked sequence
			for ( int cnt_3 = 0; cnt_3 < MOTIFLENGTH; cnt_3 ++ ) {
				motifs[MOTIFLENGTH*cnt_2 + cnt_3] = updatedMotif[cnt_3];
			}

			// Compare the scores & Update motifs or not based on the score
			int currentScore = score(motifs);
			if ( currentScore < resultsScore ) {
				for ( int cnt_3 = 0; cnt_3 < SEQNUM*MOTIFLENGTH; cnt_3 ++ ) {
					results[cnt_3] = motifs[cnt_3];
				}
				resultsScore = currentScore;
			}
		}
	}
}

// Gibbs Sampler Wrapper
void gibbsSamplerWrapper( char *bestMotifs, char *sequences ) {
	int bestScore = 0;
	
	char results[SEQNUM*MOTIFLENGTH + 1];
	for ( int cnt_1 = 0; cnt_1 < NUMSEEDS; cnt_1 ++ ) {
		gibbsSampler(results, sequences);
		if ( cnt_1 == 0 ) {
			int currentScore = score(results);
			bestScore = currentScore;
			for ( int cnt_2 = 0; cnt_2 < SEQNUM*MOTIFLENGTH; cnt_2 ++ ) {
				bestMotifs[cnt_2] = results[cnt_2];
			}
		} else {
			int currentScore = score(results);
			if ( currentScore < bestScore ) {
				bestScore = currentScore;
				for ( int cnt_2 = 0; cnt_2 < SEQNUM*MOTIFLENGTH; cnt_2 ++ ) {
					bestMotifs[cnt_2] = results[cnt_2];
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
	char *filename_sequences = "../../dataset/DATASET_PROTEIN_1.fasta";
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
	//--------------------------------------------------------------------------------------------
	// Print the results
	//--------------------------------------------------------------------------------------------
	char *filename_result = "DATASET_PROTEIN_1_result.fasta";
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
