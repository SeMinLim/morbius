#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


// Parameter setting
#define SEQNUM 56054
#define SEQLENGTH 1000
#define MOTIFNUM 5098
#define MOTIFLENGTH 115


// Main
int main() {
	// Read the sequence text file first
	int seqIdx = 0;
	char seqLine[1024];
	char *sequences = (char*)malloc(sizeof(char)*SEQNUM*SEQLENGTH);
	char *filename_sqnc = "../UPSTREAM1000.fasta";
	FILE* f_data_sqnc = fopen(filename_sqnc, "r");
	if ( f_data_sqnc == NULL ) {
		printf( "File not found: %s\n", filename_sqnc );
		exit(1);
	}
	while ( fgets(seqLine, sizeof(seqLine), f_data_sqnc) != NULL ) {
		if ( seqLine[0] != '>' ) {
			for ( int i = 0; i < 50; i ++ ) {
				sequences[seqIdx++] = seqLine[i];
			}
		}
	}
	fclose(f_data_sqnc);
	printf( "Reading the query sequence file finished!\n" );
	fflush( stdout );

	// Read the motif text file then
	int mtfCnt = 0;
	int mtfIdx = 0;
	char mtfLine[1024];
	char *motifs = (char*)malloc(sizeof(char)*MOTIFNUM*MOTIFLENGTH);
	char *filename_motif = "../MA0003.2.fasta";
	FILE* f_data_motif = fopen(filename_motif, "r");
	if ( f_data_motif == NULL ) {
		printf( "File not found: %s\n", filename_motif );
		exit(1);
	}
	while ( fgets(mtfLine, sizeof(mtfLine), f_data_motif) != NULL ) {
		if ( mtfLine[0] != '>' ) {
			if ( mtfCnt == 0 ) {
				for ( int i = 0; i < 60; i ++ ) {
					motifs[mtfIdx++] = mtfLine[i];
				}
				mtfCnt = 1;
			} else {
				for ( int i = 0; i < 55; i ++ ) {
					motifs[mtfIdx++] = mtfLine[i];
				}
				mtfCnt = 0;
			}
		}
	}
	fclose(f_data_motif);
	printf( "Reading the motif file finished!\n" );
	fflush( stdout );

	// Dataset manipulation part
	for ( int i = 0; i < SEQNUM; i ++ ) {
		int pos = rand() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			int motifCnt = SEQNUM % MOTIFNUM;
			sequences[SEQLENGTH*i + pos + j] = motifs[MOTIFLENGTH*motifCnt + j];
		}
	}
	printf( "Manipulating dataset finished!\n" );
	fflush( stdout );

	// Print the results
	char *filename_result = "../TFAP2A.fasta";
	FILE *f_data_result = fopen(filename_result, "w");
	for ( int i = 0; i < SEQNUM; i ++ ) {
		fprintf(f_data_result, ">%d\n", i);
		for ( int j = 0; j < SEQLENGTH; j ++ ) {
			fputc(toupper(sequences[SEQLENGTH*i + j]), f_data_result);
		}
		fputc('\n', f_data_result);
	}
	fclose(f_data_result);
	printf( "Writing the manipulated sequence finished!\n" );
	fflush( stdout );

	// Free
	free(sequences);
	free(motifs);

	return 0;
}
