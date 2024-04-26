#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <time.h>


// For dataset 1
#define DATASETNUM 32768
#define DATASETLENGTH 1000
// For dataset 2
//#define DATASETNUM 65536
//#define DATASETLENGTH 1000
// For dataset 3
//#define DATASETNUM 131072
//#define DATASETLENGTH 1000


// Main
int main() {
	// Dataset generation part
	char *sequences = (char*)malloc(sizeof(char)*DATASETNUM*DATASETLENGTH);
	for ( int i = 0; i < DATASETNUM; i ++ ) {
		for ( int j = 0; j < DATASETLENGTH; j ++ ) {
			int r = rand() % 20;
			if ( r == 0 ) sequences[DATASETLENGTH*i + j] = 'M';
			else if ( r == 1 ) sequences[DATASETLENGTH*i + j] = 'T';
			else if ( r == 2 ) sequences[DATASETLENGTH*i + j] = 'N';
			else if ( r == 3 ) sequences[DATASETLENGTH*i + j] = 'K';
			else if ( r == 4 ) sequences[DATASETLENGTH*i + j] = 'S';
			else if ( r == 5 ) sequences[DATASETLENGTH*i + j] = 'R';
			else if ( r == 6 ) sequences[DATASETLENGTH*i + j] = 'V';
			else if ( r == 7 ) sequences[DATASETLENGTH*i + j] = 'A';
			else if ( r == 8 ) sequences[DATASETLENGTH*i + j] = 'D';
			else if ( r == 9 ) sequences[DATASETLENGTH*i + j] = 'E';
			else if ( r == 10 ) sequences[DATASETLENGTH*i + j] = 'G';
			else if ( r == 11 ) sequences[DATASETLENGTH*i + j] = 'F';
			else if ( r == 12 ) sequences[DATASETLENGTH*i + j] = 'L';
			else if ( r == 13 ) sequences[DATASETLENGTH*i + j] = 'Y';
			else if ( r == 14 ) sequences[DATASETLENGTH*i + j] = 'C';
			else if ( r == 15 ) sequences[DATASETLENGTH*i + j] = 'W';
			else if ( r == 16 ) sequences[DATASETLENGTH*i + j] = 'P';
			else if ( r == 17 ) sequences[DATASETLENGTH*i + j] = 'H';
			else if ( r == 18 ) sequences[DATASETLENGTH*i + j] = 'Q';
			else if ( r == 19 ) sequences[DATASETLENGTH*i + j] = 'I';
		}
	}
	printf( "Generating dataset finished!\n" );
	fflush( stdout );

	// Print the results
	char *filename_result = "../DATASET_PROTEIN_1.fasta";
	FILE *f_data_result = fopen(filename_result, "w");
	for ( int i = 0; i < DATASETNUM; i ++ ) {
		fprintf(f_data_result, ">%d\n", i);
		for ( int j = 0; j < DATASETLENGTH; j ++ ) {
			fputc(sequences[DATASETLENGTH*i + j], f_data_result);
		}
		fputc('\n', f_data_result);
	}
	fclose(f_data_result);
	printf( "Writing the generated sequence finished!\n" );
	fflush( stdout );

	// Free
	free(sequences);

	return 0;
}
