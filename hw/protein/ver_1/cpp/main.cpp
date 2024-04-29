#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>

#include "bdbmpcie.h"
#include "dmasplitter.h"


// Parameter setting
#define SEQNUM 32768
#define SEQLENGTH 256
#define MOTIFLENGTH 16


// Elapsed time checker
double timespec_diff_sec( timespec start, timespec end ) {
	double t = end.tv_sec - start.tv_sec;
	t += ((double)(end.tv_nsec - start.tv_nsec)/1000000000L);
	return t;
}

// Function for representing sequences as 5-bit
void encoderSeq( uint64_t* sequencesEncoded, char* sequences ) {
	uint32_t idx = 0;
	uint64_t unit = 0;
	uint64_t tmp = 0;
	uint64_t cnt_1 = 0;
	uint64_t cnt_2 = 0;
	for ( uint32_t i = 0; i < SEQNUM; i ++ ) {
		for ( uint32_t j = 0; j < SEQLENGTH; j ++ ) {
			if ( sequences[SEQLENGTH*i + j] == 'M' ) unit = (0 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'T' ) unit = (1 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'N' ) unit = (2 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'K' ) unit = (3 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'S' ) unit = (4 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'R' ) unit = (5 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'V' ) unit = (6 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'A' ) unit = (7 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'D' ) unit = (8 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'E' ) unit = (9 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'G' ) unit = (10 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'F' ) unit = (11 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'L' ) unit = (12 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'Y' ) unit = (13 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'C' ) unit = (14 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'W' ) unit = (15 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'P' ) unit = (16 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'H' ) unit = (17 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'Q' ) unit = (18 << 5*cnt_2) | unit;
			else if ( sequences[SEQLENGTH*i + j] == 'I' ) unit = (19 << 5*cnt_2) | unit;

			if ( cnt_1 == 0 ) {
				if ( cnt_2 + 1 == 12 ) {
					tmp_1 = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 1 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 60) | tmp;
					unit = unit >> 4;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					tmp = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 2 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 61) | tmp;
					unit = unit >> 3;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					tmp = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 3 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 62) | tmp;
					unit = unit >> 2;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					tmp = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 4 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 63) | tmp;
					unit = unit >> 1;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					sequencesEncoded[idx++] = unit;
					unit = 0;
					cnt_1 = 0;
					cnt_2 = 0;
				} else cnt_2++;
			}
		}
	}
}

// Function for representing motifs as 5-bit
void encoderMtf( uint64_t* motifsEncoded, char* motifs ) {
	uint32_t idx = 0;
	uint64_t unit = 0;
	uint64_t tmp = 0;
	uint64_t cnt_1 = 0;
	uint64_t cnt_2 = 0;
	for ( uint32_t i = 0; i < SEQNUM/16; i ++ ) {
		for ( uint32_t j = 0; j < MOTIFLENGTH*16; j ++ ) {
			if ( motifs[MOTIFLENGTH*16*i + j] == 'M' ) unit = (0 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'T' ) unit = (1 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'N' ) unit = (2 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'K' ) unit = (3 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'S' ) unit = (4 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'R' ) unit = (5 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'V' ) unit = (6 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'A' ) unit = (7 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'D' ) unit = (8 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'E' ) unit = (9 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'G' ) unit = (10 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'F' ) unit = (11 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'L' ) unit = (12 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'Y' ) unit = (13 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'C' ) unit = (14 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'W' ) unit = (15 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'P' ) unit = (16 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'H' ) unit = (17 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'Q' ) unit = (18 << 5*cnt_2) | unit;
			else if ( motifs[MOTIFLENGTH*16*i + j] == 'I' ) unit = (19 << 5*cnt_2) | unit;

			if ( cnt_1 == 0 ) {
				if ( cnt_2 + 1 == 12 ) {
					tmp_1 = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 1 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 60) | tmp;
					unit = unit >> 4;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					tmp = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 2 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 61) | tmp;
					unit = unit >> 3;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					tmp = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 3 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 62) | tmp;
					unit = unit >> 2;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					tmp = unit;
					unit = 0;
					cnt_1++;
					cnt_2 = 0;
				} else cnt_2++;
			} else if ( cnt_1 == 4 ) {
				if ( cnt_2 == 0 ) {
					sequencesEncoded[idx++] = (unit << 63) | tmp;
					unit = unit >> 1;
					cnt_2++;
				} else if ( cnt_2 + 1 == 13 ) {
					sequencesEncoded[idx++] = unit;
					unit = 0;
					cnt_1 = 0;
					cnt_2 = 0;
				} else cnt_2++;
			}
		}
	}
}


// Main
int main(int argc, char** argv) {
	BdbmPcie* pcie = BdbmPcie::getInstance();
	uint64_t* dmabuf = (uint64_t*)pcie->dmaBuffer();

	unsigned int d = pcie->readWord(0);
	printf( "Magic: %x\n", d );
	fflush( stdout );
	//-------------------------------------------------------------------------------
	// Stage1: Read a file that contains the PROTEIN sequences
	//-------------------------------------------------------------------------------
	uint32_t seqIdx = 0;
	char seqLine[1024];
	char *sequences = (char*)malloc(sizeof(char)*SEQNUM*SEQLENGTH);
	char *filename_sequences = "../../../dataset/DATASET_PROTEIN_1.fasta";
	FILE *f_data_sequences = fopen(filename_sequences, "r");
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
	printf( "[SW]: Reading the sequence file finished!\n" );
	fflush( stdout );
	//-------------------------------------------------------------------------------
	// Stage2: Convert char to 5bits (encoding) 
	// 1 sequence = 256 x 5bits = 64 x 4 x 5 = 128 x 2 x 5
	//-------------------------------------------------------------------------------
	uint32_t sizeSeq = SEQNUM*20;
	uint64_t *sequencesEncoded = (uint64_t*)malloc(sizeof(uint64_t)*sizeSeq);
	encoderSeq(&sequencesEncoded[0], &sequences[0]);
	printf( "[SW]: Encoding the sequences finished!\n" );
	fflush( stdout );
	//-------------------------------------------------------------------------------
	// Stage3: Designate motif of each sequence randomly
	//-------------------------------------------------------------------------------
	char *motifs = (char*)malloc(sizeof(char)*SEQNUM*MOTIFLENGTH);
	for ( int i = 0; i < SEQNUM; i ++ ) {
		int sp = rand() % (SEQLENGTH - MOTIFLENGTH + 1);
		for ( int j = 0; j < MOTIFLENGTH; j ++ ) {
			motifs[MOTIFLENGTH*i + j] = sequences[SEQLENGTH*i + sp + j];
		}
	}
	printf( "[SW]: Designating the motifs randomly finished!\n" );
	fflush( stdout );
	//-------------------------------------------------------------------------------
	// Stage4: Do 5-bit encoding for the motifs
	// 16 motifs = 16 x 16 x 5bits = 64 x 4 x 5 = 128 x 2 x 5
	//-------------------------------------------------------------------------------
	uint32_t sizeMtf = (SEQNUM/16)*20;
	uint64_t *motifsEncoded = (uint64_t*)malloc(sizeof(uint64_t)*sizeMtf);
	encoderMtf(&motifsEncoded[0], &motifs[0]);
	printf( "[SW]: Encoding the motifs finished!\n" );
	fflush( stdout );
	//-------------------------------------------------------------------------------
	// Stage5: Send the sequences to FPGA through DMA
	// SEQNUM x 20 x 8B(64bits) = SEQNUM x 10 x 16B(128bits)
	// 640 x 8KB
	//-------------------------------------------------------------------------------
	uint32_t idx = 0;
	while ( idx < sizeSeq ) {
		// 8KB
		for ( uint32_t i = 0; i < 1024; i ++ ) {
			dmabuf[i] = sequencesEncoded[idx++];

			if ( idx == sizeSeq ) break;
		}
		if ( (idx == sizeSeq) && (sizeSeq % (512*2) != 0) ) {
			uint32_t left = idx % (512*2);
			if ( left % 2 == 0 ) pcie->userWriteWord(0, (left/2));
			else pcie->userWriteWord(0, (left/2)+1);
		} else pcie->userWriteWord(0, 512);
	}
	printf( "[SW]: Sending the encoded sequences to DMA buffer finished!\n" );
	fflush( stdout );
	//-------------------------------------------------------------------------------
	// Stage6: Send the motifs to FPGA through DMA
	// (SEQNUM/16) x 20 x 8B(64bits) = (SEQNUM/16) x 10 x 16B(128bits)
	// 40 x 8KB
	//-------------------------------------------------------------------------------
	idx = 0;
	while ( idx < sizeMtf ) {
		// 8KB
		for ( uint32_t i = 0; i < 1024; i ++ ) {
			dmabuf[i] = motifsEncoded[idx++];

			if ( idx == sizeMtf ) break;
		}

		if ( (idx == sizeMtf) && (sizeMtf % (512*2) != 0) ) {
			uint32_t left = idx % (512*2);
			if ( left % 2 == 0 ) pcie->userWriteWord(0, (left/2));
			else pcie->userWriteWord(0, (left/2)+1);
		} else pcie->userWriteWord(0, 512);
	}
	printf( "[SW]: Sending the encoded motifs finished!\n" );
	fflush( stdout );

	for ( uint32_t i = 0; i < 2*1024*1024; i ++ ) {
		pcie->userWriteWord(8, 0);
	}
	if ( pcie->userReadWord(0) == 1 ) {
		printf( "[SW]: Reading DMA buffer on HW-side successfully finished!\n" );
		fflush( stdout );
	}
	
	for ( uint32_t i = 0; i < 32*1024*1024; i ++ ) {
		pcie->userWriteWord(8, 0);
	}

	free(sequences);
	free(sequencesEncoded);
	free(motifs);
	free(motifsEncoded);

	return 0;
}
