#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>

// Parameter setting
#define SEQSNUM 26
#define SEQLENGTH 111
#define MOTIFLENGTH 11


// Elapsed time checker
double timespec_diff_sec( timespec start, timespec end ) {
	double t = end.tv_sec - start.tv_sec;
	t += ((double)(end.tv_nsec - start.tv_nsec)/1000000000L);
	return t;
}

// Function for reading benchmark file
void readfromfile( char* data, char* filename, size_t length ) {
	FILE* f_data = fopen(filename, "rb");
	if (f_data == NULL ) {
		printf("File not found: %s\n", filename);
		exit(1);
	}

	fread(data, sizeof(char), length, f_data);

	fclose(f_data);
}

// Deserializer
void deserializer( uint32_t* deserialSeq, char* sequences, size_t length ) {
	// Decide deserializer's length
	size_t deserialDiv = length / 4;
	size_t deserialRem = length % 4;
	
	// Do deserialize
	size_t cnt = 0;
	size_t deserialCnt = 0;
	size_t deserialLen = 0;
	uint32_t fourChar = 0;
	for ( size_t i = 0; i < deserialDiv*4; i ++ ) {
		uint32_t oneChar = (uint32_t)sequences[i];
		fourChar = fourChar | (oneChar << cnt*8);
		if ( cnt == 3 ) {
			deserialSeq[deserialLen++] = fourChar;
			fourChar = 0;
			cnt = 0;
		} else cnt++;
	}
	for ( size_t i = deserialDiv*4; i < length; i ++ ) {
		uint32_t oneChar = (uint32_t)sequences[i];
		fourChar = fourChar | (oneChar << cnt*8);
		if ( i == length - 1 ) {
			deserialSeq[deserialLen] = fourChar;
		} else cnt++;
	}
}

// Serializer
void serializer( char* serialSeq, uint32_t* deserialSeq, size_t length ) {
	size_t deserialLen = 0;
	size_t cnt = 0;
	char oneChar;
	for ( size_t i = 0; i < length; i ++ ) {
		if ( cnt == 0 ) {
			uint32_t oneByte_1 = deserialSeq[deserialLen] << 24;
			uint32_t oneByte_2 = oneByte_1 >> 24;
			oneChar = (char)oneByte_2;
			cnt++;
		} else if ( cnt == 1 ) {
			uint32_t oneByte_1 = deserialSeq[deserialLen] << 16;
			uint32_t oneByte_2 = oneByte_1 >> 24;
			oneChar = (char)oneByte_2;
			cnt++;
		} else if ( cnt == 2 ) {
			uint32_t oneByte_1 = deserialSeq[deserialLen] << 8;
			uint32_t oneByte_2 = oneByte_1 >> 24;
			oneChar = (char)oneByte_2;
			cnt++;
		} else if ( cnt == 3 ) {
			uint32_t oneByte = deserialSeq[deserialLen++] >> 24;
			oneChar = (char)oneByte;
			cnt = 0;
		}
		
		serialSeq[i] = oneChar;
	}
}


// Main
int main(int argc, char** argv) {
	// Read query file
	size_t sizeQuery = SEQSNUM*SEQLENGTH;
	char* filenameQuery = "./dataset/mm9Gata4MotifCollection.bin";
	char* sequences = (char*)malloc(sizeof(char)*sizeQuery);
	readfromfile(&sequences[0], filenameQuery, sizeQuery);

	// Read answer file
	size_t sizeAnswer = SEQSNUM*MOTIFLENGTH;
	char* filenameAnswer = "./dataset/mm9Gata4Solutions.bin";
	char* answer = (char*)malloc(sizeof(char)*sizeAnswer);
	readfromfile(&answer[0], filenameAnswer, sizeAnswer);

	// Deserializing to send 4 characters at once
	size_t deserialDiv = sizeQuery / 4;
	size_t deserialRem = sizeQuery % 4;
	size_t deserialLen = deserialDiv;
	if ( deserialRem != 0 ) deserialLen++;
	uint32_t* deserialSeq = (uint32_t*)malloc(sizeof(uint32_t)*deserialLen);
	deserializer(&deserialSeq[0], &sequences[0], sizeQuery);

	// Serializing to decompose one 32-bit to 4 characters
	char* serialSeq = (char*)malloc(sizeof(char)*sizeQuery);
	serializer(&serialSeq[0], &deserialSeq[0], sizeQuery);

	// Validation
	for ( size_t i = 0; i < SEQSNUM; i ++ ) {
		for ( size_t j = 0; j < SEQLENGTH; j ++ ) {
			printf( "%c", serialSeq[SEQLENGTH*i + j] );
		}
		printf( "\n" );
	}

	return 0;
}
