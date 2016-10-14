/* 
 * File:   Types.h
 * Author: galvez
 *
 * Created on 5 de septiembre de 2016, 18:03
 */

#ifndef TYPES_H
#define	TYPES_H

#include <inttypes.h>

#define MAX_FILENAME_LENGTH 1024
#define MAX_LINE_LENGTH 1024
#define NUM_THREAD_FOR_LOADING 16
#define NUM_THREAD_FOR_PROCESSING 228
#define MAX_SEQUENCE_LENGTH (1024*100)
#define VECTOR_SIZE 64

typedef struct {
	int numSequences;
	uint32_t fastaHeaderNumBytes;
	uint32_t fastaRealDataNumBytes;
	uint64_t fastaDataNumBytes;
} FileHeader;

typedef struct {
    char * name;
    char * data;
    uint32_t dataLength;
} Sequence;

typedef struct {
    char * name;
    char * data;
    uint32_t dataLength;
    char * realData;
    uint32_t realDataLength;
} SequenceDemultiplexed;



#endif	/* TYPES_H */

