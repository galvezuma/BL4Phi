/* 
 * File:   Types.h
 * Author: galvez
 *
 * Created on 5 de septiembre de 2020, 18:03
 */

#ifndef TYPES_H
#define	TYPES_H

#include <inttypes.h>
#include <pthread.h>

#define MAX_FILENAME_LENGTH 1024
#define MAX_LINE_LENGTH 1024
#define NUM_THREAD_FOR_LOADING 8
#define MAX_NUM_THREAD_FOR_PROCESSING 12
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

struct GlobalContext {
	int open_gap_cost;
	int extend_gap_cost;
	int threshold; // See http://www.biology.wustl.edu/gcg/psiblast.html
	int numThreadsForProcessing;
	Sequence * sec_ref;
	int8_t * matrix;
	uint32_t num_letters;
	uint8_t letters[128];
	struct {
		char * letras_sec_ref;
		uint16_t long_profile;
		int32_t * profile;
	} Farrar;
	pthread_mutex_t mutex_next_db_seq;
	volatile int next_db_seq_number;
	pthread_mutex_t mutex_check_best;
	volatile uint32_t best_score;
	volatile int best_databaseIdx;
	// Single byte fields
	unsigned char nearby;
	unsigned char nonExhaustive;
	unsigned char bestOnly;
};

extern struct GlobalContext Context;

#endif	/* TYPES_H */

