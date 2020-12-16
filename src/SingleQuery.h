/* 
 * File:   SingleQuery.h
 * Author: galvez
 *
 * Created on 5 de septiembre de 2020, 17:56
 */

#ifndef SINGLEQUERY_H
#define	SINGLEQUERY_H

#include "Types.h"

typedef struct {
    uint16_t id;
    Sequence * query;
    //uint32_t first, last; // Both included. 0-indexed
    struct {
    	int32_t numSequencesProcessed;
    	int32_t numLettersProcessed;
    	int32_t numHits;
    	int32_t numFarrar;
    	double total_time;
    } ret;
} paramToSingleQueryProcessThread;

void checkloadSingleFasta(Sequence * query);
void processSingleFastaWholeDatabase(Sequence * query, int * first, int * last, int numWorkers);
void processBunchSingleFastaWholeDatabase(void * params);
void freeSingleFasta(Sequence * query);

#endif	/* SINGLEQUERY_H */

