/*
 * ManageStats.h
 *
 *  Created on: 9 de sept. de 2016
 *      Author: galvez
 */

#ifndef MANAGESTATS_H_
#define MANAGESTATS_H_

void createStats(char *filename);
void checkcreateStats();

extern SequenceDemultiplexed * databaseAlignedDemultiplexed;
extern uint32_t databaseNumSequences;
typedef struct {
    char * filename;
    uint8_t id;
    SequenceDemultiplexed * arraySeqAlgDemul;
    uint32_t first, last; // Both included. 0-indexed
} paramToLoadThread;

void loadFastaAndBalanceLoad(char *filename, int * first, int * last, int numWorkers);
void loadFastaFromStats(char *filename);
void * loadBunchFromStats(void * params);
void checkloadFastaFromStats();
void checkLoadFastaAndBalanceLoad(int * first, int * last, int numWorkers);

#endif /* MANAGESTATS_H_ */
