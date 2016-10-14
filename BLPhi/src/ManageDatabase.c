/*
 * ManageDatabase.c
 *
 *  Created on: 13 de oct. de 2016
 *      Author: galvez
 */

#define __USE_LARGEFILE
#define _POSIX_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <pthread.h>
#include <sys/time.h>
#include <inttypes.h>

#include "ManageDatabase.h"
#include "GeneralFunctions.h"

extern SequenceDemultiplexed * databaseAlignedDemultiplexed;
extern uint32_t databaseNumSequences;

char * bulkFile = NULL;

void loadDatabase(char * filename){
	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);

    FILE * in = fopen(filename, "rb");
    if (in == NULL) errorAndExit(filename, "Cannot read from file.");

    fseek(in, 0, SEEK_END);
    long int size = ftell(in);
    fseek(in, 0, SEEK_SET);
    bulkFile = (char *)malloc(size+1);
    if (bulkFile == NULL) errorAndExit("", "Fatal error: Out of memory.");
    fread(bulkFile, size, 1, in);
    bulkFile[size]=0;
    fclose(in);

    // Create 4 threads to load sequences into memory if needed.
    uint8_t numLoaders = (size < 100000)? 1 : NUM_THREAD_FOR_LOADING;
    double ratio = 1.0/(double)numLoaders;
    pthread_t threads[numLoaders];
    LoaderThreadParameters params[numLoaders];
    for(int i=0; i < numLoaders; i++){
        params[i].id = i;
        params[i].bulkFile = bulkFile;
        params[i].first = size*ratio*i;
        params[i].last = size*ratio*(i+1)-1;
        pthread_create(&threads[i], NULL, loadStep1, (void *)&params[i]);
    }
    // We need to wait the threads to finish because we have passed local variables as parameters
    databaseNumSequences = 0;
    for(int i=0; i < numLoaders; i++){
        pthread_join(threads[i], NULL);
        databaseNumSequences += params[i].ret.numSequences;
    }
    printf("Num Sequences: %d\n", databaseNumSequences);
    // Now, it is time to put all the sequences in a single array
    databaseAlignedDemultiplexed = (SequenceDemultiplexed *)malloc(sizeof(SequenceDemultiplexed)*databaseNumSequences);
    uint32_t currentGlobalPos = 0;
    SequencesChunk * currentChunk = NULL;
    for(int i=0; i < numLoaders; i++){
    	// printf("Grouping %d\n", i);
        for (int counter = 0; counter < params[i].ret.numSequences; counter ++){
        	if (counter % CHUNK_SEQ_SIZE == 0) {
        		if (currentChunk == NULL) currentChunk = params[i].ret.headerChunk;
        		else {
        			SequencesChunk * aux = currentChunk;
        			currentChunk = currentChunk->ptrNext;
        			free(aux);
        		}
        	}
        	int posInChunk = counter % CHUNK_SEQ_SIZE;
        	databaseAlignedDemultiplexed[currentGlobalPos] = currentChunk->chunkSequences[posInChunk];
        	currentGlobalPos ++;
        }
        if (currentChunk != NULL) { free(currentChunk); currentChunk = NULL; }
        params[i].ret.headerChunk = NULL;
    }
    if (currentGlobalPos != databaseNumSequences) errorAndExit("", "Fatal error: Correctness in numSequences not met.");
    //
    gettimeofday(&tiempo_final, NULL);
    printf("Time to load database: %5.3f\n", (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0));
}

void * loadStep1(void * vparams){
	LoaderThreadParameters * params = (LoaderThreadParameters *)vparams;
	char *saveptr;
	uint64_t bulkDataSize = 0;
    initSequenceList(&(params->ret));

    char * currentPos;
    char sequence[MAX_SEQUENCE_LENGTH];
    char * okReading;
    SequenceDemultiplexed currentSeq;
    //
    // We go to the first occurrence of '>'
    currentPos = strtok_r(params->bulkFile+params->first, ">", &saveptr);
    if (currentPos == NULL) return NULL;
    if (*currentPos != '>') {
    	currentPos = strtok_r(NULL, "\r\n", &saveptr);
    	*(--currentPos) = '>';
    }
    while (currentPos != NULL && currentPos <= params->bulkFile + params->last) { // We are beyond our competence
    	currentSeq.name = currentPos;
    	currentPos = strtok_r(NULL, "\r\n", &saveptr);
    	currentSeq.realData = currentPos; // This points to the first line of the sequence
    	currentSeq.realDataLength = strlen(currentSeq.realData);
    	currentPos = strtok_r(NULL, "\r\n", &saveptr);
    	while(currentPos != NULL && *currentPos != '>'){
    		int auxLen = strlen(currentPos);
    		memcpy(currentSeq.realData+currentSeq.realDataLength, currentPos, auxLen);
    		currentSeq.realDataLength += auxLen;
        	currentPos = strtok_r(NULL, "\r\n", &saveptr);
    	}
//    	printf("%s\n", currentSeq.name);
//    	printf("%d --- %s\n", currentSeq.realDataLength, currentSeq.realData);
    	if (currentSeq.realDataLength < 4) {
    		currentSeq.dataLength = 0;
    	} else {
    		currentSeq.dataLength = excess64(4*currentSeq.realDataLength-12);
    	}
    	bulkDataSize += currentSeq.dataLength;
    	addSequence(&(params->ret), &currentSeq);
    }

    // Step 2: Load Data demultiplexed

    params->retBulkData = (char *) _mm_malloc(bulkDataSize, 64);
    char * bulkData = params->retBulkData;
    //
        uint32_t offset = 0;
        SequencesChunk * currentChunk = NULL;
        for (int counter = 0; counter < params->ret.numSequences; counter ++){
        	if (counter % CHUNK_SEQ_SIZE == 0) {
        		if (currentChunk == NULL) currentChunk = params->ret.headerChunk;
        		else currentChunk = currentChunk->ptrNext;
        	}
        	int posInChunk = counter % CHUNK_SEQ_SIZE;

        	currentChunk->chunkSequences[posInChunk].data = bulkData + offset;
        	offset += currentChunk->chunkSequences[posInChunk].dataLength;
        	if (currentChunk->chunkSequences[posInChunk].dataLength == 0) continue;
        	// We copy the sequence as explained in the doumentation: demultiplexing
    		// Create first copy
    				int posSeq = 0;
    				uint32_t lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength);
    				memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData, lengthToCopy);
    				// Create second copy
    				posSeq += lengthToCopy;
    				lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength - 1);
    				memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData+1, lengthToCopy);
    				// Create third copy
    				posSeq += lengthToCopy;
    				lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength - 2);
    				memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData+2, lengthToCopy);
    				// Create fourth copy
    				posSeq += lengthToCopy;
    				lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength - 3);
    				memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData+3, lengthToCopy);
    				// Fill the final gaps
    				posSeq += lengthToCopy;
    				while (posSeq%64 != 0) currentChunk->chunkSequences[posInChunk].data[posSeq++] = 0;
    				if (currentChunk->chunkSequences[posInChunk].dataLength != posSeq)  {
    					printf("Seq of length %d, name %s\n%s\n%d\n%d\n\n",
    							currentChunk->chunkSequences[posInChunk].realDataLength,
								currentChunk->chunkSequences[posInChunk].name,
								currentChunk->chunkSequences[posInChunk].realData,
    							posSeq,
								currentChunk->chunkSequences[posInChunk].dataLength);
    					errorAndExit("", "Internal error: inconsistency in alignments.");
    				}
            // Finish of demultiplexing
    	}
    return NULL;
}

void initSequenceList(SequencesList *sl){
	if (sl == NULL) errorAndExit("", "Internal error: NULL data.");
	sl->numSequences = 0;
	sl->headerChunk = NULL;
}

void * addSequence(SequencesList *sl, SequenceDemultiplexed * s){
	if (sl == NULL || s == NULL) errorAndExit("", "Internal error: NULL data.");
	if (sl->numSequences % CHUNK_SEQ_SIZE == 0) {
		SequencesChunk * aux = (SequencesChunk *) malloc(sizeof(SequencesChunk));
	    if (aux == NULL) errorAndExit("", "Fatal error: Out of memory.");
		aux->ptrNext = sl->headerChunk;
		sl->headerChunk = aux;
	}
	sl->headerChunk->chunkSequences[sl->numSequences % CHUNK_SEQ_SIZE] = *s;
	sl->numSequences++;
	return (void *) sl->headerChunk;
}
