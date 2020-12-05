/*
 * ManageDatabase.c
 *
 *  Created on: 13 de oct. de 2020
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
uint8_t numLoaders = 0;
LoaderThreadParameters paramsLoad[NUM_THREAD_FOR_LOADING];

/*
 * Because calling malloc is costly, we load at once the database file into memory and
 * we use this memory as an "already allocated" block to where the pointers must point to.
 * The pointers are made to point the correct locations in several threads to go faster. In
 * a second step, the data demultiplexed is allocated in a single _mm_malloc operation and, then
 * data demultiplexed is calculated and pointed to.
 * Each thread returns a linked list of Sequences. However, because allocating a single node is costly,
 * nodes are allocated in blocks of 2000 so the number of mallocs is reduced in a factor of 2000.
 * Finally, the linked lists are rearranged in an array.
 */
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
    numLoaders = (size < 100000)? 1 : NUM_THREAD_FOR_LOADING;
    double ratio = 1.0/(double)numLoaders;
    pthread_t threads[numLoaders];

    for(int i=0; i < numLoaders; i++){
        paramsLoad[i].id = i;
        paramsLoad[i].bulkFile = bulkFile;
        paramsLoad[i].first = size*ratio*i;
        paramsLoad[i].last = size*ratio*(i+1)-1;
        pthread_create(&threads[i], NULL, loadStep1, (void *)&paramsLoad[i]);
    }
    // We need to wait the threads to finish because we have passed local variables as parameters
    databaseNumSequences = 0;
    for(int i=0; i < numLoaders; i++){
        pthread_join(threads[i], NULL);
        databaseNumSequences += paramsLoad[i].ret.numSequences;
    }
    printf("Num Sequences: %d\n", databaseNumSequences);
    // Now, it is time to put all the sequences in a single array
    databaseAlignedDemultiplexed = (SequenceDemultiplexed *)malloc(sizeof(SequenceDemultiplexed)*databaseNumSequences);
    uint64_t total_letters = 0;
    uint32_t currentGlobalPos = 0;
    SequencesChunk * currentChunk = NULL;
    for(int i=0; i < numLoaders; i++){
    	 
    	currentChunk = paramsLoad[i].ret.headerChunk;
        for (int counter = paramsLoad[i].ret.numSequences - 1; counter >= 0 ; counter --){
        	int posInChunk = counter % CHUNK_SEQ_SIZE;
        	databaseAlignedDemultiplexed[currentGlobalPos] = currentChunk->chunkSequences[posInChunk];
        	total_letters += currentChunk->chunkSequences[posInChunk].realDataLength;
        	currentGlobalPos ++;
        	if (counter % CHUNK_SEQ_SIZE == 0) {
				SequencesChunk * aux = currentChunk;
				currentChunk = currentChunk->ptrNext;
				free(aux);
        	}
        }
        if (currentChunk != NULL) { free(currentChunk); currentChunk = NULL; }
        paramsLoad[i].ret.headerChunk = NULL;
    }
    printf("Total aa: %ld\n", total_letters);
    if (currentGlobalPos != databaseNumSequences) errorAndExit("", "Fatal error: Correctness in numSequences not met.");
    //
    gettimeofday(&tiempo_final, NULL);
    printf("Time to load database: %5.3f\n", (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0));
}

void freeDatabase(){
	if (bulkFile != NULL) free(bulkFile);
	for(int i=0; i< numLoaders; i++){
		_mm_free(paramsLoad[i].retBulkData);
	}
}

void * loadStep1(void * vparams){
	LoaderThreadParameters * params = (LoaderThreadParameters *)vparams;
	char *saveptr;
	uint64_t bulkDataSize = 0;
    initSequenceList(&(params->ret));
    char * currentPos;
    char * okReading;
    SequenceDemultiplexed currentSeq;
    //
    // We go to the first occurrence of '>'
    currentPos = params->bulkFile+params->first;
    if (*currentPos == '>'){
    	currentPos = strtok_r(params->bulkFile+params->first, "\r\n", &saveptr);
    } else {
		currentPos = strtok_r(params->bulkFile+params->first, ">", &saveptr);
		if (currentPos == NULL) return NULL;
		if (*currentPos != '>') {
			currentPos = strtok_r(NULL, "\r\n", &saveptr);
			*(--currentPos) = '>';
		}
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
 
    	if (currentSeq.realDataLength < 4) {
    		currentSeq.dataLength = 0;
    	} else if (Context.nonExhaustive){
            // We skip the 2nd and 4th copies shifted
           
    		currentSeq.dataLength = 2 * excess64(currentSeq.realDataLength);
    	} else {
    		
    		currentSeq.dataLength = 4 * excess64(currentSeq.realDataLength);
    	}
    	bulkDataSize += currentSeq.dataLength;
    	addSequence(&(params->ret), &currentSeq);
    }

    // Step 2: Load Data demultiplexed

    params->retBulkData = (char *) _mm_malloc(bulkDataSize, 64);
    char * bulkData = params->retBulkData;
    //
        uint32_t offset = 0;
        SequencesChunk * currentChunk = currentChunk = params->ret.headerChunk;
        for (int counter = params->ret.numSequences - 1; counter >= 0; counter --){
        	int posInChunk = counter % CHUNK_SEQ_SIZE;

        	currentChunk->chunkSequences[posInChunk].data = bulkData + offset;
        	offset += currentChunk->chunkSequences[posInChunk].dataLength;
        	if (currentChunk->chunkSequences[posInChunk].dataLength == 0) continue;
        	// We copy the sequence as explained in the documentation: demultiplexing
    				// Create first copy
    				int posSeq = 0;
    				uint32_t lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength);
    				memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData, lengthToCopy);
    				posSeq += lengthToCopy;
    				uint32_t lengthToZero = excess64(currentChunk->chunkSequences[posInChunk].realDataLength) - lengthToCopy;
    				memset(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), 0, lengthToZero);
    				posSeq += lengthToZero;
    				//
    				if (! Context.nonExhaustive) {
						// Create second copy
						lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength - 1);
						memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData+1, lengthToCopy);
						posSeq += lengthToCopy;
						lengthToZero = excess64(currentChunk->chunkSequences[posInChunk].realDataLength) - lengthToCopy;
	    				memset(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), 0, lengthToZero);
	    				posSeq += lengthToZero;
    				}
    				//
    				// Create third copy
    				lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength - 2);
    				memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData+2, lengthToCopy);
    				posSeq += lengthToCopy;
    				lengthToZero = excess64(currentChunk->chunkSequences[posInChunk].realDataLength) - lengthToCopy;
    				memset(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), 0, lengthToZero);
    				posSeq += lengthToZero;
    				//
    				if (! Context.nonExhaustive) {
    					// Create fourth copy
						lengthToCopy = trunc4(currentChunk->chunkSequences[posInChunk].realDataLength - 3);
						memcpy(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), currentChunk->chunkSequences[posInChunk].realData+3, lengthToCopy);
						posSeq += lengthToCopy;
	    				lengthToZero = excess64(currentChunk->chunkSequences[posInChunk].realDataLength) - lengthToCopy;
	    				memset(&(currentChunk->chunkSequences[posInChunk].data[posSeq]), 0, lengthToZero);
	    				posSeq += lengthToZero;
    				}
    				// Shrink data
    				for(int i=0; i<currentChunk->chunkSequences[posInChunk].dataLength; i++)
    					currentChunk->chunkSequences[posInChunk].data[i] = shrinkLetter(currentChunk->chunkSequences[posInChunk].data[i]);
    				//
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
			if (counter % CHUNK_SEQ_SIZE == 0) {
				currentChunk = currentChunk->ptrNext;
			}
    	}
    return NULL;
}

int interleaveZero3(char * dest, char * orig, int len) {
	if (len % 3 != 0){
		errorAndExit("", "Internal error: inconsistency in 3 module.");
	}
	int posDest = 0, posOrig = 0;
	while (posOrig < len) {
		dest[posDest++] = orig[posOrig++];
		dest[posDest++] = orig[posOrig++];
		dest[posDest++] = orig[posOrig++];
		dest[posDest++] = 0;
	}
	return posDest;
}

void initSequenceList(SequencesList *sl){
	if (sl == NULL) errorAndExit("", "Internal error: NULL data.");
	sl->numSequences = 0;
	sl->headerChunk = NULL;
}

void addSequence(SequencesList *sl, SequenceDemultiplexed * s){
	if (sl->numSequences % CHUNK_SEQ_SIZE == 0) {
		SequencesChunk * aux = (SequencesChunk *) malloc(sizeof(SequencesChunk));
	    if (aux == NULL) errorAndExit("", "Fatal error: Out of memory.");
		aux->ptrNext = sl->headerChunk;
		sl->headerChunk = aux;
	}
	sl->headerChunk->chunkSequences[sl->numSequences % CHUNK_SEQ_SIZE] = *s;
	sl->numSequences++;
}
