/*
 * Threads.c
 *
 *  Created on: Dec 3, 2020
 *      Author: galvez
 */

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "Threads.h"
#include "Types.h"
#include "SingleQuery.h"


extern uint32_t databaseNumSequences;

typedef struct {
	pthread_mutex_t semaphoreForNewQuery;
	pthread_mutex_t semaphoreForThreadFinishQuery;
	uint8_t stop;
	paramToSingleQueryProcessThread params;
} SemaphoresOfThreads;

uint8_t numThreads = 0;
pthread_t *threads = NULL;
SemaphoresOfThreads * semaphores = NULL;

void * runThread(void * vparams);

uint8_t numberOfThreads() {
	return numThreads;
}

void allocateThreads() {
    // Create threads to process sequences.
    numThreads = (databaseNumSequences < 10 * Context.numThreadsForProcessing)? 1 : Context.numThreadsForProcessing;
    threads = (pthread_t *) malloc(numThreads * sizeof(pthread_t));
    semaphores = (SemaphoresOfThreads *) malloc(numThreads * sizeof(SemaphoresOfThreads));
}

void deallocateThreads() {
	free(threads);
	threads = NULL;
	free(semaphores);
	semaphores = NULL;
}

void startAndWaitThreads() {
    for(int i=0; i < numThreads; i++){
    	// Initialize data of threads
    	pthread_mutex_init(& (semaphores[i].semaphoreForNewQuery), NULL);
    	pthread_mutex_init(& (semaphores[i].semaphoreForThreadFinishQuery), NULL);
    	pthread_mutex_lock(& (semaphores[i].semaphoreForNewQuery));
    	pthread_mutex_lock(& (semaphores[i].semaphoreForThreadFinishQuery));
    	semaphores[i].stop = 0;
    	// Initialize some params
    	semaphores[i].params.id = i;
    	pthread_create(&(threads[i]), NULL, runThread, (void *)&semaphores[i]);
    }

}
void restartThreads(Sequence * query) {
    for(int i=0; i < numThreads; i++){
    	semaphores[i].params.query = query;
    	pthread_mutex_unlock(& (semaphores[i].semaphoreForNewQuery));
    }
}

paramToSingleQueryProcessThread * waitForThread(int id) {
	pthread_mutex_lock(& (semaphores[id].semaphoreForThreadFinishQuery));
	return &(semaphores[id].params);
}

void stopThreads() {
    for(int i=0; i < numThreads; i++){
    	semaphores[i].stop = 1;
    	pthread_mutex_unlock(& (semaphores[i].semaphoreForNewQuery));
    }
}


void * runThread(void * vparams) {
	SemaphoresOfThreads * mySemaphores = (SemaphoresOfThreads *) vparams;
	while (! mySemaphores->stop) {
		// Wait for a new query to be ready
		pthread_mutex_lock(&(mySemaphores->semaphoreForNewQuery));
		if (mySemaphores->stop) continue;
		// Call the main function of a Thread.
		processBunchSingleFastaWholeDatabase(&(mySemaphores->params));
		pthread_mutex_unlock(&(mySemaphores->semaphoreForThreadFinishQuery));
	}
	pthread_exit(NULL);
}
