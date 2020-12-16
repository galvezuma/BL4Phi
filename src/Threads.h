/*
 * Threads.h
 *
 *  Created on: Dec 3, 2020
 *      Author: galvez
 */

#ifndef THREADS_H_
#define THREADS_H_


#include "Types.h"
#include "SingleQuery.h"

uint8_t numberOfThreads();
void allocateThreads();
void deallocateThreads();
void startAndWaitThreads();
void restartThreads(Sequence * query);
paramToSingleQueryProcessThread * waitForThread(int id);
void stopThreads();

#endif /* THREADS_H_ */
