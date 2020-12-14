#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <malloc.h>
#include <limits.h>
#include <sys/time.h>
#include <immintrin.h>

#include "Types.h"
#include "GeneralFunctions.h"
#include "SingleQuery.h"
#include "Farrar.h"

extern SequenceDemultiplexed * databaseAlignedDemultiplexed;
extern uint32_t databaseNumSequences;


// The last three parameters may be NULL.
// If NULL then the load is balanced based on the number of sequences
// If not NULL then they contain the size of each work (based on the number of bases).
void processSingleFastaWholeDatabase(Sequence * query, int * first, int * last, int numWorkers){
	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);
    // Create threads to process sequences.
    uint8_t numThreads = (databaseNumSequences < 10 * Context.numThreadsForProcessing)? 1 : Context.numThreadsForProcessing;
 
    double ratio = 1.0/(double)numThreads;
    pthread_t threads[numThreads];
    paramToSingleQueryProcessThread params[numThreads];
    pthread_mutex_init(& Context.mutex_next_db_seq, NULL);
    if (Context.bestOnly) {
    	pthread_mutex_init(& Context.mutex_check_best, NULL);
    	Context.best_score = 0;
    	Context.best_databaseIdx = -1;
    }
    Context.next_db_seq_number = 0;
    fprintf(stdout, "\nResults for: %s\n", query->name);
    for(int i=0; i < numThreads; i++){
        params[i].id = i;
        params[i].query = query;
 
        pthread_create(&threads[i], NULL, processBunchSingleFastaWholeDatabase, (void *)&params[i]);
    }
    // We need to wait the threads to finish because we have passed local variables as parameters
    pthread_join(threads[0], NULL);
    for(int i=1; i < numThreads; i++){
        pthread_join(threads[i], NULL);
        // We use index 0 position to accumulate statistical data
        params[0].ret.numFarrar += params[i].ret.numFarrar;
        params[0].ret.numHits +=  params[i].ret.numHits;
        params[0].ret.numLettersProcessed += params[i].ret.numLettersProcessed;
        params[0].ret.numSequencesProcessed += params[i].ret.numSequencesProcessed;
        params[0].ret.total_time = (params[0].ret.total_time < params[i].ret.total_time) ? params[0].ret.total_time : params[i].ret.total_time;
    }
    if (Context.bestOnly) pthread_mutex_destroy(& Context.mutex_check_best);
    pthread_mutex_destroy(& Context.mutex_next_db_seq);
	gettimeofday(&tiempo_final, NULL);
	double tiempo_total = (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0);
	if (Context.bestOnly)
		if (Context.best_databaseIdx != -1)
			fprintf(stdout, "Hit. Score: %d. Sequence: %s\n", Context.best_score, databaseAlignedDemultiplexed[Context.best_databaseIdx].name);
		else
			fprintf(stdout, "No hits found.\n");
	fprintf(stdout, "Time elapsed: %9.1fms. \n", tiempo_total);
	fprintf(stdout, "GCUPS: %5.3f. \n", ((double)params[0].ret.numLettersProcessed * (double)query->dataLength)/(tiempo_total/1000.0)/1000000000.0);
	fprintf(stdout, "Farrar executions %d.\n", params[0].ret.numFarrar);
	fprintf(stdout, "Hits found %d.\n", params[0].ret.numHits);
	fprintf(stdout, "Amino acids processed %d.\n", params[0].ret.numLettersProcessed);
	fprintf(stdout, "Query length %d.\n", query->dataLength);
	fprintf(stdout, "Proteins processed %d.\n", params[0].ret.numSequencesProcessed);
	fprintf(stdout, "Minimum execution time: %9.1fms.\n", params[0].ret.total_time);
}


void * processBunchSingleFastaWholeDatabase(void * vparams) {
	int databaseIdx, startIdx;
	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);
    paramToSingleQueryProcessThread * params = (paramToSingleQueryProcessThread *)vparams;
    FarrarObject o;
    prepareFarrarObject(&o);
    uint16_t neighbourhood[MAX_SEQUENCE_LENGTH];
    // Shrink data
    char shrinkedQuery[params->query->dataLength];
	for(int i=0; i<params->query->dataLength; i++)
		shrinkedQuery[i] = shrinkLetter(params->query->data[i]);
	//

	params->ret.numSequencesProcessed = 0;
	params->ret.numLettersProcessed = 0;
	params->ret.numHits = 0;
	params->ret.numFarrar = 0;
	params->ret.total_time = 0.0;
	//
	const int step = 10;
	while(1) {
		//
				pthread_mutex_lock(& Context.mutex_next_db_seq);
					startIdx = Context.next_db_seq_number;
					Context.next_db_seq_number += step;
				pthread_mutex_unlock(& Context.mutex_next_db_seq);
				if (startIdx >= databaseNumSequences) break;
		//
		for(databaseIdx=startIdx; databaseIdx < startIdx+step; databaseIdx++){
			if (databaseIdx >= databaseNumSequences) break;
			params->ret.numSequencesProcessed ++;
			params->ret.numLettersProcessed += databaseAlignedDemultiplexed[databaseIdx].realDataLength;

			uint16_t nearbyShifter = 0;
			// For each block of 4 consecutive letters in the query
			for(int queryIdx=0; queryIdx < params->query->dataLength - 3; queryIdx++){
				uint32_t block __attribute__((aligned(4))) = * (uint32_t *)&(shrinkedQuery [queryIdx]); // & 0x00FFFFFF;
				__m512i queryVector = _mm512_extload_epi32 (&block, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 1);
				__mmask16 res = 0;
				#pragma unroll(4)
				for(int targetIdx=0; targetIdx < databaseAlignedDemultiplexed[databaseIdx].dataLength ; targetIdx+=VECTOR_SIZE){
	 
					__m512i targetVector = _mm512_load_epi32 ((__m512i const*) (databaseAlignedDemultiplexed[databaseIdx].data + targetIdx));
	//                view512iAsChar(targetVector);
					res |= _mm512_cmpeq_epi32_mask  (queryVector, targetVector);
 
				}
				if (res != 0) {

					if ((nearbyShifter != 0) && (__builtin_popcount (nearbyShifter) >= Context.nearby)) {
							params->ret.numFarrar ++;
							uint32_t score = smith_waterman_farrar(&o, databaseAlignedDemultiplexed[databaseIdx].realData, (int16_t)(databaseAlignedDemultiplexed[databaseIdx].realDataLength));
							if (score >= Context.threshold) {
								params->ret.numHits ++;
								if (! Context.bestOnly) {
									// Show hit with sequence name
									fprintf(stdout, "Hit. Score: %d. Sequence: %s\n", score, databaseAlignedDemultiplexed[databaseIdx].name);
								} else {
									pthread_mutex_lock(& Context.mutex_check_best);
										if (   (score > Context.best_score) // The obtained score is better OR
										    || (  (score == Context.best_score) // The score is the same but the length of the sequence is more similar
										    	&& (abs(params->query->dataLength - databaseAlignedDemultiplexed[databaseIdx].realDataLength)
										    		<
													abs(params->query->dataLength - databaseAlignedDemultiplexed[Context.best_databaseIdx].realDataLength))
											   )
										   ) {
											Context.best_score = score;
											Context.best_databaseIdx = databaseIdx;
										}
									pthread_mutex_unlock(& Context.mutex_check_best);
								}
							}
							queryIdx = INT_MAX - 1;
//                    	}
					}
					nearbyShifter = (nearbyShifter | 0x01);
				}
				nearbyShifter = nearbyShifter << 1;
			}
		}
    }
    freeFarrarObject(&o);
    gettimeofday(&tiempo_final, NULL);
    params->ret.total_time = (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0);
}

void freeSingleFasta(Sequence * query) {
    free(query->name);
    free(query->data);
    free(query);
}

void checkloadSingleFasta(Sequence * query){
    printf("%s\n", query->name);
    printf("%d\n", query->dataLength);
    printf("%s\n", query->data);
}
