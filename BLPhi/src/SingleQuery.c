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
extern struct {
	int open_gap_cost;
	int extend_gap_cost;
	int threshold; // See http://www.biology.wustl.edu/gcg/psiblast.html
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
} Context;

Sequence * loadSingleFasta(char * filename) {
    Sequence * seqReturned;
    FILE * in = fopen(filename, "rb");
    if (in == NULL) errorAndExit(filename, "Cannot open file.");
    
    char line[MAX_LINE_LENGTH];
    char sequence[MAX_SEQUENCE_LENGTH];
    char * okReading;
    
    // Reading the name of the sequence
    okReading = fgets(line, MAX_LINE_LENGTH, in);
    if ((okReading == NULL) || (line[0] != '>')) errorAndExit(line, "Not a fasta header.");
    
    seqReturned = (Sequence *)malloc(sizeof(Sequence));
    if (seqReturned == NULL) errorAndExit("SingleFasta", "Fatal error: Out of memory.");
    line[strcspn(line, "\r\n")] = 0;
    uint32_t lengthText = strlen(line);
    seqReturned->name = (char *) malloc(lengthText + 1);
    if (seqReturned->name == NULL) errorAndExit("SingleFasta.name", "Fatal error: Out of memory.");
    strcpy(seqReturned->name, line);
    
    // Reading data of the sequence
    uint32_t posSeq = 0;
    okReading = fgets(line, MAX_LINE_LENGTH, in);
    while((line[0] != '>') && okReading != NULL){
        line[strcspn(line, "\r\n")] = 0;
        lengthText = strlen(line);
        if (posSeq+lengthText > MAX_SEQUENCE_LENGTH) errorAndExit(seqReturned->name, "Fatal error. Sequence too long.");
        strcpy(sequence+posSeq, line);
        posSeq += lengthText;
        okReading = fgets(line, MAX_LINE_LENGTH, in);
    }
    seqReturned->dataLength = posSeq;
    seqReturned->data = (char *) malloc(seqReturned->dataLength);
    if (seqReturned->data == NULL) errorAndExit("SingleFasta.data", "Fatal error: Out of memory.");
    memcpy(seqReturned->data, sequence, seqReturned->dataLength);
    if (seqReturned->dataLength < 4) errorAndExit("SingleFasta.length", "The minimum length is 4.");
    
    fclose(in);    
    return seqReturned;
}

// The last three parameters may be NULL.
// If NULL then the load is balanced based on the number of sequences
// If not NULL the they contain the size of each work (based on the number of bases).
void processSingleFastaWholeDatabase(Sequence * query, int * first, int * last, int numWorkers){
	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);
    // Create threads to process sequences.
    uint8_t numThreads = (databaseNumSequences < 5 * NUM_THREAD_FOR_PROCESSING)? 1 : NUM_THREAD_FOR_PROCESSING;
    if (first == NULL || numWorkers != numThreads) printf("Using non balanced load for each thread.\n");
    else printf("Using balanced workload for each thread.\n");
    double ratio = 1.0/(double)numThreads;
    pthread_t threads[numThreads];
    paramToSingleQueryProcessThread params[numThreads];
    pthread_mutex_init(& Context.mutex_next_db_seq, NULL);
    Context.next_db_seq_number = 0;
    fprintf(stdout, "\nResults for: %s\n", query->name);
    for(int i=0; i < numThreads; i++){
        params[i].id = i;
        params[i].query = query;
        if (first == NULL || numWorkers != numThreads) {
			params[i].first = databaseNumSequences*ratio*i;
			params[i].last = databaseNumSequences*ratio*(i+1)-1;
        } else {
        	params[i].first = first[i];
        	params[i].last = last[i];
        	// fprintf(stdout, "Worker %d: \t%d \t%d\n", i, first[i], last[i]);
        }
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
    pthread_mutex_destroy(& Context.mutex_next_db_seq);
	gettimeofday(&tiempo_final, NULL);
	double tiempo_total = (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0);
	fprintf(stdout, "Tiempo transcurrido: %9.1fms. \n", tiempo_total);
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

	params->ret.numSequencesProcessed = 0;
	params->ret.numLettersProcessed = 0;
	params->ret.numHits = 0;
	params->ret.numFarrar = 0;
	params->ret.total_time = 0.0;
//    printf("Thread Id %d: start %d end %d\n", params->id, params->first, params->last);
    // For each database sequence in our range
	databaseIdx=params->first - 1;
	const int step = 10;
	while(1) {
//		if (params->query->dataLength > 1000){
			pthread_mutex_lock(& Context.mutex_next_db_seq);
				startIdx = Context.next_db_seq_number;
				Context.next_db_seq_number += step;
			pthread_mutex_unlock(& Context.mutex_next_db_seq);
			if (startIdx >= databaseNumSequences) break;
//		} else {
//			databaseIdx++;
//			if (databaseIdx > params->last) break;
//		}
	for(databaseIdx=startIdx; databaseIdx < startIdx+step; databaseIdx++){
		if (databaseIdx >= databaseNumSequences) break;
//    for(int databaseIdx=params->first; databaseIdx <= params->last; databaseIdx++){
    	params->ret.numSequencesProcessed ++;
    	params->ret.numLettersProcessed += databaseAlignedDemultiplexed[databaseIdx].realDataLength;
//            printf("%s\n", databaseAlignedDemultiplexed[databaseIdx].data);

    	uint16_t nearbyShifter = 0;
        // For each block of 4 consecutive letters in the query
        for(int queryIdx=0; queryIdx < params->query->dataLength - 3; queryIdx++){
//        	printf("Loading queryy %c%c%c%c\n", params->query->data [queryIdx],
//        			params->query->data [queryIdx+1],
//					params->query->data [queryIdx+2],
//					params->query->data [queryIdx+3]);
        	uint32_t block __attribute__((aligned(4))) = * (uint32_t *)&(params->query->data [queryIdx]);
            __m512i queryVector = _mm512_extload_epi32 (&block, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 1);
//            printf("--- "); view512iAsChar(queryVector);
            // We check every block of VECTOR_SIZE bytes
            for(int targetIdx=0; targetIdx < databaseAlignedDemultiplexed[databaseIdx].dataLength ; targetIdx+=VECTOR_SIZE){
//            	printf("TargetIdx %d\n", targetIdx);
            	__m512i targetVector = _mm512_load_epi32 ((__m512i const*) (databaseAlignedDemultiplexed[databaseIdx].data + targetIdx));
//                view512iAsChar(targetVector);
            	__mmask16 res = _mm512_cmpeq_epi32_mask  (queryVector, targetVector);
//                printf("Res %d\n", res);
                if (res != 0) {
//                	printf("%s contains %c%c%c%c\n", databaseAlignedDemultiplexed[databaseIdx].name,
//                        params->query->data[queryIdx],
//                        params->query->data[queryIdx+1],
//                        params->query->data[queryIdx+2],
//                        params->query->data[queryIdx+3]
//                        );
                    if ((nearbyShifter != 0) && (__builtin_popcount (nearbyShifter) >= 3)) {
                    	params->ret.numFarrar ++;
                    	uint32_t score = smith_waterman_farrar(&o, databaseAlignedDemultiplexed[databaseIdx].realData, (int16_t)(databaseAlignedDemultiplexed[databaseIdx].realDataLength));
                    	if (score >= Context.threshold) {
                    		params->ret.numHits ++;
                    		// Show hit with sequence name
                    		fprintf(stdout, "Hit. Score: %d. Sequence: %s\n", score, databaseAlignedDemultiplexed[databaseIdx].name);
                    	}
                        targetIdx = (queryIdx = INT_MAX - 1);
                    }
                    nearbyShifter = (nearbyShifter | 0x01);
                }
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
