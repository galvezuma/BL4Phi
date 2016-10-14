/*
 * ManageExtended.c
 *
 *  Created on: 3 de oct. de 2016
 *      Author: galvez
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>
#include <inttypes.h>

#include "Types.h"
#include "GeneralFunctions.h"
#include "ManageExtended.h"

extern SequenceDemultiplexed * databaseAlignedDemultiplexed;
extern uint32_t databaseNumSequences;

char * bulkNames = NULL;
char * bulkRealData = NULL;
extern char * bulkData = NULL;

/*
 * The idea of this function is to reduce the execution time of the load database.
 * It seems that Xeon-Phi takes a lot of time to execute malloc so here we use a file that
 * contains a structure very similar to that of memory so we do not use time to malloc but to
 * make the pointers to point to the correct position.
 */
SequenceDemultiplexed * loadDatabaseExtended(char * filename){
	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);

    char filenameExtended[MAX_FILENAME_LENGTH];
    sprintf(filenameExtended, "%s.extended", filename);
    FILE * in = fopen(filenameExtended, "rb");
    if (in == NULL) errorAndExit(filenameExtended, "Cannot read from file.");

    FileHeader header;
    fread(&header, sizeof(header), 1, in);
    databaseNumSequences = header.numSequences;
    //printf("%d, %d, %d, %ld\n", header.numSequences, header.fastaHeaderNumBytes, header.fastaRealDataNumBytes, header.fastaDataNumBytes);
    databaseAlignedDemultiplexed = (SequenceDemultiplexed *)malloc(sizeof(SequenceDemultiplexed)*databaseNumSequences);
    //
    bulkNames = (char *) malloc(header.fastaHeaderNumBytes);
    fread(bulkNames, header.fastaHeaderNumBytes, 1, in);
    //
    bulkRealData = (char *) malloc(header.fastaRealDataNumBytes);
    fread(bulkRealData, header.fastaRealDataNumBytes, 1, in);
    //
    bulkData = (char *) _mm_malloc(header.fastaDataNumBytes, 64);
    //fread(bulkData, header.fastaDataNumBytes, 1, in);
    ////
    int offset = 0;
    uint32_t data;
    for(int i=0;i<databaseNumSequences;i++){
    	data = *(uint32_t *)(bulkNames + offset);
    	databaseAlignedDemultiplexed[i].name = bulkNames + offset + sizeof(data);
    	offset += data + sizeof(data);
	}
    //
    offset = 0;
    for(int i=0;i<databaseNumSequences;i++){
        data = *(uint32_t *)(bulkRealData + offset);
    	databaseAlignedDemultiplexed[i].realDataLength = data;
    	databaseAlignedDemultiplexed[i].dataLength = *(uint32_t *)(bulkRealData + offset + sizeof(data));
    	databaseAlignedDemultiplexed[i].realData = bulkRealData + offset + sizeof(data) + sizeof(databaseAlignedDemultiplexed[i].dataLength);
    	offset += data + sizeof(data) + sizeof(databaseAlignedDemultiplexed[i].dataLength);
	}
    //
    offset = 0;
    for(int i=0;i<databaseNumSequences;i++){
    	databaseAlignedDemultiplexed[i].data = bulkData + offset;
    	offset += databaseAlignedDemultiplexed[i].dataLength;
    	if (databaseAlignedDemultiplexed[i].dataLength == 0) continue;
    	// We copy the sequence as explained in the doumentation: demultiplexing
		// Create first copy
				int posSeq = 0;
				uint32_t lengthToCopy = trunc4(databaseAlignedDemultiplexed[i].realDataLength);
				memcpy(&(databaseAlignedDemultiplexed[i].data[posSeq]), databaseAlignedDemultiplexed[i].realData, lengthToCopy);
				// Create second copy
				posSeq += lengthToCopy;
				lengthToCopy = trunc4(databaseAlignedDemultiplexed[i].realDataLength - 1);
				memcpy(&(databaseAlignedDemultiplexed[i].data[posSeq]), databaseAlignedDemultiplexed[i].realData+1, lengthToCopy);
				// Create third copy
				posSeq += lengthToCopy;
				lengthToCopy = trunc4(databaseAlignedDemultiplexed[i].realDataLength - 2);
				memcpy(&(databaseAlignedDemultiplexed[i].data[posSeq]), databaseAlignedDemultiplexed[i].realData+2, lengthToCopy);
				// Create fourth copy
				posSeq += lengthToCopy;
				lengthToCopy = trunc4(databaseAlignedDemultiplexed[i].realDataLength - 3);
				memcpy(&(databaseAlignedDemultiplexed[i].data[posSeq]), databaseAlignedDemultiplexed[i].realData+3, lengthToCopy);
				// Fill the final gaps
				posSeq += lengthToCopy;
				while (posSeq%64 != 0) databaseAlignedDemultiplexed[i].data[posSeq++] = 0;
				if (databaseAlignedDemultiplexed[i].dataLength != posSeq)  {
					printf("Seq of length %d, name %s\n%s\n%d\n%d\n\n",
							databaseAlignedDemultiplexed[i].realDataLength,
							databaseAlignedDemultiplexed[i].name,
							databaseAlignedDemultiplexed[i].realData,
							posSeq,
							databaseAlignedDemultiplexed[i].dataLength);
					errorAndExit("", "Internal error: inconsistency in alignments.");
				}
        // Finish of demultiplexing
	}

    printf("%s\n%d\n%d\n%p\n", databaseAlignedDemultiplexed[125678].name, databaseAlignedDemultiplexed[125678].realDataLength, databaseAlignedDemultiplexed[125678].dataLength, databaseAlignedDemultiplexed[125678].data);

    fclose(in);
    gettimeofday(&tiempo_final, NULL);
    printf("Time to load database: %5.3f\n", (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0));
    return NULL;
}


void freeDatabaseExtended(){
	_mm_free(bulkData); bulkData = NULL;
	free(bulkRealData); bulkRealData = NULL;
	free(bulkNames); bulkNames = NULL;
	free(databaseAlignedDemultiplexed); databaseAlignedDemultiplexed = NULL;
}


void checkDatabaseExtended(){
    for(int i=0;i<databaseNumSequences;i++){
    	uint32_t data = excess64(databaseAlignedDemultiplexed[i].realDataLength * 4 -12);
    	if (data != databaseAlignedDemultiplexed[i].dataLength) {
    		printf("Bad data lengths (%s): %d - %d\n",databaseAlignedDemultiplexed[i].name, databaseAlignedDemultiplexed[i].realDataLength, databaseAlignedDemultiplexed[i].dataLength);
    		exit(1);
    	}
    	if ((uint64_t)databaseAlignedDemultiplexed[i].data % 64 != 0) {
    		printf("Bad data alignment (%d)(%s).\n",i, databaseAlignedDemultiplexed[i].name);
    		exit(1);
    	}
    	if (databaseAlignedDemultiplexed[i].name[0] != '>') {
    		printf("Bad header (%d)(%s).\n",i, databaseAlignedDemultiplexed[i].name);
    		exit(1);
    	}
	}
}
