/*
 * ManageStats.c
 *
 *  Created on: 9 de sept. de 2016
 *      Author: galvez
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <malloc.h>
#include <sys/time.h>

#include "Types.h"
#include "GeneralFunctions.h"
#include "SingleQuery.h"
#include "ManageStats.h"

//#define NON_EXHAUSTIVE

/*
 * This function is called to load a Fasta database and, at the same time, to split the database into numWorkers
 * blocks each of them with a similar number of letters to process.
 * The result of the splitting is stored into the arrays first and last.
 * Note that numWorkers is the number of workers that will work later in the search algorithm.
 */
void loadFastaAndBalanceLoad(char *filename, int * first, int * last, int numWorkers) {
	loadFastaFromStats(filename);

    char filenameStats[MAX_FILENAME_LENGTH];
    sprintf(filenameStats, "%s.stats", filename);
    FILE * inx = fopen(filenameStats, "rb");
    if (inx == NULL) errorAndExit(filenameStats, "Cannot open file.");
    uint32_t data, lengthSequence;
    fread(&data, sizeof(data), 1, inx); // There are data sequences to read

    uint64_t numLetters = 0;
    fread(&numLetters, sizeof(numLetters), 1, inx); // There are numLetters letters to work with
    printf("Letters in database: %ld\n", numLetters);

    double ratio = (double)numLetters/(double)numWorkers;
    int currentWorker = 0, currentLength = 0;
    first[0] = 0;
    for(int i=0; i<databaseNumSequences; i++){
    	fread(&data, sizeof(data), 1, inx); // Skip posHeader
    	fread(&data, sizeof(data), 1, inx); // Skip lengthHeader
    	fread(&data, sizeof(data), 1, inx); // Skip posSequence
    	fread(&lengthSequence, sizeof(lengthSequence), 1, inx); // Read lengthSequence
    	currentLength += lengthSequence;
    	if (currentLength > ratio * (currentWorker + 1)) {
    		last[currentWorker] = i;
    		currentWorker ++;
    		if (currentWorker >= numWorkers) break; // To avoid problems at the end of the loop.
    		first[currentWorker] = i+1;
    	}
    }
    last[numWorkers-1] = databaseNumSequences - 1;

    fclose(inx);
}

void checkLoadFastaAndBalanceLoad(int * first, int * last, int numWorkers) {
	checkloadFastaFromStats();
	for (int i=0; i<numWorkers; i++){
		fprintf(stdout, "Worker %d: \t%d \t%d\n", i, first[i], last[i]);
	}
}

void loadFastaFromStats(char *filename) {
	struct timeval tiempo_inicio, tiempo_final;
	gettimeofday(&tiempo_inicio, NULL);

    char filenameStats[MAX_FILENAME_LENGTH];
    sprintf(filenameStats, "%s.stats", filename);
    FILE * inx = fopen(filenameStats, "rb");
    if (inx == NULL) errorAndExit(filenameStats, "Cannot open file.");
    fread(&databaseNumSequences, sizeof(databaseNumSequences), 1, inx); // There are data sequences to read
    printf("Sequences in database: %d\n", databaseNumSequences);
    fclose(inx);

    // We create the structure that will be fulfilled by the threads
    databaseAlignedDemultiplexed = (SequenceDemultiplexed *)malloc(sizeof(SequenceDemultiplexed)*databaseNumSequences);

    // Create 4 threads to load sequences into memory if needed.
    uint8_t numThreads = (databaseNumSequences < 1000)? 1 : NUM_THREAD_FOR_LOADING;
    pthread_t threads[numThreads];
    paramToLoadThread params[numThreads];
    for(int i=0; i < numThreads; i++){
        params[i].filename = filename;
        params[i].id = i;
        params[i].arraySeqAlgDemul = databaseAlignedDemultiplexed;
        double ratio = 1.0/(double)numThreads;
        params[i].first = databaseNumSequences*ratio*i;
        params[i].last = databaseNumSequences*ratio*(i+1)-1;
        pthread_create(&threads[i], NULL, loadBunchFromStats, (void *)&params[i]);
    }
    // We need to wait the threads to finish because we have passed local variables as parameters
    for(int i=0; i < numThreads; i++){
        pthread_join(threads[i], NULL);
        printf("Joining...\n");
    }
    gettimeofday(&tiempo_final, NULL);
    printf("Time to load database: %5.3f\n", (double) (tiempo_final.tv_sec - tiempo_inicio.tv_sec) * 1000 + ((double) (tiempo_final.tv_usec - tiempo_inicio.tv_usec) / 1000.0));
}

void * loadBunchFromStats(void * vparams){
    paramToLoadThread * params = (paramToLoadThread *)vparams;
    char filenameStats[MAX_FILENAME_LENGTH];
    sprintf(filenameStats, "%s.stats", params->filename);
    FILE * inx = fopen(filenameStats, "rb");
    if (inx == NULL) errorAndExit(filenameStats, "Cannot open file.");
    uint32_t data;
    uint64_t big_data;
    fread(&data, sizeof(data), 1, inx); // There are data sequences to read, from first to last
    fread(&big_data, sizeof(big_data), 1, inx); // There are big_data letters to read, from first to last
    printf("Thread Id %d: start %d end %d\n", params->id, params->first, params->last);

    // Let's process the sequences
    fseek(inx, sizeof(uint32_t)*(3+4*params->first), SEEK_SET); // We put at the first index
    fread(&data, sizeof(data), 1, inx);
    fclose(inx);
    // We close the index because we need it only to positionate at the correct position in the fasta file.
    FILE * in = fopen(params->filename, "r");
    if (in == NULL) errorAndExit(params->filename, "File not found.");
    fseek(in, data, SEEK_SET);
    printf("Positioning at: %d\n", data);
    char line[MAX_LINE_LENGTH];
    char sequence[MAX_SEQUENCE_LENGTH];
    char * okReading;
    okReading = fgets(line, MAX_LINE_LENGTH, in);
    for(int i=params->first; i <= params->last; i++) {
        //// Read sequence and put it in the array
                if ((okReading == NULL) || (line[0] != '>')) errorAndExit(line, "Not a fasta header.");
                line[strcspn(line, "\r\n")] = 0;
                uint32_t lengthText = strlen(line);
                params->arraySeqAlgDemul[i].name = (char *) malloc(lengthText + 1);
                if (params->arraySeqAlgDemul[i].name == NULL) errorAndExit("", "Fatal error: Out of memory.");
                strcpy(params->arraySeqAlgDemul[i].name, line);
                uint32_t posSeq = 0;
                okReading = fgets(line, MAX_LINE_LENGTH, in);
                while((line[0] != '>') && okReading != NULL){
                    line[strcspn(line, "\r\n")] = 0;
                    lengthText = strlen(line);
                    if (posSeq+lengthText > MAX_SEQUENCE_LENGTH) errorAndExit(params->arraySeqAlgDemul[i].name, "Fatal error. Sequence too long.");
                    strcpy(sequence+posSeq, line);
                    posSeq += lengthText;
                    okReading = fgets(line, MAX_LINE_LENGTH, in);
                }
                lengthText = posSeq;
                if (lengthText < 4) {
                    params->arraySeqAlgDemul[i].dataLength = 0;
                    params->arraySeqAlgDemul[i].data = NULL;
                    continue;
                }
                params->arraySeqAlgDemul[i].realDataLength = lengthText;
                params->arraySeqAlgDemul[i].realData = (char *) malloc(sizeof(char)*lengthText);
                memcpy(params->arraySeqAlgDemul[i].realData, sequence, lengthText);
#ifdef NON_EXHAUSTIVE
                // We skip the 2nd and 4th copies shifted
                int lengthFirst = trunc4(lengthText);
                int lengthThird = trunc4(lengthText - 2);
                params->arraySeqAlgDemul[i].dataLength = excess64(lengthFirst + lengthThird);
                params->arraySeqAlgDemul[i].data = (char *) _mm_malloc(params->arraySeqAlgDemul[i].dataLength, 64);
                if (params->arraySeqAlgDemul[i].data == NULL) errorAndExit("", "Fatal error: Out of memory.");
                // We copy the sequence as explained in the doumentation: demultiplexing
                // but copies 1st and 3rd only.
                	// Create first copy
                	posSeq = 0;
					uint32_t lengthToCopy = trunc4(lengthText);
					memcpy(&(params->arraySeqAlgDemul[i].data[posSeq]), sequence, lengthToCopy);
                    // Create third copy
                    posSeq += lengthToCopy;
                    lengthToCopy = trunc4(lengthText - 2);
                    memcpy(&(params->arraySeqAlgDemul[i].data[posSeq]), sequence+2, lengthToCopy);
                    // Fill the final gaps
                    posSeq += lengthToCopy;
                    while (posSeq%64 != 0) params->arraySeqAlgDemul[i].data[posSeq++] = 0;
                    if (params->arraySeqAlgDemul[i].dataLength != posSeq)  {
                        printf("Seq %d (%d), name %s\n%s\n%d\n%d\n\n", i, lengthText, params->arraySeqAlgDemul[i].name, sequence, posSeq, params->arraySeqAlgDemul[i].dataLength);
                        errorAndExit("", "Internal error: inconsistency in alignments.");
                    }
                // Finish of demultiplexing
#else
                params->arraySeqAlgDemul[i].dataLength = excess64(4*lengthText-12);
                params->arraySeqAlgDemul[i].data = (char *) _mm_malloc(params->arraySeqAlgDemul[i].dataLength, 64);
                if (params->arraySeqAlgDemul[i].data == NULL) errorAndExit("", "Fatal error: Out of memory.");
                // We copy the sequence as explained in the doumentation: demultiplexing
                    // Create first copy
                    posSeq = 0;
                    uint32_t lengthToCopy = trunc4(lengthText);
                    memcpy(&(params->arraySeqAlgDemul[i].data[posSeq]), sequence, lengthToCopy);
                    // Create second copy
                    posSeq += lengthToCopy;
                    lengthToCopy = trunc4(lengthText - 1);
                    memcpy(&(params->arraySeqAlgDemul[i].data[posSeq]), sequence+1, lengthToCopy);
                    // Create third copy
                    posSeq += lengthToCopy;
                    lengthToCopy = trunc4(lengthText - 2);
                    memcpy(&(params->arraySeqAlgDemul[i].data[posSeq]), sequence+2, lengthToCopy);
                    // Create fourth copy
                    posSeq += lengthToCopy;
                    lengthToCopy = trunc4(lengthText - 3);
                    memcpy(&(params->arraySeqAlgDemul[i].data[posSeq]), sequence+3, lengthToCopy);
                    // Fill the final gaps
                    posSeq += lengthToCopy;
                    while (posSeq%64 != 0) params->arraySeqAlgDemul[i].data[posSeq++] = 0;
                    if (params->arraySeqAlgDemul[i].dataLength != posSeq)  {
                        printf("Seq %d (%d), name %s\n%s\n%d\n%d\n\n", i, lengthText, params->arraySeqAlgDemul[i].name, sequence, posSeq, params->arraySeqAlgDemul[i].dataLength);
                        errorAndExit("", "Internal error: inconsistency in alignments.");
                    }
                // Finish of demultiplexing
#endif
        ////
    }
    fclose(in);
    //pthread_exit(NULL);
    return NULL;
}

void checkloadFastaFromStats(){
    for(int i=0; i < 10; i++) {
        printf("%s\n", databaseAlignedDemultiplexed[i].name);
    }
    for(int i=0; i < 10; i++) {
        printf("%s\n", databaseAlignedDemultiplexed[databaseNumSequences-1-i].name);
    }
}

void createStats(char *filename) {
    FILE * in = fopen(filename, "r");
    if (in == NULL) errorAndExit(filename, "File not found.");
    char filenameStats[MAX_FILENAME_LENGTH];
    // TODO Break if the file size is greater than 32 bits (4 GiB))
    sprintf(filenameStats, "%s.stats", filename);
    FILE * out = fopen(filenameStats, "wb");
    if (out == NULL) errorAndExit(filenameStats, "Cannot write on file.");

    char line[MAX_LINE_LENGTH];
    uint32_t numSequences = 0;
    uint64_t numBases = 0;
    uint32_t outOffset = 0;
    char * okReading;
    okReading = fgets(line, MAX_LINE_LENGTH, in);
    // Write something to leave space for the number of sequences.
    fwrite(&numSequences, sizeof(numSequences), 1, out);
    // Write something to leave space for the number of letters.
    fwrite(&numBases, sizeof(numBases), 1, out);
    do { // While there are lines in the fasta file
        if ((okReading == NULL) || (line[0] != '>')) errorAndExit(line, "Not a fasta header.");
        numSequences++;
        uint32_t lengthText = strlen(line);
        fwrite(&outOffset, sizeof(outOffset), 1, out); // Write posHeader
        fwrite(&lengthText, sizeof(lengthText), 1, out); // Write lengthHeader
        outOffset = (uint32_t) ftell(in);
        fwrite(&outOffset, sizeof(outOffset), 1, out); // Write posSequence
        lengthText = 0;
        okReading = fgets(line, MAX_LINE_LENGTH, in);
        while((line[0] != '>') && okReading != NULL){
            line[strcspn(line, "\r\n")] = 0;
            lengthText += strlen(line);
            outOffset = (uint32_t) ftell(in);
            okReading = fgets(line, MAX_LINE_LENGTH, in);
        }
        fwrite(&lengthText, sizeof(lengthText), 1, out); // Write lengthSequence
        numBases += lengthText;
    } while(okReading != NULL);
    fclose(out);
    fclose(in);
    out = fopen(filenameStats, "r+b");
    fwrite(&numSequences, sizeof(numSequences), 1, out); // Write num of sequences at the beginning of the file
    fwrite(&numBases, sizeof(numBases), 1, out); // and, then, num of letters
    printf("Num sequences: %d\n", numSequences);
    printf("Num bases: %ld\n", numBases);
    fclose(out);

}

void checkcreateStats() {

    FILE * in = fopen("/mic0fs/blphi/uniprot_sprot.fasta", "r");
    FILE * inx = fopen("/mic0fs/blphi/uniprot_sprot.fasta.stats", "r");
    fseek(inx, 12+10040*16, SEEK_SET);
    uint32_t data;
    fread(&data, sizeof(data), 1, inx);
    printf("%d\n", data);
    fseek(in, data, SEEK_SET);
    char c;
    fread(&c, sizeof(c), 1, in);
    // The character > must appear
    printf("%c\n", c);
    fread(&data, sizeof(data), 1, inx);
    printf("name length %d\n", data);
    fread(&data, sizeof(data), 1, inx);
    fread(&data, sizeof(data), 1, inx);
    printf("sequence length %d\n", data);
    fclose(inx);
    fclose(in);
}
