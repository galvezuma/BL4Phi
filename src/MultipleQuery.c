/*
 * MultipleQuery.c
 *
 *  Created on: 13 de sept. de 2020
 *      Author: galvez
 */

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

Sequence * loadNextQueryFromFasta(char * filename) {
	static FILE * in = NULL;
	static int end_of_quieries = 0;
    static char line[MAX_LINE_LENGTH];

    Sequence * seqReturned;
    char sequence[MAX_SEQUENCE_LENGTH];
    char * okReading;

    // Open file for the very first time
    if (in == NULL) {
    	if (end_of_quieries)
    		return NULL;
    	in = fopen(filename, "rb");
    	if (in == NULL) errorAndExit(filename, "Cannot open file.");
        // Reading the name of the sequence
        okReading = fgets(line, MAX_LINE_LENGTH, in);
        if ((okReading == NULL) || (line[0] != '>')) errorAndExit(line, "Not a fasta header.");
    } else { // The file was already open
    	;
    }

    seqReturned = (Sequence *)malloc(sizeof(Sequence));
    if (seqReturned == NULL) errorAndExit("MultipleQuery", "Fatal error: Out of memory.");
    line[strcspn(line, "\r\n")] = 0;
    uint32_t lengthText = strlen(line);
    seqReturned->name = (char *) malloc(lengthText + 1);
    if (seqReturned->name == NULL) errorAndExit("MultipleQuery.name", "Fatal error: Out of memory.");
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

    // We close definitely the file when there is no more query sequence to read
    if (line[0] != '>') {
    	fclose(in);
    	in = NULL;
    	end_of_quieries = 1;
    }

    return seqReturned;
}
