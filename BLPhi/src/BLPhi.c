/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.c
 * Author: galvezs
 *
 * Created on 01 June 2016, 14:03
 */

/*
 * TO DO: Optimizations to do:
 *
 * 1.- OTHER SOLUTION GIVEN. Mixed use of threads for several queries. This should avoid the effect of very busy threads
 * at the end of the processing.
 *
 * 2.- Usage of different types of nearby shifters.
 *
 * 3.- KNL. Usage of amino acids letters in 4 bits instead of 8 bits. This should divide by two the execution time.
 * Intrinsic of epi16 must be used.
 *
 * 4.- Save pre-built structures to load it quickly at processing time. Actually, this is not a pre indexation
 * but a memory dump.
 *
 * 5.- Partition of a big database into parts.
 *
 * 6.- A good wish. Utilization of several Xeon-Phi cards. I do not like this because it forces to coordinate the cards from
 * outside, whereas we have programmed the cards in native mode.
 *
 * 7.- Future one. Provide the actual alignment including data: start, end, etc.
 *
 * 8.- Usage like blastx.
 *
 * Please, do not forget to use the latest release of Swiss Prot.
 *
 * 9.- Transformation to use nucleotides instead of amino acids. I think that the use of nibbles for this may be
 * the best approach: we require a match of 8 letters and we may require three or four matches in a
 * particular nearby region.
 *
 */

/*
 * Release 2016_08 of 07-Sep-16 of UniProtKB/Swiss-Prot contains 551987 sequence entries,
 * comprising 197275398 amino acids abstracted from 246580 references.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <malloc.h>
#include <immintrin.h>

#include "Types.h"
#include "GeneralFunctions.h"
#include "SingleQuery.h"
#include "MultipleQuery.h"
#include "ManageStats.h"
#include "ManageExtended.h"
#include "Farrar.h"


SequenceDemultiplexed * databaseAlignedDemultiplexed = NULL;
uint32_t databaseNumSequences;

// Read only in the threads
struct {
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

void initContext(){
	Context.open_gap_cost = 10;
	Context.extend_gap_cost = 1;
	Context.threshold = 120;
	Context.sec_ref = NULL;
	Context.matrix = NULL;
	Context.num_letters = 0;
	Context.Farrar.letras_sec_ref = NULL;
	Context.Farrar.long_profile = 0;
	Context.Farrar.profile = NULL;
}

/*
 *
 */
int main2(int argc, char** argv) {
//    createStats("/mic0fs/blphi/uniprot_sprot.fasta");
//    checkcreateStats();

	loadDatabaseExtended("/mic0fs/blphi/uniprot_sprot.fasta");
	checkDatabaseExtended();

    return (EXIT_SUCCESS);

	int numWorkers = NUM_THREAD_FOR_PROCESSING;
	int first[numWorkers], last[numWorkers];
	loadFastaAndBalanceLoad("/mic0fs/blphi/uniprot_sprot.fasta", first, last, numWorkers);
    printf("Ending...\n");
//    checkLoadFastaAndBalanceLoad(first, last, numWorkers);

    Sequence * query = loadSingleFasta("/mic0fs/blphi/SingleQuery.fasta");
    checkloadSingleFasta(query);
    processSingleFastaWholeDatabase(query, first, last, numWorkers);
//    freeSingleFasta(query);
//
//    printf("Freeing memory");
//    for(int i=0;i<databaseNumSequences; i++) {
//        free(databaseAlignedDemultiplexed[i].name);
//        _mm_free(databaseAlignedDemultiplexed[i].data);
//    }
//    free(databaseAlignedDemultiplexed);
    return (EXIT_SUCCESS);
}

int main(int argc, char** argv) {
	  char matrix_filename[MAX_LINE_LENGTH] = "BLOSUM60";
	  initContext();
	  // PARAMETERS MANAGEMENT
	  int pos = 1;
	  while(pos < argc - 2){
	    if (!strcmp(argv[pos], "-t")) { // Read threshold
	      Context.threshold = atoi(argv[pos+1]);
	      pos += 2;
	    } else if (!strcmp(argv[pos], "-m")) { // Read matrix
	      strcpy(matrix_filename, argv[pos+1]);
	      pos += 2;
	    } else if (!strcmp(argv[pos], "-g")) { // Read costs
	      Context.open_gap_cost = atoi(argv[pos+1]);
	      Context.extend_gap_cost = atoi(argv[pos+2]);
	      fprintf(stdout, "Open gap set to: %d\nExtend gap set to: %d\n", Context.open_gap_cost, Context.extend_gap_cost);
	      pos += 3;
	    } else {
	        fprintf(stderr, "Invalid parameter: %s\n", argv[pos]);
	        return 1;
	    }
	  }
	  if (argc < pos + 2) {
	    fprintf(stderr, "Usage: BLPhi [-t threshold] [-m matrix] [-g open_gap_cost extend_gap_cost] query.fasta database.fasta\n");
	    return 1;
	  }
	  if (((uint8_t)Context.open_gap_cost) > 0xFF) { fprintf(stderr, "Open gap cost (%d) too big. The size is a byte.\n", Context.open_gap_cost); return 1; }
	  if (((uint8_t)Context.extend_gap_cost) > 0xFF) { fprintf(stderr, "Extend gap cost (%d) too big. The size is a byte.\n", Context.extend_gap_cost); return 1; }
	  if (Context.threshold > 0xFFFF) { fprintf(stderr, "Threshold (%d) too big. Maximum value is 0xFFFF.\n", Context.threshold); return 1; }

	  Context.matrix = load_matrix(matrix_filename); // This assignment is redundant
	  if (Context.matrix == NULL) { return 1; }
	  int numWorkers = NUM_THREAD_FOR_PROCESSING;
	  int first[numWorkers], last[numWorkers];
	  loadDatabaseExtended(argv[pos+1]);
	  //loadFastaAndBalanceLoad(argv[pos+1], first, last, numWorkers);


	  while ((Context.sec_ref = loadNextQueryFromFasta(argv[pos])) != NULL) {
		  if (Context.sec_ref == NULL) { return 1; }
		  // OPTIMIZATIONS AND TRANSFORMATIONS FOR FARRAR
		  transformLettersIntoPositions();
		  calculateFarrarProfile();

		  /*
		   * Main loop that iterates over the Query Fasta file for each query sequence
		   */
		  //processSingleFastaWholeDatabase(Context.sec_ref, first, last, numWorkers);
		  processSingleFastaWholeDatabase(Context.sec_ref, NULL, NULL, numWorkers);
		  freeSingleFasta(Context.sec_ref);
		  prepareForNextQuery();
	  }


	  //TAREAS DE FINALIZACIÓN
	  printf("Freeing memory\n");
//	  for(int i=0;i<databaseNumSequences; i++) {
//	      free(databaseAlignedDemultiplexed[i].name);
//	      free(databaseAlignedDemultiplexed[i].realData);
//	      _mm_free(databaseAlignedDemultiplexed[i].data);
//	  }
//	  free(databaseAlignedDemultiplexed);
	  freeDatabaseExtended();
	  fflush(stdout);
	  fflush(stderr);
	  free(Context.matrix);

	  return (EXIT_SUCCESS);
}

