/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.c
 * Author: Sergio Gálvez Rojas
 *
 * Created on 01 June 2016, 14:03
 */

/*
 * TO DO: Optimizations to do:
 *
 * 1.- Partition of a big database into parts.
 *
 * 2.- Future one. Provide the actual alignment including data: start, end, etc.
 *
 * 3.- Usage like blastx.
 *
 * Please, do not forget to use the latest release of Swiss Prot.
 *
 * 4.- Transformation to use nucleotides instead of amino acids. I think that the use of nibbles for this may be
 * the best approach: we require a match of 8 letters and we may require three or four matches in a
 * particular nearby region.
 *
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
#include "ManageDatabase.h"
#include "Farrar.h"


SequenceDemultiplexed * databaseAlignedDemultiplexed = NULL;
uint32_t databaseNumSequences;

// Context is mainly Read Only in the threads (but -best-)
// This alignment is critical for performance
//#pragma pack(1)
__declspec(align(64)) struct GlobalContext Context;

void initContext(){
	loadCluster(21);
	Context.open_gap_cost = 10;
	Context.extend_gap_cost = 1;
	Context.threshold = 80;
	Context.nearby = 4;
	Context.nonExhaustive = 0;
	Context.bestOnly = 0;
	Context.best_score = 0;
	Context.best_databaseIdx = -1;
	Context.sec_ref = NULL;
	Context.matrix = NULL;
	Context.num_letters = 0;
	Context.Farrar.letras_sec_ref = NULL;
	Context.Farrar.long_profile = 0;
	Context.Farrar.profile = NULL;
}

int main(int argc, char** argv) {
	  char matrix_filename[MAX_LINE_LENGTH] = "BLOSUM62";
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
	    } else if (!strcmp(argv[pos], "-n")) { // Read nearby number of matches
	    	Context.nearby = atoi(argv[pos+1]);
		    pos += 2;
	    } else if (!strcmp(argv[pos], "-c")) { // Load cluster
	    	loadCluster(atoi(argv[pos+1]));
		    pos += 2;
		} else if (!strcmp(argv[pos], "-f")) { // Set non exhaustive
	    	Context.nonExhaustive = 1;
		    pos += 1;
		} else if (!strcmp(argv[pos], "-b")) { // Set display the best hit only
	    	Context.bestOnly = 1;
		    pos += 1;
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
	    fprintf(stderr, "Usage: BLVector [-t threshold] [-c cluster] [-m matrix] [-n nearby_per_16] [-f(ast:non exhaustive)] [-b(est hit only)] [-g open_gap_cost extend_gap_cost] query.fasta database.fasta\n");
	    return 1;
	  }
	  if (((uint8_t)Context.open_gap_cost) > 0xFF) { fprintf(stderr, "Open gap cost (%d) too big. The size is a byte.\n", Context.open_gap_cost); return 1; }
	  if (((uint8_t)Context.extend_gap_cost) > 0xFF) { fprintf(stderr, "Extend gap cost (%d) too big. The size is a byte.\n", Context.extend_gap_cost); return 1; }
	  if (Context.threshold > 0xFFFF) { fprintf(stderr, "Threshold (%d) too big. Maximum value is 0xFFFF.\n", Context.threshold); return 1; }

	  Context.matrix = load_matrix(matrix_filename); // This assignment is redundant
	  if (Context.matrix == NULL) { return 1; }
	  int numWorkers = NUM_THREAD_FOR_PROCESSING;
	  int first[numWorkers], last[numWorkers];
	  loadDatabase(argv[pos+1]);

	  while ((Context.sec_ref = loadNextQueryFromFasta(argv[pos])) != NULL) {
		  if (Context.sec_ref == NULL) { return 1; }
		  // OPTIMIZATIONS AND TRANSFORMATIONS FOR FARRAR
		  transformLettersIntoPositions();
		  calculateFarrarProfile();

		  /*
		   * Main loop that iterates over the Query Fasta file for each query sequence
		   */
		  processSingleFastaWholeDatabase(Context.sec_ref, NULL, NULL, numWorkers);
		  freeSingleFasta(Context.sec_ref);
		  prepareForNextQuery();
	  }


	  // TASKS TO FINISH
	  printf("Freeing memory\n");
	  freeDatabase();
	  fflush(stdout);
	  fflush(stderr);
	  free(Context.matrix);

	  return (EXIT_SUCCESS);
}

