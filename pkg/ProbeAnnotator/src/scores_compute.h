#include <R.h>
#ifdef HAVE_PCRE_H
#include <pcre.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <limits.h>
#include "uthash.h"


void scores_compute(probe_t **probes, gene_t **genes, score_t **scores, double *w) {
	int i;
	probe_t * PROBE;
	score_t * SCORE;
	gene_t  * GENE;

	/* iterate each probe */
	for(PROBE=*probes; PROBE != NULL; PROBE=PROBE->hh.next) {
		/* if has mapped genes -> compute score */
		/* else -> no score */
		if(PROBE->nb_genes > 0) { 
			/* initialise best_score and best_index (set to -1 for 'PROBE->best_score < temp_score')*/ 
			PROBE->best_score = -1.0;
			PROBE->best_index = -1; 
			/* compute score for each mapped gene */
			for(i = 0; i < PROBE->nb_genes; ++i) {
				/* get score_t item */
				SCORE = PROBE->scores[i];
		
				/* compute score */
				double temp_score = 0.0;
				for(int j = 0; j < PROBE->nb_probes; ++j) {
					const double temp_d = (SCORE->scores[j] == SCOREMAX) ? 0.0 : w[ (SCORE->scores[j]) ];
					temp_score +=  temp_d;
				}
				
				temp_score = temp_score / (double)(PROBE->nb_probes);
				/* if score is higher than actual score -> set has best score */
				if(PROBE->best_score < temp_score) {
					PROBE->best_score = temp_score;
					PROBE->best_score_s[0] = SCORE;
				}
				
			}
			/* find gene associated with best score */
			HASH_FIND(hh, *genes, &((PROBE->best_score_s[0])->id.gene_id), sizeof(gene_key), GENE); /* struct already in the hash? */
			if(GENE == NULL) /* SHOULD NEVER HAPPEN -> all genes are supposed to be in the hashtable */
			{
			  Rprintf("ERROR:\tCannot find gene %s in the hashtable\n", ((PROBE->best_score_s[0])->id.gene_id).gene_name);
			  return;
			  //exit(0);
			}
			
			if(GENE->best < PROBE->best_score) {
				GENE->best = PROBE->best_score;
				GENE->best_count = 1;
			} 
			else if(GENE->best == PROBE->best_score) {
				GENE->best_count++;
			}
		}
		
	} 
  
}
