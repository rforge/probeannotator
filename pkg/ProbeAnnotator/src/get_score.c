#include <R.h>
#ifdef HAVE_PCRE_H
#include <pcre.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "probes_get.h"
#include "scores_get.h"
#include "scores_compute.h"
#include "scores_print.h"



void score(char **fasta_file, char **bowtie_map_file, int* nb_map, char **pattern_gene, char **pattern_probe, int *hasProbeSet, double *w, int *maxMissmatch, char **file_out) {
#ifndef HAVE_PCRE_H
Rprintf("PCRE library not found during compilation. Function empty");
return;
#endif
  
#ifdef HAVE_PCRE_H
  /* There are 8 steps in this function  */
  
  
  /* STEP 0 - INITIALISE VARIABLES */
  int max_missmatch = (maxMissmatch[0]);
  int hps = (hasProbeSet[0]);
  unsigned int num_users;
  
  Rprintf("GET PROBES MAPPING\n|Is Probeset: %d\n|Max Missmatch: %d\n|Weights: ", hps, max_missmatch);
  for(int i = 0; i < max_missmatch-1; ++i) Rprintf("%f-", w[i]); 
  Rprintf("%f\n", w[max_missmatch-1]);
  Rprintf("|Pattern Probe: %s\n|Pattern Gene: %s\n", pattern_probe[0], pattern_gene[0]);
  
  
  /* HASHTABLES (see 'structs.h' for details)*/
  /* probes */
  probe_t *probes = NULL;
  /* genes */
  gene_t *genes = NULL;
  /* scores */
  score_t *scores= NULL;
  
  
  /* STEP 1 - GET PROBES */
  Rprintf("- Get Probes\n");
  if(hps == 0) { /* when there are no probeset */
    get_probes(fasta_file, &probes);
  } 
  else if(hps == 1) { /* when there are probeset */
    get_probes_set(fasta_file, pattern_probe, &probes);
  }
  
  /* STEP 2 - COMPUTE SCORES 1 (get)*/
  Rprintf("- Compute Scores\n");
  if(hps == 0) { /* when there are no probeset */
    scores_get(bowtie_map_file, nb_map[0], &probes, &genes, &scores, max_missmatch, pattern_gene);
  } 
  else if(hps == 1) { /* when there are probeset */
    scores_get_set(bowtie_map_file, nb_map[0], &probes, &genes, &scores, max_missmatch, pattern_gene, pattern_probe);
  }
  
  
  
  /* STEP 3 - COMPUTE SCORES 2 (~compute quality string) s*/
  Rprintf("- Compute Quality String\n");
  scores_compute(&probes, &genes, &scores, w);
  
  /* STEP 4 - PRINT OUTPUT FILE */
  Rprintf("- Print Output File\n");
  scores_print(&probes, &genes, &scores, file_out);
  
  /* STEP 5 - print result info in console */
  /* print STEP 1 */
  num_users = HASH_COUNT(probes);
  Rprintf("|Nb probes: %u.\n", num_users);
  /* print STEP 2 */
  num_users = HASH_COUNT(genes);
  Rprintf("|Nb genes: %u.\n", num_users);
  /* print STEP 3 */
  num_users = HASH_COUNT(scores);
  Rprintf("|Nb scores: %u.\n", num_users);
  
  
  /* STEP 6- FREE */
  /* free probes */
  probe_t *current, *temp;
  HASH_ITER(hh, probes, current, temp) {
    HASH_DEL(probes,current);  /* delete; users advances to next */
    free(current);            /* optional- if you want to free  */
  }
  /* free genes */
  gene_t *current_g, *temp_g;
  HASH_ITER(hh, genes, current_g, temp_g) {
    HASH_DEL(genes,current_g);  /* delete; users advances to next */
    free(current_g);            /* optional- if you want to free  */
  }
  /* free scores */
  score_t *current_t, *temp_t;
  HASH_ITER(hh, scores, current_t, temp_t) {
    HASH_DEL(scores,current_t);  /* delete; users advances to next */
    free(current_t);            /* optional- if you want to free  */
  }
#endif
}

