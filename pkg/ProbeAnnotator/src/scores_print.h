#include <R.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "uthash.h"


void scores_print(probe_t **probes, gene_t **genes, score_t **scores, char **file_out) {
FILE *fp2;

gene_t * GENE;
probe_t * PROBE;

fp2=fopen(*file_out, "w");
fprintf(fp2, "probeset_id\tbest_match\tfinal_score\tbest_match_best_score\twell_mapped\tgene_best_map\tgene_unique_best_map\tgene_unique_map\n"); 
for(PROBE=*probes; PROBE != NULL; PROBE=PROBE->hh.next) {
	if(PROBE->nb_genes > 0) {
	HASH_FIND(hh, *genes, &((PROBE->best_score_s[0])->id.gene_id), sizeof(gene_key), GENE);
	//g = s->genes[s->best_index];
	fprintf(fp2, "%s\t%s\t%f\t%f\t%s\t%s\t%s\t%s\n", 
	  PROBE->name,((PROBE->best_score_s[0])->id.gene_id).gene_name, PROBE->best_score, GENE->best, PROBE->best_score == 1 ? "TRUE" : "FALSE", PROBE->nb_genes == 1 ? "TRUE" : "FALSE", GENE->count == 1 ? "TRUE": "FALSE", GENE->best_count == 1 ? "TRUE" : "FALSE"); 
	} else {
	fprintf(fp2, "%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", 
	  PROBE->name); 
	}	
}
fclose(fp2);  
}
