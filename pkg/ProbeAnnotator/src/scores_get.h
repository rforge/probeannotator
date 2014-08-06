#include <R.h>
#include "pcre.h"
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "uthash.h"
#include <limits.h>


void scores_get(char **bowtie_file, int nb_map, probe_t **probes, gene_t **genes, score_t **scores, int max_missmatch, char **pattern_gene) {

  
/* PCRE+variables */
pcre *re_gene;
const char *pcreErrorStr;
int pcreErrorOffset;
int subStrVec[30];
int pcreExecRet;
const char *match_gene;
/* compile PCRE */
re_gene = pcre_compile(pattern_gene[0], 0, &pcreErrorStr, &pcreErrorOffset, NULL);

//probe
probe_t *temp1;
score_t *temp3;
gene_t *temp2;

//connection+variables
FILE *fp1;
char line[LINEMAX];
//for strtol
char *pch;

for(int map = 0; map < nb_map; map++) {
fp1=fopen(bowtie_file[map], "r");

while(fgets (line, LINEMAX, fp1)!=NULL ) {
	size_t linelen = strlen(line);
	
	// remove end of line character
	if (line[linelen - 1] == '\n') line[linelen - 1] = '\0';

	char *str_probe  = strtok(line, "\t");
	char *str_strand = strtok(NULL, "\t");
	char *str_gene   = strtok(NULL, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");
	char *str_missmatch = strtok(NULL, "\t");
	int missmatch = strtol(str_missmatch, &pch, 10);  	
	int number = 1;

	
	/* if missmatch is lower than max_missmatch */ 
	if(missmatch < max_missmatch) {	
		/* gene PCRE */
		pcreExecRet = pcre_exec(re_gene, NULL, str_gene, strlen(str_gene), 0, 0, subStrVec, 30);  
		if(pcreExecRet>0) {
			pcre_get_substring(str_gene, subStrVec, pcreExecRet, 0, &(match_gene));
		} else {
			Rprintf("NO1 %d %s %d\n", number, str_gene, index);
			break;
		}
		  
		
		/* FIND probe */
		HASH_FIND_STR(*probes, str_probe, temp1);  /* name already in the hash? */
		if(temp1 == NULL) /* SHOULD NEVER HAPPEN -> all probes are supposed to be in the hashtable */
		{
			Rprintf("ERROR:\tCannot find probe %s in the hashtable.\nExiting...\n", str_probe);
			fclose(fp1);
		    return;
		}

		/* FIND or CREATE gene */
		gene_key temp5;
		strncpy(temp5.gene_name, match_gene, sizeof(char)*PROBE_NAMEMAX);  
		temp5.strand = str_strand[0];
		HASH_FIND(hh, *genes, &temp5, sizeof(gene_key), temp2); /* struct already in the hash? */
		if(temp2 == NULL) {
			temp2 = (gene_t*)malloc(sizeof(gene_t)); 
			temp2->id = temp5;
			temp2->best = 0;  
			temp2->best_count = 0;  
			temp2->count = 0;  
			
			HASH_ADD(hh, *genes, id, sizeof(gene_key), temp2);
		}
		
		/* FIND or CREATE score */
		score_key temp4;
		strncpy(temp4.probe_name, str_probe, sizeof(char)*PROBE_NAMEMAX); 
		temp4.gene_id = temp5;  
		HASH_FIND(hh, *scores, &temp4, sizeof(score_key), temp3); /* struct already in the hash? */
		if(temp3 == NULL) {
			temp3 = (score_t*)malloc(sizeof(score_t));
			if(temp3 == NULL) {
			  Rprintf("It appears there is not enought memory\n"); 
		  	  fclose(fp1);
		  	  return;
			}
			temp3->id = temp4;
			for(int i = 0; i < PROBE_LEN; ++i)
			  temp3->scores[i] = SCOREMAX;
			HASH_ADD(hh, *scores, id, sizeof(score_key), temp3);
			
			if(temp1->nb_genes < GENE_LEN) { 
				temp1->scores[temp1->nb_genes] = temp3;
				temp1->nb_genes++;
			}
			
		}
		
		if(temp3 == NULL) {
		  Rprintf("It appears there is not enought memory\n"); 
		  fclose(fp1);
		  return;
		
		}
		
		
// 		Rprintf(" %d - %d - %d\n", index, number, missmatch);
// 		for(int i=0; i <= 50;++i) {
// 		  int v = temp3->scores[i];
// 		   temp3->scores[i] = v;
// 		}
		//Rprintf(" %d - %d - %d\n", index, number, missmatch);
		//for(int i=T2; i <= T2;++i) {
		//  int v = temp3->scores[i];
		//  if(i == number-1) {
		//    Rprintf("%d %d %d\n", i, number-1, temp3->scores[i]);
		//    temp3->scores[i] = M2; 
		//  }
		    //temp3->scores[T2] = v;
		//}
		//temp3->scores[number-1] = missmatch;
		//if(temp3->scores[number-1] < missmatch) {
		// Rprintf("%d < %d",  temp3->scores[number-1], missmatch); 
		//}
		
		//temp3->scores[0] = 450;
		temp3->scores[number-1] = (temp3->scores[number-1] < missmatch) ?  temp3->scores[number-1] : missmatch;
		//if(0 == 1) {
		    //t
		//}


		pcre_free_substring(match_gene);

	}
}

// free connection
fclose(fp1);
}
}



void scores_get_set(char **bowtie_file, int nb_map, probe_t **probes, gene_t **genes, score_t **scores, int max_missmatch, char **pattern_gene, char **pattern_probe) {

//probe
probe_t *temp1;
score_t *temp3;
gene_t *temp2;

//connection+variables
FILE *fp1;
char line[100000];
//for strtol
char *pch;


/* PCRE+variables */
pcre *re_gene, *re_probe;
const char *pcreErrorStr;
int pcreErrorOffset;
int subStrVec[30];
int pcreExecRet;
const char *match_probe, *match_gene;
/* compile PCRE */
re_probe = pcre_compile(pattern_probe[0], 0, &pcreErrorStr, &pcreErrorOffset, NULL);
re_gene = pcre_compile(pattern_gene[0], 0, &pcreErrorStr, &pcreErrorOffset, NULL);

for(int map = 0; map < nb_map; map++) {
  Rprintf("map %d (%s)\n", map, bowtie_file[map]);
  fp1=fopen(bowtie_file[map], "r");
  int index = 0;
  while(fgets (line, 100000, fp1)!=NULL ) {
	size_t linelen = strlen(line);
	index++;
	// remove end of line character
	if (line[linelen - 1] == '\n') line[linelen - 1] = '\0';	

	char *str_probe  = strtok(line, "\t");
	char *str_strand = strtok(NULL, "\t");
	char *str_gene   = strtok(NULL, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");
	strtok(NULL, "\t");
	char *str_missmatch = strtok(NULL, "\t");
	int missmatch = strtol(str_missmatch, &pch, 10);  	
	//Rprintf("b1 %d\n", index);
		
	

	/* if missmatch is lower tan max_missmatch */ 
	if(missmatch < max_missmatch) {
		/* probe PCRE */
		pcreExecRet = pcre_exec(re_probe, NULL, str_probe, strlen(str_probe), 0, 0, subStrVec, 30);    
		if(pcreExecRet>0) {
			pcre_get_substring(str_probe, subStrVec, pcreExecRet, 0, &(match_probe));
		} else { /* SHOULD NEVER HAPPEN -> users must make sure the regex will work */
			Rprintf("ERROR:\tPattern '%s' cannot be found in '%s'.\nExiting...\n", pattern_probe[0], str_probe);
			
		  fclose(fp1);
		  return;
		}
	  //Rprintf("Probe: %s\n", str_probe);
	  //Rprintf("Probe extract: %s\n", match_probe);
	  //Rprintf("Nb: %s\n", str_probe+subStrVec[1]); 

		int number = strtol(str_probe+subStrVec[1], &pch, 10); 
		if((number < 1) || (number >= PROBE_LEN)) {
		  Rprintf("Number out of range: %d\n", number); 
		  return; 
		}
		/* gene PCRE */
		pcreExecRet = pcre_exec(re_gene, NULL, str_gene, strlen(str_gene), 0, 0, subStrVec, 30);  
		if(pcreExecRet>0) {
			pcre_get_substring(str_gene, subStrVec, pcreExecRet, 0, &(match_gene));
		} else {
			Rprintf("NO1 %d %s %d\n", number, str_gene, index);
			break;
		}
		  
		
		/* FIND probe */
		HASH_FIND_STR(*probes, match_probe, temp1);  /* name already in the hash? */
		if(temp1 == NULL) /* SHOULD NEVER HAPPEN -> all probes are supposed to be in the hashtable */
		{
			Rprintf("ERROR:\tCannot find probe %s in the hashtable.\nExiting...\n", match_probe);
			
		  fclose(fp1);
		  return;
		}

		/* FIND or CREATE gene */
		gene_key temp5;
		strncpy(temp5.gene_name, match_gene, sizeof(char)*PROBE_NAMEMAX);  
		temp5.strand = str_strand[0];
		HASH_FIND(hh, *genes, &temp5, sizeof(gene_key), temp2); /* struct already in the hash? */
		if(temp2 == NULL) {
			temp2 = (gene_t*)malloc(sizeof(gene_t)); 
			temp2->id = temp5;
			temp2->best = 0;  
			temp2->best_count = 0;  
			temp2->count = 0;  
			
			HASH_ADD(hh, *genes, id, sizeof(gene_key), temp2);
		}
		
		/* FIND or CREATE score */
		score_key temp4;
		strncpy(temp4.probe_name, match_probe, sizeof(char)*PROBE_NAMEMAX); 
		temp4.gene_id = temp5;  
		HASH_FIND(hh, *scores, &temp4, sizeof(score_key), temp3); /* struct already in the hash? */
		if(temp3 == NULL) {
			temp3 = (score_t*)malloc(sizeof(score_t));
			if(temp3 == NULL) {
		  Rprintf("It appears there is not enought memory\n"); 
		  fclose(fp1);
		  return;
			}
			temp3->id = temp4;
			for(int i = 0; i < PROBE_LEN; ++i)
			  temp3->scores[i] = SCOREMAX;
			HASH_ADD(hh, *scores, id, sizeof(score_key), temp3);
			
			if(temp1->nb_genes < GENE_LEN) { 
				temp1->scores[temp1->nb_genes] = temp3;
				temp1->nb_genes++;
			}
			
		}
		
		if(temp3 == NULL) {
		  Rprintf("It appears there is not enought memory\n"); 
		  fclose(fp1);
		  return;
		  
		}
		
		
// 		Rprintf(" %d - %d - %d\n", index, number, missmatch);
// 		for(int i=0; i <= 50;++i) {
// 		  int v = temp3->scores[i];
// 		   temp3->scores[i] = v;
// 		}
		//Rprintf(" %d - %d - %d\n", index, number, missmatch);
		//for(int i=T2; i <= T2;++i) {
		//  int v = temp3->scores[i];
		//  if(i == number-1) {
		//    Rprintf("%d %d %d\n", i, number-1, temp3->scores[i]);
		//    temp3->scores[i] = M2; 
		//  }
		    //temp3->scores[T2] = v;
		//}
		//temp3->scores[number-1] = missmatch;
		//if(temp3->scores[number-1] < missmatch) {
		// Rprintf("%d < %d",  temp3->scores[number-1], missmatch); 
		//}
		
		//temp3->scores[0] = 450;
		temp3->scores[number-1] = (temp3->scores[number-1] < missmatch) ?  temp3->scores[number-1] : missmatch;
		//if(0 == 1) {
		    //t
		//}
		
		pcre_free_substring(match_probe);
		pcre_free_substring(match_gene);


	}
	//Rprintf("e1 %d\n", index);
}
// free connection
fclose(fp1);
}
}
