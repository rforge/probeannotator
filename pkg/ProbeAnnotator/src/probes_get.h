#include <R.h>
#ifdef HAVE_PCRE_H
#include <pcre.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "structs.h"

//R CMD SHLIB --preclean -I./ get_score.c

void get_probes(char **fasta_file, probe_t **probes) {

/* STEP 0 - initialise variables */
/* - probe hashtable (temp) */
probe_t *temp;
/* - file connection and associated variables */
FILE *fp1;
char line[5000];

/* STEP 1 - open fasta file */
fp1=fopen(fasta_file[0], "r");


while(fgets (line, 5000, fp1)!=NULL ) {
	size_t linelen = strlen(line);
	/*IF: line is probe name when it begins by '>' in fasta files */	
	if(linelen > 0 && line[0] == '>') {
	
	/*FORMAT: remove end of line character */
	if (line[linelen - 1] == '\n') line[linelen - 1] = '\0';

  	/* copy to sized buffer (overflow safe): */
	char name[PROBE_NAMEMAX];
  	strncpy ( name, line+1, sizeof(name) );	
	
  	/* STEP x - initialise probe's struct */
	temp = (probe_t*)malloc(sizeof(probe_t));
 	
  	/* copy to sized buffer (overflow safe): */
  	strncpy ( temp->name, line+1, sizeof(temp->name) );
	temp->nb_probes = 1;
	temp->nb_genes = 0;
	//for(int i = 0; i < PROBE_LEN; ++i)
	//  temp->scores[i] = 0;
	
  	/* STEP x - add probe's struct */
	HASH_ADD_STR( *probes, name, temp );  /* id: name of key field */
	#ifdef VERBOSE
	Rprintf("%s\n",temp->name);
	#endif
	}
}

/* STEP x - free connection */
fclose(fp1);
}




void get_probes_set(char **fasta_file, char **pattern, probe_t **probes) {
#ifndef HAVE_PCRE_H
Rprintf("PCRE library not found during compilation. Function empty");
return;
#endif
  
#ifdef HAVE_PCRE_H
//probes
probe_t *temp;

//PCRE+variables
pcre *reCompiled;
const char *pcreErrorStr;
int pcreErrorOffset;
int subStrVec[30];
int pcreExecRet;

//connection+variables
FILE *fp1;
char line[5000];

//open connection
fp1=fopen(*fasta_file, "r");
//compile pattern
reCompiled = pcre_compile(pattern[0], 0, &pcreErrorStr, &pcreErrorOffset, NULL);
int index = 0;
while(fgets (line, 5000, fp1)!=NULL ) {
	++index;
	
	size_t linelen = strlen(line);
	
	if(line[0] == '>') {
	
	//remove end of line character
	if (line[linelen - 1] == '\n') line[linelen - 1] = '\0';
	//Rprintf(".l ");

	//Rprintf("%s\n", line);
	pcreExecRet = pcre_exec(reCompiled, NULL, line+1, strlen(line+1), 0, 0, subStrVec, 30);    
	
	// get first match
	
 	const char *psubStrMatchStr;
 	if(pcreExecRet>0) {
 		pcre_get_substring(line+1, subStrVec, pcreExecRet, 0, &(psubStrMatchStr));
 	} else 
	{
		Rprintf("ERROR:\tPattern '%s' cannot be found in '%s'.\nExiting...\n", pattern[0], line);
		fclose(fp1);
		return;
		//exit(0);
	}


  	/* copy to sized buffer (overflow safe): */
	char name[PROBE_NAMEMAX];
  	strncpy ( name, psubStrMatchStr, sizeof(name) );	
	
	
	#ifdef VERBOSE
	Rprintf("%s\n",name);
	#endif
	HASH_FIND_STR(*probes, name, temp);  /* name already in the hash? */
     	if (temp==NULL) {
 	  	/* create probe struct */
 		temp = (probe_t*)malloc(sizeof(probe_t));
 	  	/* copy to sized buffer (overflow safe): */
 	  	strncpy ( temp->name, psubStrMatchStr, sizeof(temp->name) ); 
 		temp->nb_probes = 1; 
 		temp->nb_genes = 0;
 		HASH_ADD_STR( *probes, name, temp );  /* id: name of key field */

 		#ifdef VERBOSE
 		Rprintf("%s\n",temp->name);
 		#endif
 	} else {
 		temp->nb_probes++;
 	}
	
	/* free the substring */
	pcre_free_substring(psubStrMatchStr);
	}
}

/* free connection */
fclose(fp1);
#endif
}
