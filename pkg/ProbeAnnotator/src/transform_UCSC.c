#include <R.h>
#ifdef HAVE_PCRE_H
//#include <pcre.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stddef.h>

#define LINEMAX 600


void get_UCSC_refSeq_count(char **fasta_file, int *result_len) {

    
/* connection+variables */
FILE *fp1;
char line[LINEMAX];

/* counter */
int counter = 0;


fp1=fopen(fasta_file[0], "r");
while(fgets (line, LINEMAX, fp1)!=NULL ) {
	//counter += (line[0] == '>') ? 1 : 0;
	(line[0] == '>') ? counter++ : 0;
}
/* free connection */
fclose(fp1);

/* set result */
result_len[0] = counter;

}

void get_UCSC_refSeq(char **fasta_file, char **result, int *result_len) {
#ifndef HAVE_PCRE_H
Rprintf("PCRE library not found during compilation. Function empty");
return;
#endif
  
#ifdef HAVE_PCRE_H
  
/* PCRE+variables */
pcre *re_1;
const char *pcreErrorStr;
int pcreErrorOffset;
int subStrVec[2];
int pcreExecRet;
const char *match_gene;
/* PCRE patterns */
const char *pat1 = "refGene_";

/* compile PCRE */
re_1 = pcre_compile(pat1, 0, &pcreErrorStr, &pcreErrorOffset, NULL);

  
/* connection+variables */
FILE *fp1;
char line[LINEMAX];


int counter = 0;
fp1=fopen(fasta_file[0], "r");


//example of line: >hg19_refGene_NM_032291 range=chr1:66998325-66999824 5'pad=0 3'pad=0 strand=+ repeatMasking=none
while(fgets (line, LINEMAX, fp1)!=NULL ) {
	size_t linelen = strlen(line);
	if(linelen > 0 && line[0] == '>') {
		/* remove end of line character */
		if (line[linelen - 1] == '\n') line[linelen - 1] = '\0';	
		
		char *str_gene  = strtok(line, " ");
		
		/* 1: gene */
		pcreExecRet = pcre_exec(re_1, NULL, str_gene, strlen(str_gene), 0, 0, subStrVec, 2); 
		pcre_get_substring(str_gene, subStrVec, pcreExecRet, 0, &(match_gene));
		int gene_start = subStrVec[1];
		
		
		strncpy(result[counter], str_gene+gene_start, 30);
		pcre_free_substring(match_gene);
		counter++;
		if(counter >= result_len[0]) break;
		
	}
}
/* free connection */
fclose(fp1);
/* free PCRE */
free(re_1);
#endif
}

void print_transformed_UCSC(char **fasta_file, char **file_out, char **new_names) {
  
/* connection+variables */
FILE *fp1, *fp2;
char line[LINEMAX];

/* counter */
int counter = -1;

char *p1 = ">%s\n";
char *p2 = "%.0s%s";
char *p3;


fp1=fopen(fasta_file[0], "r");
fp2=fopen(file_out[0], "w");
while(fgets (line, LINEMAX, fp1)!=NULL ) {
	
	(line[0] == '>') ? counter++, (p3 = p1) : (p3 = p2);	
	fprintf(fp2, p3, new_names[counter], line);
	/*if(line[0] == '>') {
		//Rprintf(">%s\n",new_names[counter]);	
		fprintf (fp2, ">%s\n",new_names[counter]);	
		counter++;
	} else {
		fprintf (fp2, "%s",line);	
	}*/
}
/* free connection */
fclose(fp1);
fclose(fp2);


/* set result */
}

