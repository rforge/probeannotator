#include <R.h>
#include <pcre.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>

#define LINEMAX 600

void transform_UCSC(char **fasta_file, char **file_out) {

  
/* PCRE+variables */
pcre *re_0, *re_1, *re_2, *re_3, *re_4;
const char *pcreErrorStr;
int pcreErrorOffset;
int subStrVec[30];
int pcreExecRet;
const char *match_build, *match_gene, *match_chr, *match_loc1, *match_loc2;
/* PCRE patterns */
const char *pat0 = "hg[[:alnum:]]{2}";
const char *pat1 = "refGene_";
const char *pat2 = "chr[[:alnum:]]{1,2}";
const char *pat3 = "\\:[[:alnum:]]+";
const char *pat4 = "\\-[[:alnum:]]+";

/* compile PCRE */
re_0 = pcre_compile(pat0, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
re_1 = pcre_compile(pat1, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
re_2 = pcre_compile(pat2, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
re_3 = pcre_compile(pat3, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
re_4 = pcre_compile(pat4, 0, &pcreErrorStr, &pcreErrorOffset, NULL);

  
/* connection+variables */
FILE *fp1;
char line[LINEMAX];
/* for strtol */
char *pch;
int max = 0;
fp1=fopen(fasta_file[0], "r");
Rprintf("Build\tTranscript\tChromosome\tSeq_Start\tSeq_End\tStrand\n");
		
//example of line: >hg19_refGene_NM_032291 range=chr1:66998325-66999824 5'pad=0 3'pad=0 strand=+ repeatMasking=none
while(fgets (line, LINEMAX, fp1)!=NULL ) {
	size_t linelen = strlen(line);
	if(linelen > 0 && line[0] == '>') {
		/* remove end of line character */
		if (line[linelen - 1] == '\n') line[linelen - 1] = '\0';	
		
		char *str_gene  = strtok(line, " ");
		char *str_loc = strtok(NULL, " ");
		strtok(NULL, " ");
		strtok(NULL, " ");
		char *str_strand   = strtok(NULL, " ");
		
		/* 0: build */
		pcreExecRet = pcre_exec(re_0, NULL, str_gene, strlen(str_gene), 0, 0, subStrVec, 30); 
		pcre_get_substring(str_gene, subStrVec, pcreExecRet, 0, &(match_build));
		/* 1: gene */
		pcreExecRet = pcre_exec(re_1, NULL, str_gene, strlen(str_gene), 0, 0, subStrVec, 30); 
		pcre_get_substring(str_gene, subStrVec, pcreExecRet, 0, &(match_gene));
		int gene_start = subStrVec[1];
		
		/* 2: chr */
		pcreExecRet = pcre_exec(re_2, NULL, str_loc, strlen(str_loc), 0, 0, subStrVec, 30); 
		pcre_get_substring(str_loc, subStrVec, pcreExecRet, 0, &(match_chr));
		/* 3: loc1 */
		pcreExecRet = pcre_exec(re_3, NULL, str_loc, strlen(str_loc), 0, 0, subStrVec, 30); 
		pcre_get_substring(str_loc, subStrVec, pcreExecRet, 0, &(match_loc1));
		/* 4: loc2 */
		pcreExecRet = pcre_exec(re_4, NULL, str_loc, strlen(str_loc), 0, 0, subStrVec, 30); 
		pcre_get_substring(str_loc, subStrVec, pcreExecRet, 0, &(match_loc2));
			
		/* : strand */
		char strand = str_strand[strlen(str_strand)-1];
		
		
		Rprintf("%s\t%s\t%s\t%s\t%s\t%c\n", match_build, str_gene+gene_start, match_chr, match_loc1+1, match_loc2+1, strand);
		
		
		pcre_free_substring(match_build);
		pcre_free_substring(match_gene);

		pcre_free_substring(match_chr);
		pcre_free_substring(match_loc1);
		pcre_free_substring(match_loc2);

	max++;
	if(max == 10) {break;}
	}
}
/* free connection */
fclose(fp1);
/* free PCRE */
free(re_0);
free(re_1);
free(re_2);
free(re_3);
free(re_4);


}

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

