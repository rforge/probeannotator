#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "uthash.h"


#define LINEMAX 600
#define SCOREMAX 500
#define PROBE_NAMEMAX 40
#define PROBE_LEN 500
#define GENE_LEN 200

typedef struct gene_key_stu {           
	char gene_name[PROBE_NAMEMAX];    
	char strand;
	UT_hash_handle hh; /* makes this structure hashable */
} gene_key;


typedef struct score_key_stu {
	char probe_name[PROBE_NAMEMAX];             
	gene_key gene_id;      
	UT_hash_handle hh; /* makes this structure hashable */
} score_key;


typedef struct score_stu {
	score_key id;		
	int scores[PROBE_LEN];	
	UT_hash_handle hh;		/* makes this structure hashable */
} score_t;	

typedef struct gene_stu {
    gene_key id;             /* key (string is WITHIN the structure) == gene name*/
    double best;
    int best_count;
    int count;
    UT_hash_handle hh;         /* makes this structure hashable */
} gene_t;



typedef struct probe_stu {
    char name[PROBE_NAMEMAX];             			/* key (string is WITHIN the structure) == probe name */
    int nb_probes;		/* nb of probes in probeset */
    int nb_genes;		/* number of genes pointed */
    double best_score;
    int  best_index;
    score_t *best_score_s[1];			/* pointer to gene structure */
    score_t *scores[GENE_LEN];			/* pointer to gene structure */
    UT_hash_handle hh;         			/* makes this structure hashable */
} probe_t;





