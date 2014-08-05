#include <R.h>
#include <pcre.h>
#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include "structs.h"


#define LINEMAX 600

void transform_TAB(char **fasta_file, char **file_out, char **scan_pattern, int* args) {
 //args=str_first? + has_probes?
  
 Rprintf("|Args: %d-%d\n|Scan Pattern: %s\n", args[0], args[1], scan_pattern[0]);

/* connection+variables */
FILE *fp1,*fp2;
char line[LINEMAX];

char *pch;
int max = 0;
fp1=fopen(fasta_file[0], "r");
fp2=fopen(file_out[0], "w");

char str_probe[PROBE_NAMEMAX];
char str_sequence[500];

// read header
if(fgets (line, LINEMAX, fp1)==NULL)
{
  return;
}


if(args[1] == 0) {
// - considering the index of 'str_probe' < index of 'str_sequence'
    if(args[0] == 0) {
      while(fgets (line, LINEMAX, fp1)!=NULL ) {
	      size_t linelen = strlen(line);
	      sscanf(line, scan_pattern[0], str_probe, str_sequence);
	      fprintf(fp2, ">%s\n%s\n", str_probe, str_sequence);
      }
    } else {
      // - considering the index of 'str_probe' > index of 'str_sequence'
	while(fgets (line, LINEMAX, fp1)!=NULL ) {
	      size_t linelen = strlen(line);
	      sscanf(line, scan_pattern[0], str_sequence, str_probe);
	      fprintf(fp2, ">%s\n%s\n", str_probe, str_sequence);
      }
    }
  } else {
    char str_lastProbe[PROBE_NAMEMAX]; 
    int probeNumber = 1; 
   
    if(args[0] == 0) {
      while(fgets (line, LINEMAX, fp1)!=NULL ) {
	      size_t linelen = strlen(line);
	      sscanf(line, scan_pattern[0], str_probe, str_sequence);
	      
	      if(strcmp(str_lastProbe, str_probe) == 0) {
		probeNumber++;
	      } else {
		probeNumber=1;
	      }
	      fprintf(fp2, ">%s-%d\n%s\n", str_probe, probeNumber, str_sequence);
	      
	      strncpy(str_lastProbe, str_probe, PROBE_NAMEMAX);
      }
    } else {
      // - considering the index of 'str_probe' > index of 'str_sequence'
	while(fgets (line, LINEMAX, fp1)!=NULL ) {
	      size_t linelen = strlen(line);
	      sscanf(line, scan_pattern[0], str_sequence, str_probe);

	      if(strcmp(str_lastProbe, str_probe) == 0) {
		probeNumber++;
	      } else {
		probeNumber=1;
	      }
	      fprintf(fp2, ">%s-%d\n%s\n", str_probe, probeNumber, str_sequence);
	      
	      strncpy(str_lastProbe, str_probe, PROBE_NAMEMAX);
      }
    }

  }

/* free connection */
fclose(fp1);
fclose(fp2);


}
