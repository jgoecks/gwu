/*
  John Gaspar
  June 2014

  Header file for removeDups.c.
*/

#define MAX_SIZE    1024   // maximum length for input line
#define HEADER      60     // maximum header length

#define HELP        "-h"
#define INFILE      "-i"
#define OUTFILE     "-o"
#define DUPFILE     "-d"
#define CSV         ","
#define CHILD       5      // number of children nodes

#define ERROPEN     0
#define MERROPEN    "cannot open file for reading"
#define ERRCLOSE    1
#define MERRCLOSE   "cannot close file"
#define ERROPENW    2
#define MERROPENW   "cannot open file for writing"
#define ERRUNK      3
#define MERRUNK     "unknown file type (not fasta or fastq)"
#define ERRMEM      4
#define MERRMEM     "cannot allocate memory"
#define ERRSEQ      5
#define MERRSEQ     "cannot load sequence"

/*typedef struct read {
  char* header;
  struct read* next;
} Read;*/

typedef struct node {
  char* seq;
  uint16_t st;
  uint16_t end;
  //Read* first;
  char* head;
  //struct node* next;
  struct node* child[CHILD];
} Node;
