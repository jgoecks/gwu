/*
  John Gaspar
  July 2014

  Header file for removePrimer.c.
*/

#define MAX_SIZE    1024   // maximum length for input line
#define HEADER      60     // maximum header length

#define HELP        "-h"
#define INFILE      "-i"
#define OUTFILE     "-o"
#define PRIMFILE    "-p"
#define LOGFILE     "-l"
#define WASTEFILE   "-w"
#define MISALLOW    "-e"
#define REVLENGTH   "-r"
#define REVMIS      "-er"
#define REVOPT      "-rq"
#define CORRFILE    "-c"
#define CSV         ",\t"
#define DEL         "\n"

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
#define ERRPRIM     6
#define MERRPRIM    "cannot load primer sequence"
#define ERRPREP     7
#define MERRPREP    ": cannot repeat primer name"
#define ERRINT      8
#define MERRINT     ": cannot convert to int"

typedef struct primer {
  char* name;
  char* fwd;
  char* rev;
  char* frc;
  char* rrc;
  int fcount;
  int rcount;
  int fcountr;
  int rcountr;
  struct primer* next;
} Primer;
