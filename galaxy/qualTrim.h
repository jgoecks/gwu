/*
  John Gaspar
  July 2014

  Header file for qualTrim.c.
*/

#define MAX_SIZE    1024   // maximum length for input line

#define HELP        "-h"
#define INFILE      "-i"
#define OUTFILE     "-o"
#define WINDOWLEN   "-l"   // length for sliding window of quality scores
#define WINDOWAVG   "-q"   // average quality score for sliding window
#define QUALAVG     "-t"   // average quality score for entire read
#define MINLEN      "-n"   // minimum length of a read

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
#define ERRFLOAT    6
#define MERRFLOAT   ": cannot convert to float"
#define ERRINT      7
#define MERRINT     ": cannot convert to int"
#define ERRLEN      8
#define MERRLEN     ": sequence too short for quality window"
