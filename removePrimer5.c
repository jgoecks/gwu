/*
  John Gaspar
  July 2014

  Removing primers from a fasta/fastq file.
  Version 5: reattaches correct primers.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "removePrimer5.h"

// globals
static char* line;
static char* hline;
static Primer* primo;

// usage
void usage(void) {
  fprintf(stderr, "Usage: ./removePrimer  %s <input file>  ", INFILE);
  fprintf(stderr, "%s <output file>  ", OUTFILE);
  fprintf(stderr, "%s <primer file>\n", PRIMFILE);
  fprintf(stderr, "Optional commands:\n");
  fprintf(stderr, "  %s  <file>   Log file for counts\n", LOGFILE);
  fprintf(stderr, "  %s  <file>   File containing non-matched reads\n", WASTEFILE);
  fprintf(stderr, "  %s  <int>    Mismatches to the forward primer to allow (def. 0)\n", MISALLOW);
  fprintf(stderr, "  %s  <int>    Length of reverse primer to check (def. full length)\n", REVLENGTH);
  fprintf(stderr, "  %s <int>    Mismatches to the reverse primer to allow (def. 0)\n", REVMIS);
  fprintf(stderr, "  %s          Option to require reverse primer be found\n", REVOPT);
  fprintf(stderr, "  %s  <file>   Output file containing reads with corrected primers attached\n", CORRFILE);
  fprintf(stderr, "                 (should only be used if specifying %s)\n", REVOPT);
  exit(-1);
}

// error
int error(char* msg, int err) {
  char* msg2;
  if (err == ERROPEN)
    msg2 = MERROPEN;
  else if (err == ERRCLOSE)
    msg2 = MERRCLOSE;
  else if (err == ERROPENW)
    msg2 = MERROPENW;
  else if (err == ERRUNK)
    msg2 = MERRUNK;
  else if (err == ERRMEM)
    msg2 = MERRMEM;
  else if (err == ERRSEQ)
    msg2 = MERRSEQ;
  else if (err == ERRPRIM)
    msg2 = MERRPRIM;
  else if (err == ERRPREP)
    msg2 = MERRPREP;
  else if (err == ERRINT)
    msg2 = MERRINT;

  fprintf(stderr, "Error! %s: %s\n", msg, msg2);
  return -1;
}

// freeMemory
void freeMemory(void) {
  Primer* temp;
  for (Primer* p = primo; p != NULL; ) {
    free(p->name);
    free(p->fwd);
    free(p->rev);
    free(p->frc);
    free(p->rrc);
    temp = p;
    p = p->next;
    free(temp);
  }
  free(line);
  free(hline);
}

// memalloc
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char** endptr = NULL;
  int ans = (int) strtol(in, endptr, 10);
  if (endptr != '\0')
    exit(error(in, ERRINT));
  return ans;
}

// fastaOrQ -- is file fasta or fastq
int fastaOrQ(FILE* in) {
  while (fgets(line, MAX_SIZE, in) != NULL)
    if (line[0] == '>')
      return 0;
    else if (line[0] == '@')
      return 1;
    else if (line[0] != '#')
      break;
  exit(error("", ERRUNK));
}

/* int ambig(char, char)
 * Checks ambiguous DNA bases.
 */
int ambig(char x, char y) {
  if (x == 'N' ||
      (x == 'W' && (y == 'A' || y == 'T')) ||
      (x == 'S' && (y == 'C' || y == 'G')) ||
      (x == 'M' && (y == 'A' || y == 'C')) ||
      (x == 'K' && (y == 'G' || y == 'T')) ||
      (x == 'R' && (y == 'A' || y == 'G')) ||
      (x == 'Y' && (y == 'C' || y == 'T')) ||
      (x == 'B' && (y == 'C' || y == 'G' || y == 'T')) ||
      (x == 'D' && (y == 'A' || y == 'G' || y == 'T')) ||
      (x == 'H' && (y == 'A' || y == 'C' || y == 'T')) ||
      (x == 'V' && (y == 'A' || y == 'C' || y == 'G')))
    return 0;
  return 1;
}

/* Primer* findPrim(char*)
 * Finds a primer match to the given sequence.
 */
Primer* findPrim(char* seq, int misAllow, int* st, int* f) {
  for (Primer* p = primo; p != NULL; p = p->next) {
    for (int i = 0; i < 2; i++) {
      char* prim = (i ? p->rrc : p->fwd);
      int mis = misAllow;
      int j;
      for (j = 0; prim[j] != '\0'; j++)
        if (seq[j] == '\0' || (prim[j] != seq[j] &&
            (prim[j] == 'A' || prim[j] == 'C' ||
            prim[j] == 'G' || prim[j] == 'T' ||
            ambig(prim[j], seq[j])) && --mis < 0))
          break;

      if (prim[j] == '\0') {
        if (seq[j] != '\0') {
          *st = j;
          i ? p->rcount++ : p->fcount++;
          *f = i;
          return p;
        } else
          return NULL;
      }
    }
  }
  return NULL;
}

/* int checkRev()
 * Checks the seq for a match of the reverse primer.
 */
int checkRev(char* seq, int st, Primer* p, int f,
    int misAllow, int revLen) {
  char* rev = (f ? p->frc : p->rev);
  int last = strlen(seq) - revLen + 1;
  for (int i = st + 1; i < last; i++) {
    int mis = misAllow;
    int j;
    for (j = 0; j < revLen; j++)
      if (rev[j] != seq[i + j] &&
          (rev[j] == 'A' || rev[j] == 'C' ||
          rev[j] == 'G' || rev[j] == 'T' ||
          ambig(rev[j], seq[i + j])) && --mis < 0)
        break;
    if (j == revLen) {
      f ? p->rcountr++ : p->fcountr++;
      return i;
    }
  }
  return 0;
}

/* int readFile()
 * Parses the input file. Produces the output file(s).
 */
int readFile(FILE* in, FILE* out, int misAllow, int* match,
    int* rcmatch, FILE* waste, int revLen, int revMis,
    int revOpt, FILE* corr) {
  int aorq = fastaOrQ(in);
  rewind(in);
  int count = 0;
  while (fgets(hline, HEADER, in) != NULL) {
    count++;
    if (fgets(line, MAX_SIZE, in) == NULL)
      exit(error("", ERRSEQ));
    int len = strlen(line) - 1;
    if (line[len] == '\n')
      line[len] = '\0';

    int st = 0, end = 0, f = 0;
    Primer* p = findPrim(line, misAllow, &st, &f);
    if (p != NULL) {
      (*match)++;
      int rcLen = (f ? strlen(p->frc) : strlen(p->rev));
      if (revLen && rcLen > revLen)
        rcLen = revLen;
      end = checkRev(line, st, p, f, revMis, rcLen); // search for reverse primer
      if (revOpt && !end) {
        if (waste != NULL)
          fprintf(waste, "%s%s\n", hline, line);
      } else {
        // print header
        for (int i = 0; hline[i] != '\0' && hline[i] != '\n'; i++)
          fprintf(out, "%c", hline[i]);
        fprintf(out, " %s%s\n", p->name, end ? " rev" : "");
        if (corr != NULL) {
          for (int i = 0; hline[i] != '\0' && hline[i] != '\n'; i++)
            fprintf(corr, "%c", hline[i]);
          fprintf(corr, " %s%s\n", p->name, end ? " rev" : "");
        }
        // print sequence
        if (!end)
          end = len;
        else
          (*rcmatch)++;
        for (int i = st; i < end; i++)
          fprintf(out, "%c", line[i]);
        fprintf(out, "\n");
        if (corr != NULL) {
          // reattach primers
          fprintf(corr, "%s", f ? p->rrc : p->fwd);
          for (int i = st; i < end; i++)
            fprintf(corr, "%c", line[i]);
          fprintf(corr, "%s\n", f ? p->frc : p->rev);
        }
      }
    } else if (waste != NULL)
      fprintf(waste, "%s%s\n", hline, line);

    // read next 2 lines if fastq
    if (aorq) {
      for (int i = 0; i < 2; i++)
        if (fgets(line, MAX_SIZE, in) == NULL)
          exit(error("", ERRSEQ));
        else if (p != NULL) {
          if (revOpt && !end) {
            if (waste != NULL)
              fprintf(waste, "%s", line);
          } else if (i) {
            for (int j = st; j < end; j++)
              fprintf(out, "%c", line[j]);
            fprintf(out, "\n");
            if (corr != NULL) {
              for (int j = 0; j < strlen(f ? p->rrc : p->fwd); j++)
                fprintf(corr, "I");
              for (int j = st; j < end; j++)
                fprintf(corr, "%c", line[j]);
              for (int j = 0; j < strlen(f ? p->frc : p->rev); j++)
                fprintf(corr, "I");
              fprintf(corr, "\n");
            }
          } else {
            fprintf(out, "%s", line);
            if (corr != NULL)
              fprintf(corr, "%s", line);
          }
        } else if (waste != NULL)
          fprintf(waste, "%s", line);
    }
  }
  //printf("  reads: %d\n", count);
  return count;
}

// openWrite
FILE* openWrite(char* outFile) {
  FILE* out = fopen(outFile, "r");
  if (out != NULL) {
    if (fclose(out))
      exit(error("", ERRCLOSE));
    exit(error(outFile, ERROPENW));
  } else
    out = fopen(outFile, "w");
  if (out == NULL)
    exit(error(outFile, ERROPENW));
  return out;
}

// openFiles
void openFiles(char* outFile, FILE** out,
    char* primFile, FILE** prim, char* inFile, FILE** in,
    char* logFile, FILE** log, char* wasteFile, FILE** waste,
    char* corrFile, FILE** corr) {
  *out = openWrite(outFile);
  *prim = fopen(primFile, "r");
  if (*prim == NULL)
    exit(error(primFile, ERROPEN));
  *in = fopen(inFile, "r");
  if (*in == NULL)
    exit(error(inFile, ERROPEN));
  if (logFile != NULL)
    *log = openWrite(logFile);
  if (wasteFile != NULL)
    *waste = openWrite(wasteFile);
  if (corrFile != NULL)
    *corr = openWrite(corrFile);
}

/* char rc(char)
 * Returns the complement of the given base.
 */
char rc(char in) {
  char out;
  if (in == 'A')
    out = 'T';
  else if (in == 'T')
    out = 'A';
  else if (in == 'C')
    out = 'G';
  else if (in == 'G')
    out = 'C';
  else if (in == 'Y')
    out = 'R';
  else if (in == 'R')
    out = 'Y';
  else if (in == 'W')
    out = 'W';
  else if (in == 'S')
    out = 'S';
  else if (in == 'K')
    out = 'M';
  else if (in == 'M')
    out = 'K';
  else if (in == 'B')
    out = 'V';
  else if (in == 'V')
    out = 'B';
  else if (in == 'D')
    out = 'H';
  else if (in == 'H')
    out = 'D';
  else if (in == 'N')
    out = 'N';
  else
    exit(error("", ERRPRIM));
  return out;
}

/* char* revComp(char*)
 * Reverse-complements the given sequence.
 */
char* revComp(char* seq) {
  int i = strlen(seq) - 1;
  char* out = (char*) memalloc(2 + i);
  int j;
  for (j = 0; i > -1; j++) {
    char nuc = seq[i--];
    out[j] = rc(nuc);
  }
  out[j] = '\0';
  return out;
}

/* int loadSeqs(FILE*)
 * Loads the primers from the given file.
 */
int loadSeqs(FILE* prim) {

  Primer* prev = NULL;
  int count = 0;
  while (fgets(line, MAX_SIZE, prim) != NULL) {

    if (line[0] == '#')
      continue;

    // load name and sequence
    char* name = strtok(line, CSV);
    char* seq = strtok(NULL, CSV);
    char* rev = strtok(NULL, DEL);
    if (name == NULL || seq == NULL || rev == NULL) {
      error("", ERRPRIM);
      continue;
    }

    // check for duplicate
    for (Primer* pc = primo; pc != NULL; pc = pc->next)
      if (!strcmp(pc->name, name))
        exit(error(name, ERRPREP));

    // create primer
    Primer* p = (Primer*) memalloc(sizeof(Primer));
    p->name = (char*) memalloc(1 + strlen(name));
    p->fwd = (char*) memalloc(1 + strlen(seq));
    p->rev = (char*) memalloc(1 + strlen(rev));
    strcpy(p->name, name);
    strcpy(p->fwd, seq);
    strcpy(p->rev, rev);

    // save sequence rc's
    p->frc = revComp(p->fwd);
    p->rrc = revComp(p->rev);

    p->fcount = p->rcount = p->fcountr = p->rcountr = 0;
    p->next = NULL;
    if (primo == NULL)
      primo = p;
    else
      prev->next = p;
    prev = p;
    count++;
  }

  return count;
}


// getParams
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile = NULL, *primFile = NULL,
    *logFile = NULL, *wasteFile = NULL, *corrFile = NULL;
  int misAllow = 0, revLen = 0, revMis = 0, revOpt = 0;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (!strcmp(argv[i], REVOPT))
      revOpt = 1;
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], INFILE))
        inFile = argv[++i];
      else if (!strcmp(argv[i], PRIMFILE))
        primFile = argv[++i];
      else if (!strcmp(argv[i], LOGFILE))
        logFile = argv[++i];
      else if (!strcmp(argv[i], WASTEFILE))
        wasteFile = argv[++i];
      else if (!strcmp(argv[i], MISALLOW))
        misAllow = getInt(argv[++i]);
      else if (!strcmp(argv[i], REVLENGTH))
        revLen = getInt(argv[++i]);
      else if (!strcmp(argv[i], REVMIS))
        revMis = getInt(argv[++i]);
      else if (!strcmp(argv[i], CORRFILE))
        corrFile = argv[++i];
    }
    else
      usage();
  }

  if (outFile == NULL || inFile == NULL || primFile == NULL)
    usage();

  FILE* out = NULL, *prim = NULL, *in = NULL, *log = NULL,
    *waste = NULL, *corr = NULL;
  openFiles(outFile, &out, primFile, &prim, inFile, &in,
    logFile, &log, wasteFile, &waste, corrFile, &corr);
  int pr = loadSeqs(prim);

  int match = 0, rcmatch = 0;
  int count = readFile(in, out, misAllow, &match, &rcmatch,
    waste, revLen, revMis, revOpt, corr);
  if (log != NULL) {
    fprintf(log, "Primer pairs: %d\nRead count: %d\n", pr, count);
    fprintf(log, "Primer matches: %d\nRev primer matches: %d\n", match, rcmatch);
    for (Primer* p = primo; p != NULL; p = p->next)
      fprintf(log, "%s\t%d\t%d\t%d\t%d\n", p->name, p->fcount, p->rcount,
        p->fcountr, p->rcountr);
  }

  if (fclose(out) || fclose(prim) || fclose(in) ||
      (log != NULL && fclose(log)) || (waste != NULL && fclose(waste))
      || (corr != NULL && fclose(corr)))
    exit(error("", ERRCLOSE));
}

// main
int main(int argc, char* argv[]) {
  line = (char*) memalloc(MAX_SIZE);
  hline = (char*) memalloc(HEADER);
  primo = NULL;
  getParams(argc, argv);
  freeMemory();
  return 0;
}
