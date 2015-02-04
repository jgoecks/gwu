/*
  John Gaspar
  July 2014

  Quality trimming a fastq file at the ends.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qualTrim.h"

/* void usage()
 * Print usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./qualTrim  %s <input file>  ", INFILE);
  fprintf(stderr, "%s <output file>", OUTFILE);
  fprintf(stderr, "\nOptional parameters:\n");
  fprintf(stderr, "  %s <int>    Window length\n", WINDOWLEN);
  fprintf(stderr, "  %s <float>  Minimum avg. quality in the window\n", WINDOWAVG);
  fprintf(stderr, "  %s <float>  Minimum avg. quality for the full read\n", QUALAVG);
  fprintf(stderr, "                (after any window truncations)\n");
  fprintf(stderr, "  %s <int>    Minimum length of a read\n", MINLEN);
  fprintf(stderr, "                (after any window truncations)\n");
  exit(-1);
}

/* int error()
 * Print an error message.
 */
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
  else if (err == ERRINT)
    msg2 = MERRINT;
  else if (err == ERRFLOAT)
    msg2 = MERRFLOAT;
  else if (err == ERRLEN)
    msg2 = MERRLEN;

  fprintf(stderr, "Error! %s: %s\n", msg, msg2);
  return -1;
}

/* void memalloc()
 * Allocate memory from the heap.
 */
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* float getFloat(char*)
 * Converts the given char* to a float.
 */
float getFloat(char* in) {
  char** endptr = NULL;
  float ans = strtof(in, endptr);
  if (endptr != '\0')
    exit(error(in, ERRFLOAT));
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

/* int getEnd()
 * Determine the 3' end.
 */
int getEnd(char* line, int len, float qual, int end) {
  int last = 0;
  int i;
  for (i = 1; i < len; i++) {
    float sum = 0.0f;
    int j = end - 1;
    for (int k = 0; k < i; k++)
      sum += line[j - k] - 33;
    if (sum / i < qual)
      last = j - i + 1;
    else if (last)
      return last;
  }
  for (i = end - 1; i > len - 2; i--) {
    float sum = 0.0f;
    for (int j = 0; j < len; j++)
      sum += line[i - j] - 33;
    if (sum / len < qual)
      last = i - len + 1;
    else if (last)
      return last;
    else
      break;
  }
  return i == len - 2 ? last : end;
}

/* int getStart()
 * Determine the 5' end.
 */
int getStart(char* line, int len, float qual, int end) {
  int st = 0;
  int i;
  for (i = 1; i < len; i++) {
    float sum = 0.0f;
    for (int k = 0; k < i; k++)
      sum += line[k] - 33;
    if (sum / i < qual)
      st = i;
    else if (st)
      return st;
  }
  for (i = 0; i < end - len + 1; i++) {
    float sum = 0.0f;
    for (int j = 0; j < len; j++)
      sum += line[i + j] - 33;
    if (sum / len < qual)
      st = i + len;
    else if (st)
      return st;
    else
      break;
  }
  return st;
}

/* int trimQual()
 * Determine ends of the read.
 */
int trimQual(char* line, int len, float qual, int* end) {
  *end = getEnd(line, len, qual, *end);
  int st = getStart(line, len, qual, *end);
  return st;
}

/* int checkQual()
 * Check average of all qual scores (from st to end).
 * Return 1 if OK, else 0.
 */
int checkQual(char* line, int st, int end, float avg) {
  float sum = 0.0f;
  for (int i = st; i < end; i++)
    sum += line[i] - 33;
  return sum / (end - st) < avg ? 1 : 0;
}

/* void readFile()
 * Control the I/O.
 */
void readFile(FILE* in, FILE* out, int len, float qual,
    float avg, int minLen) {
  char* head = (char*) memalloc(MAX_SIZE);
  char* seq = (char*) memalloc(MAX_SIZE);
  char* line = (char*) memalloc(MAX_SIZE);

  int count = 0, elim = 0;
  while (fgets(head, MAX_SIZE, in) != NULL) {
    if (head[0] != '@')
      continue;

    if (fgets(seq, MAX_SIZE, in) == NULL)
      exit(error("", ERRSEQ));

    for (int i = 0; i < 2; i++)
      if (fgets(line, MAX_SIZE, in) == NULL)
        exit(error("", ERRSEQ));

    int end = strlen(line) - 1;
    if (line[end] == '\n')
      line[end] = '\0';
    if (end < len) {
      elim++;
      continue;
    }

    int st = 0;
    if (len)
      st = trimQual(line, len, qual, &end);
    if (avg && checkQual(line, st, end, avg)) {
      elim++;
      continue;
    }

    if (st < end && end - st >= minLen) {
      fprintf(out, "%s", head);
      for (int i = st; i < end; i++)
        fprintf(out, "%c", seq[i]);
      fprintf(out, "\n+\n");
      for (int i = st; i < end; i++)
        fprintf(out, "%c", line[i]);
      fprintf(out, "\n");
      count++;
    } else
      elim++;
  }
  printf("Reads printed: %d\nReads eliminated: %d\n", count, elim);

  free(seq);
  free(line);
  free(head);
}

/* FILE* openWrite()
 * Open a file for writing.
 */
FILE* openWrite(char* outFile) {
  FILE* out = fopen(outFile, "w");
  if (out == NULL)
    exit(error(outFile, ERROPENW));
  return out;
}

/* void openFiles()
 * Open input and output files.
 */
void openFiles(char* outFile, FILE** out, char* inFile, FILE** in) {
  *in = fopen(inFile, "r");
  if (*in == NULL)
    exit(error(inFile, ERROPEN));

  *out = openWrite(outFile);
}

/* void getParams()
 * Get command-line parameters.
 */
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile = NULL;
  int windowLen = 0;
  float windowAvg = 0.0f, qualAvg = 0.0f;
  int minLen = 0;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], INFILE))
        inFile = argv[++i];
      else if (!strcmp(argv[i], WINDOWLEN))
        windowLen = getInt(argv[++i]);
      else if (!strcmp(argv[i], WINDOWAVG))
        windowAvg = getFloat(argv[++i]);
      else if (!strcmp(argv[i], QUALAVG))
        qualAvg = getFloat(argv[++i]);
      else if (!strcmp(argv[i], MINLEN))
        minLen = getInt(argv[++i]);
    } else
      usage();
  }

  if (outFile == NULL || inFile == NULL)
    usage();

  FILE* out = NULL, *in = NULL;
  openFiles(outFile, &out, inFile, &in);
  readFile(in, out, windowLen, windowAvg, qualAvg, minLen);

  if (fclose(out) || fclose(in))
    exit(error("", ERRCLOSE));
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getParams(argc, argv);
  return 0;
}