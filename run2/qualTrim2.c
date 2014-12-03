/*
  John Gaspar
  July 2014

  Quality trimming a fastq file. Prints the longest
    substring of a sequence without a bad quality window.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qualTrim.h"

// globals

// usage
void usage(void) {
  fprintf(stderr, "Usage: ./qualTrim  %s <input file>  ", INFILE);
  fprintf(stderr, "%s <output file>  ", OUTFILE);
  fprintf(stderr, "%s <window length>  ", WINDOWLEN);
  fprintf(stderr, "%s <avg. quality>\n", WINDOWAVG);
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
  else if (err == ERRINT)
    msg2 = MERRINT;
  else if (err == ERRFLOAT)
    msg2 = MERRFLOAT;
  else if (err == ERRLEN)
    msg2 = MERRLEN;

  fprintf(stderr, "Error! %s: %s\n", msg, msg2);
  return -1;
}

// memalloc
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


// int getRegion()
// Determine endpoints of good quality window.
int getRegion(char* line, int len, float qual, int st, int* last, int end) {
  int start = -1;
  int i = st;
  if (!i) {
    for (i = 1; i < len; i++) {
      float sum = 0.0f;
      for (int k = 0; k < i; k++)
//{
        sum += line[k] - 33;
//printf("%c", line[k]);}
//printf("\t%.1f\n", sum / i);
      if (sum / i < qual) {
        i = 1;
        break;
      }
    }
    start = 0;
    i = 0;
  }
  for ( ; i < end - len + 1; i++) {
    float sum = 0.0f;
    for (int j = 0; j < len; j++)
//{
      sum += line[i + j] - 33;
//printf("%c", line[i + j]);}
//printf("\t%.1f\n", sum / len);
//while (!getchar()) ;
    if (sum / len < qual) {
      if (start != -1) {
        *last = i;
        return start;
      }
    } else if (start == -1) {
      start = i + len - 1;
    }
  }

  for ( ; i < end; i++) {
    float sum = 0.0f;
    int j;
    for (j = 0; j < end - i; j++)
//{
      sum += line[i+j] - 33;
//printf("%c", line[i+j]);}
//printf("\t%.1f\n", sum / j);
//printf("i is %d, j is %d\n", i, j);
//printf("at length %d, pos %d, window is %.1f\n", i, j - i + 1, sum / i);
    if (sum / j < qual) {
      if (start != -1) {
        *last = i;
        return start;
      }
    }
  }

  *last = end;
  return start;
}

// trimQual
// Finds longest good quality window.
int trimQual(char* line, int len, float qual, int* end) {
  int maxLen = 0, maxSt = 0, maxEnd = 0;
  for (int i = 0; i < *end; ) {
//printf("trimQual %d\n", i);
    int last = 0;
    int start = getRegion(line, len, qual, i, &last, *end);

    if (start == -1)
      break;

/*printf("st: %d, end: %d\n", start, last);
for (int j = start; j < last; j++)
  printf("%c", line[j]);
printf("\n");
//while (!getchar()) ;*/

    if (last - start > maxLen) {
      maxLen = last - start;
      maxSt = start;
      maxEnd = last;
    }
    if (last > *end - len)
      break;
    i = last ? last + 1 : i + 1;
//printf("maxSt is %d, maxEnd is %d\n", maxSt, maxEnd);
  }
  *end = maxEnd;
  return maxSt;
}


// readFile
void readFile(FILE* in, FILE* out, int len, float qual) {
  char* head = (char*) memalloc(HEADER);
  char* seq = (char*) memalloc(MAX_SIZE);
  char* line = (char*) memalloc(MAX_SIZE);

  int count = 0;
  while (fgets(head, HEADER, in) != NULL) {
    if (head[0] != '@')
      continue;

    count++;
    if (fgets(seq, MAX_SIZE, in) == NULL)
      exit(error("", ERRSEQ));

    for (int i = 0; i < 2; i++)
      if (fgets(line, MAX_SIZE, in) == NULL)
        exit(error("", ERRSEQ));

    int end = strlen(line) - 1;
    if (line[end] == '\n')
      line[end] = '\0';
//printf("length %d\n%s%s\n", end, seq, line);
//while (!getchar()) ;
    if (end < len)
      exit(error(head, ERRLEN));

    int st = trimQual(line, len, qual, &end);
/*printf("st = %d, end = %d\n", st, end);
printf("the winner is:\n");
for (int j = st; j < end; j++)
  printf("%c", line[j]);
printf("\n");
while (!getchar()) ;*/

    if (st < end) {
      fprintf(out, "%s", head);
      for (int i = st; i < end; i++)
        fprintf(out, "%c", seq[i]);
      fprintf(out, "\n+\n");
      for (int i = st; i < end; i++)
        fprintf(out, "%c", line[i]);
      fprintf(out, "\n");
    }
  }
  printf("Reads analyzed: %d\n", count);

  free(seq);
  free(line);
  free(head);
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
void openFiles(char* outFile, FILE** out, char* inFile, FILE** in) {
  *in = fopen(inFile, "r");
  if (*in == NULL)
    exit(error(inFile, ERROPEN));

  *out = openWrite(outFile);
}

// getParams
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile = NULL;
  int windowLen = 0;
  float windowAvg = 0.0f;

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
    } else
      usage();
  }

  if (outFile == NULL || inFile == NULL ||
      !windowLen || !windowAvg)
    usage();

  FILE* out = NULL, *in = NULL;
  openFiles(outFile, &out, inFile, &in);
  readFile(in, out, windowLen, windowAvg);

  if (fclose(out) || fclose(in))
    exit(error("", ERRCLOSE));
}

// main
int main(int argc, char* argv[]) {
  getParams(argc, argv);
  return 0;
}
