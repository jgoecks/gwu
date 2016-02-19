/*
  John Gaspar
  June 2014

  Removing duplicate (and subset) sequences
    from an input fasta or fastq file.
  Version 2: takes multiple input files.
  Version 3: keeping 1st read as rep.
  Version 4: minimizing memory -- not keeping track of dups
  Version 5: continuing to min. memory, execution time
  Version 6: again adjusting the Node** structure
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "removeDups6.h"

// globals
static char* line;
static char* hline;
static Node* root;
//int ncount;

void usage(void) {
  fprintf(stderr, "Usage: ./removeDups  %s <input file>  ", INFILE);
  fprintf(stderr, "%s <output file>\n", OUTFILE);
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

  fprintf(stderr, "Error! %s: %s\n", msg, msg2);
  return -1;
}

// freeNode
void freeNode(Node* n) {
  free(n->seq);
  if (n->head != NULL)
    free(n->head);
  for (int i = 0; i < CHILD; i++)
    if (n->child[i] != NULL)
      freeNode(n->child[i]);
    else
      break;
  free(n);
}

// freeMemory
void freeMemory(void) {
  for (int i = 0; i < CHILD; i++)
    if (root->child[i] != NULL)
      freeNode(root->child[i]);
    else
      break;
  free(root);
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

// printTrie()
void printTrie(FILE* out, FILE* dup, Node* n, int pos,
    int* leaf, int* rem) {
  if (n == NULL)
    return;

  for (int i = n->st; i < n->end; i++)
    line[pos++] = n->seq[i];
  if (n->child[0] == NULL) {
    line[pos] = '\0';
    fprintf(out, ">%s\n%s\n", n->head, line);
    (*leaf)++;
  }
  for (int i = 0; i < CHILD; i++)
    if (n->child[i] != NULL)
      printTrie(out, dup, n->child[i], pos, leaf, rem);
    else
      break;
}


// dumpTrie -- for debugging
/*void dumpTrie(Node* first, int level) {
  for (Node* n = first; n != NULL; n = n->next) {
    for (int i = 0; i < level; i++)
      printf(" ");
    for (int i = n->st; i < n->end; i++)
      printf("%c", n->seq[i]);
    //for (Read* r = n->first; r != NULL; r = r->next)
      //printf(" ->%s", r->header);
    printf("\n");
    dumpTrie(n->child, level + n->end - n->st);
  }
}*/

// makeRead
/*Read* makeRead(void) {
  Read* r = (Read*) memalloc(sizeof(Read));
  r->header = (char*) memalloc(HEADER);
  int i;
  for (i = 1; line[i] != '\0' && line[i] != '\n' && i < HEADER; i++)
    r->header[i - 1] = line[i];
  r->header[i - 1] = '\0';
  return r;
}*/

// char* copySeq()
char* copySeq(char* src, int pos, int size) {
  char* seq = (char*) memalloc(size);
  int i;
  for (i = 0; i < size; i++)
    seq[i] = src[pos + i];
  return seq;
}

// Node* makeNode()
Node* makeNode(Node* par, int ch, int pos, int size) {
  /*int i;
  for (i = 0; i < CHILD; i++)
    if (par->child[i] == NULL)
      break;*/
  Node* n = (Node*) memalloc(sizeof(Node));
  n->head = NULL;
  par->child[ch] = n;
  if (++ch < CHILD)
    par->child[ch] = NULL;
  n->child[0] = NULL;

  n->seq = copySeq(line, pos, size);
  n->st = 0;
  n->end = size;
//ncount++;
  return n;
}

/* void splitNode()
 * Splits a node by creating a new child.
 */
Node* splitNode(Node* par, Node* gp, int ch, int size) {
  Node* n = (Node*) memalloc(sizeof(Node));
//ncount++;
  gp->child[ch] = n;
  n->child[0] = par;
  n->child[1] = NULL;
  n->head = NULL;
/*printf("splitNode\n");
for (int i = par->st; i < par->end; i++)
  printf("%c", par->seq[i]);
printf("\n");*/

  // parent gets new seq
  if (size > (par->end - par->st) / 2.0) {
/*printf("splitNode\n");
for (int i = par->st; i < par->end; i++)
  printf("%c", par->seq[i]);
printf("\n");
printf("parent gets new\n");*/
    n->seq = par->seq;
    n->st = par->st;
    n->end = par->st + size;

    par->st = 0;
    par->end = par->end - par->st - size;
    par->seq = copySeq(par->seq, par->st + size, par->end);

/*printf("results in:\nchild ");
for (int i = n->st; i < n->end; i++)
  printf("%c", n->seq[i]);
printf("\nparent ");
for (int i = par->st; i < par->end; i++)
  printf("%c", par->seq[i]);
printf("\n");
while (!getchar()) ;*/

  } else {
//printf("child gets new\n");
    // child gets new seq
    n->st = 0;
    n->end = size;
    n->seq = copySeq(par->seq, par->st, size);
    par->st = par->st + size;
  }

/*printf("results in:\nchild ");
for (int i = n->st; i < n->end; i++)
  printf("%c", n->seq[i]);
printf("\nparent ");
for (int i = par->st; i < par->end; i++)
  printf("%c", par->seq[i]);
printf("\n");
while (!getchar()) ;*/
  return n;
}


// declaration for recursive call
Node* checkNode(Node*, int, int);

// moreSeq -- check remaining sequence of node
Node* moreSeq(Node* n, Node* par, int ch, int pos, int len) {
  int i;
  for (i = 1; i < n->end - n->st; i++) {

    // read is short
    if (pos + i == len) {
      return splitNode(n, par, ch, i);
      //return n;
    }

    // mismatch at position i
    if (line[pos + i] != n->seq[n->st + i]) {
      Node* m = splitNode(n, par, ch, i);
      return makeNode(m, 1, pos + i, len - pos - i);
    }
  }

  // recurse on child node
  return checkNode(n, pos + i, len);
}

// checkNode
Node* checkNode(Node* n, int pos, int len) {

  // base case
  if (pos == len)
    return n;

  int i;
  for (i = 0; i < CHILD; i++) {
    Node* m = n->child[i];
    if (m == NULL)
      break;
    if (line[pos] == m->seq[m->st])
      return moreSeq(m, n, i, pos, len);
  }
/*if (i == CHILD) {
  printf("that's a problem\n");
  for (int j = 0; j < CHILD; j++) {
    Node* m = n->child[j];
    printf("%d: %c\n", j, m->seq[m->st]);
  }
  printf("doesn't match %c\n", line[pos]);
  while (!getchar()) ;
}
/*if (i == CHILD - 1) {
  printf("borderline\n");
  for (int j = 0; j < CHILD; j++) {
    Node* m = n->child[j];
    if (m == NULL)
      break;
    printf("%d: %c\n", j, m->seq[m->st]);
  }
  printf("doesn't match %c\n", line[pos]);
  for (int j = 0; j < pos + 1; j++)
    printf("%c", line[j]);
  printf("\n");
  while (!getchar()) ;
}*/

  return makeNode(n, i, pos, len - pos);
}

// readFile
int readFile(FILE* in) {
  int aorq = fastaOrQ(in);
  rewind(in);
  int count = 0;
  while (fgets(hline, MAX_SIZE, in) != NULL) {
    if (hline[0] == '#')
      continue;

    count++;
    if (fgets(line, MAX_SIZE, in) == NULL)
      exit(error("", ERRSEQ));
    int len = strlen(line) - 1;
    if (line[len] == '\n')
      line[len] = '\0';

    // add read to trie
    Node* n = checkNode(root, 0, len);
    if (n->child[0] == NULL && n->head == NULL) {
      n->head = (char*) memalloc(HEADER);
      int i;
      for (i = 1; hline[i] != '\0' && hline[i] != '\n' && i < HEADER; i++)
        n->head[i - 1] = hline[i];
      n->head[i - 1] = '\0';
    }

    // skip next 2 lines if fastq
    if (aorq)
      for (int i = 0; i < 2; i++)
        if (fgets(line, MAX_SIZE, in) == NULL)
          exit(error("", ERRSEQ));
  }
  printf("  reads: %d\n", count);
  return count;
}

// parseFiles -- loop through input files
void parseFiles(char* inFile) {
  int total = 0;
  char* file = strtok(inFile, CSV);
  while (file != NULL) {
    printf("Reading %s\n", file);
    FILE* in = fopen(file, "r");
    if (in == NULL)
      exit(error(file, ERROPEN));

    total += readFile(in);
    if (fclose(in))
      exit(error("", ERRCLOSE));

    file = strtok(NULL, CSV);
  }
  printf("Total reads:  %10d\n", total);
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
void openFiles(char* outFile, FILE** out, char* dupFile, FILE** dup) {
  *out = openWrite(outFile);
  if (dupFile != NULL)
    *dup = openWrite(dupFile);
}

// getParams
void getParams(int argc, char** argv) {

  char* outFile = NULL, *inFile = NULL, *dupFile = NULL;

  // parse argv
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], HELP))
      usage();
    else if (i < argc - 1) {
      if (!strcmp(argv[i], OUTFILE))
        outFile = argv[++i];
      else if (!strcmp(argv[i], INFILE))
        inFile = argv[++i];
      else if (!strcmp(argv[i], DUPFILE))
        dupFile = argv[++i];
    } else
      usage();
  }

  if (outFile == NULL || inFile == NULL)
    usage();

  FILE* out = NULL, *dup = NULL;
  openFiles(outFile, &out, dupFile, &dup);

  parseFiles(inFile);

  int leaf = 0, rem = 0;
  for (int i = 0; i < CHILD; i++)
    printTrie(out, dup, root->child[i], 0, &leaf, &rem);
  printf("Leaf reads:   %10d\n", leaf);
  if (dup != NULL)
    printf("Removed reads:%10d\n", rem);

  if (fclose(out) || (dup != NULL && fclose(dup)))
    exit(error("", ERRCLOSE));
}

// main
int main(int argc, char* argv[]) {
  root = (Node*) memalloc(sizeof(Node));
  root->child[0] = NULL;
  line = (char*) memalloc(MAX_SIZE);
  hline = (char*) memalloc(MAX_SIZE);
  getParams(argc, argv);
  freeMemory();
  return 0;
}
