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
//  if (n == NULL)
  //  return;
  free(n->seq);
/*  Read* r = n->first;
  Read* tempr;
  while (r != NULL) {
    free(r->header);
    tempr = r;
    r = r->next;
    free(tempr);
  }*/
  if (n->head != NULL)
    free(n->head);
  for (int i = 0; i < CHILD && n->child[i] != NULL; i++)
    freeNode(n->child[i]);
//  free(n->child);
  free(n);
}

// freeMemory
void freeMemory(void) {
  for (int i = 0; i < CHILD && root->child[i] != NULL; i++)
    freeNode(root->child[i]);
//  free(root->child);
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

//  int oldPos = pos;
  for (int i = n->st; i < n->end; i++)
    line[pos++] = n->seq[i];
//  if (n->first != NULL) {
    //Read* r = n->first;
    if (n->child[0] == NULL) {
      line[pos] = '\0';
      fprintf(out, ">%s\n%s\n", n->head, line);
      (*leaf)++;
      //r = r->next;
    }
    /*if (dup != NULL)
      for ( ; r != NULL; r = r->next) {
        fprintf(dup, ">%s\n%s\n", r->header, line);
        (*rem)++;
      }*/
//  }
  for (int i = 0; i < CHILD && n->child[i] != NULL; i++)
    printTrie(out, dup, n->child[i], pos, leaf, rem);
  //printTrie(out, dup, n->child, pos, leaf, rem);
  //printTrie(out, dup, n->next, oldPos, leaf, rem);
}


// dumpTrie -- for debugging
void dumpTrie(Node* n, int level) {
  if (n->seq != NULL) {
    for (int i = 0; i < level; i++)
      printf(" ");
    for (int i = n->st; i < n->end; i++)
      printf("%c", n->seq[i]);
    //for (Read* r = n->first; r != NULL; r = r->next)
      //printf(" ->%s", r->header);
    printf("\n");
  }
  for (int i = 0; i < CHILD; i++)
    if (n->child[i] != NULL)
      dumpTrie(n->child[i], level + n->end - n->st);
    else
      break;
}

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
printf("makeNode added %d\n", size);
  return n;
}

/* void splitNode()
 * Splits a node by creating a new child.
 */
void splitNode(Node* par, int size) {
  Node* n = (Node*) memalloc(sizeof(Node));
//ncount++;
  //n->next = NULL;
  //n->child = (Node**) memalloc(CHILD * sizeof(Node*));
  int i;
  for (i = 0; i < CHILD && par->child[i] != NULL; i++)
    n->child[i] = par->child[i];
  if (i < CHILD)
    n->child[i] = NULL;

//  n->child = par->child;
//printf("splitnode: n->child is %08x\n%08x\n", (int)n->child[0], (int)n->child[1]);
  n->head = par->head;

  // parent gets new seq
  if (size < (par->end - par->st) / 2.0) {
    n->seq = par->seq;
    n->st = par->st + size;
    n->end = par->end;

    par->seq = copySeq(par->seq, par->st, size);
    par->st = 0;
    par->end = size;

  } else {
    // child gets new seq
    n->st = 0;
    n->end = par->end - par->st - size;
    n->seq = copySeq(par->seq, par->st + size, n->end);
    par->end = par->st + size;
  }

  par->head = NULL;
  par->child[0] = n;
  par->child[1] = NULL;
}


// declaration for recursive call
Node* checkNode(Node*, int, int);

// moreSeq -- check remaining sequence of node
Node* moreSeq(Node* n, int pos, int len) {
  int i;
  for (i = 1; i < n->end - n->st; i++) {

    // read is short
    if (pos + i == len) {
      splitNode(n, i);
      return n;
    }

    // mismatch at position i
    if (line[pos + i] != n->seq[n->st + i]) {
      splitNode(n, i);
      return makeNode(n, 1, pos + i, len - pos - i);
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
      return moreSeq(m, pos, len);
  }

  //return makeNode(n, i, pos, len);
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
dumpTrie(root, 0);
while (!getchar()) ;
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
