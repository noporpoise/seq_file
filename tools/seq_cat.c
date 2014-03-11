/*
https://github.com/noporpoise/seq_file
Isaac Turner <turner.isaac@gmail.com>
Jan 2014, Public Domain
*/

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <limits.h> // INT_MAX
#include <ctype.h> // toupper tolower
#include "seq_file.h"

#if !defined(FASTA) && !defined(FASTQ) && !defined(PLAIN) && !defined(REVCMP)
#error You must define one of: FASTA FASTQ PLAIN REVCMP
#endif

#define UPPERCASE 1
#define LOWERCASE 2

char *cmdstr;

char parse_entire_uint(char *str, uint32_t *result)
{
  char *tmp_str = str;
  long tmp = strtol(str, &tmp_str, 10);

  if(tmp > UINT_MAX || tmp < 0 || tmp_str != str+strlen(str)) return 0;

  *result = (uint32_t)tmp;
  return 1;
}

void print_usage(const char *err) __attribute__((noreturn));

void print_usage(const char *err)
{
  if(err != NULL) fprintf(stderr, "%s: %s\n", cmdstr, err);
  else
  {
    fprintf(stderr, "Usage: %s [OPTIONS] [file1] [file2] ..\n", cmdstr);

    #if defined(FASTA)
    fprintf(stderr, "  Print files in FASTA format\n");
    #elif defined(FASTQ)
    fprintf(stderr, "  Print files in FASTQ format\n");
    #elif defined(PLAIN)
    fprintf(stderr, "  Print files in 'plain' format -- one sequence per line\n");
    #elif defined(REVCMP)
    fprintf(stderr, "  Print files with reads reverse complemented\n");
    #endif

    fprintf(stderr, "\n  OPTIONS:\n");
    #if defined(FASTA)
    fprintf(stderr, "   -w <n>  wrap lines by <n> characters [default: 80]\n");
    #elif defined(FASTQ) || defined(REVCMP)
    fprintf(stderr, "   -w <n>  wrap lines by <n> characters [default: 0 (off)]\n");
    #endif
    fprintf(stderr, "   -uc     convert sequence to uppercase\n");
    fprintf(stderr, "   -lc     convert sequence to lowercase\n");
    fprintf(stderr, "   -h      show this help text\n");
  }
  exit(EXIT_FAILURE);
}

static void seq_cat(const char *file, read_t *r,
                    uint32_t change_case, uint32_t linewrap)
{
  seq_file_t *f;

  if((f = seq_open(file)) == NULL) {
    fprintf(stderr, "Cannot open file %s\n", file);
    exit(EXIT_FAILURE);
  }

  while(seq_read(f,r) > 0)
  {
    if(change_case == UPPERCASE) seq_read_to_uppercase(r);
    if(change_case == LOWERCASE) seq_read_to_lowercase(r);

    #if defined(REVCMP)
      seq_read_reverse_complement(r);
      if(seq_is_fastq(f) || seq_is_sam(f) || seq_is_bam(f))
        seq_print_fastq(r, stdout, linewrap);
      else if(seq_is_plain(f)) puts(r->seq.b);
      else seq_print_fasta(r, stdout, linewrap);
    #elif defined(FASTQ)
      seq_print_fastq(r, stdout, linewrap);
    #elif defined(FASTA)
      seq_print_fasta(r, stdout, linewrap);
    #elif defined(PLAIN)
      // sequence only, one per line, no wrap option - don't print if empty
      (void)linewrap;
      if(r->seq.end > 0) printf("%s\n", r->seq.b);
    #endif
  }

  seq_close(f);
}

int main(int argc, char **argv)
{
  cmdstr = argv[0];

  // linewrap: 0 => don't
  // change case: 0 [don't], 1 [UPPERCASE], 2 [LOWERCASE]
  uint32_t change_case = 0, linewrap = 0;
  int argi;

  #if defined(FASTA)
  linewrap = 80;
  #endif

  for(argi = 1; argi < argc; argi++)
  {
    #if defined(FASTA) || defined(FASTQ) || defined(REVCMP)
    if(strcasecmp(argv[argi], "-w") == 0)
    {
      if(argi == argc-1) print_usage("-w <n> requires an argument");
      if(!parse_entire_uint(argv[++argi], &linewrap))
        print_usage("invalid -w argument");
    } else
    #endif
    if(strcasecmp(argv[argi], "-uc") == 0) change_case = UPPERCASE;
    else if(strcasecmp(argv[argi], "-lc") == 0) change_case = LOWERCASE;
    else if(argv[argi][0] == '-' && argv[argi][1] != '\0') print_usage(NULL);
    else break;
  }

  read_t r;
  seq_read_alloc(&r);

  if(argi == argc)
    seq_cat("-", &r, change_case, linewrap);
  else {
    for(; argi < argc; argi++)
      seq_cat(argv[argi], &r, change_case, linewrap);
  }

  seq_read_dealloc(&r);

  return EXIT_SUCCESS;
}
