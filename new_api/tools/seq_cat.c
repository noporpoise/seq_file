#include <limits.h> // INT_MAX
#include <ctype.h> // toupper tolower
#include "seq_file.h"

#if defined(FASTA)
#define CMDSTR "facat"
#elif defined(FASTQ)
#define CMDSTR "fqcat"
#elif defined(PLAIN)
#define CMDSTR "seqcat"
#else
#error You must define one of: FASTA FASTQ PLAIN
#endif

char parse_entire_uint(char *str, size_t *result)
{
  char *tmp_str = str;
  long tmp = strtol(str, &tmp_str, 10);

  if(tmp > UINT_MAX || tmp < 0 || tmp_str != str+strlen(str)) return 0;

  *result = (int)tmp;
  return 1;
}

void print_usage(const char *err)
{
  if(err != NULL) fprintf(stderr, "%s: %s\n", CMDSTR, err);
  else
  {
    fprintf(stderr, "Usage: %s [OPTIONS] <file1> [file2] ..\n", CMDSTR);

    #if defined(FASTA)
    fprintf(stderr, "  Print files in FASTA format\n");
    #elif defined(FASTQ)
    fprintf(stderr, "  Print files in FASTQ format\n");
    #else
    fprintf(stderr, "  Print files in 'plain' format -- one sequence per line\n");
    #endif

    fprintf(stderr, "\n  OPTIONS:\n");
    #if defined(FASTA)
    fprintf(stderr, "   -w <n>  wrap lines by <n> characters [default: 80]\n");
    #elif defined(FASTQ)
    fprintf(stderr, "   -w <n>  wrap lines by <n> characters [default: 0 (off)]\n");
    #endif
    fprintf(stderr, "   -uc     convert sequence to uppercase\n");
    fprintf(stderr, "   -lc     convert sequence to lowercase\n");
    fprintf(stderr, "   -h      show this help text\n");
  }
  exit(EXIT_FAILURE);
}

#define fixcase(str,len,i,func) do { \
    for((i)=0; (i)<(len); (i)++) (str)[(i)] = func((str)[(i)]); \
  } while(0)

int main(int argc, char **argv)
{
  // linewrap: 0 => don't
  // change case: 0 => don't, 1 => uppercase, 2 => lowercase
  int argi, change_case = 0;

  #if defined(FASTA)
  size_t linewrap = 80;
  #elif defined(FASTQ)
  size_t linewrap = 0;
  #endif

  for(argi = 1; argi < argc; argi++)
  {
    #if defined(FASTA) || defined(FASTQ)
    if(strcasecmp(argv[argi], "-w") == 0)
    {
      if(argi == argc-1) print_usage("-w <n> requires an argument");
      if(!parse_entire_uint(argv[++argi], &linewrap) || linewrap <= 0)
        print_usage("invalid -w argument");
    } else
    #endif
    if(strcasecmp(argv[argi], "-uc") == 0) change_case = 1;
    else if(strcasecmp(argv[argi], "-lc") == 0) change_case = 2;
    else if(strcasecmp(argv[argi], "-h") == 0) print_usage(NULL);
    else break;
  }

  if(argi == argc) print_usage(NULL);

  read_t *r = seq_read_alloc();
  seq_file_t *f;

  for(; argi < argc; argi++)
  {
    if((f = seq_open(argv[argi])) == NULL)
    {
      fprintf(stderr, "Cannot open file %s\n", argv[argi]);
      exit(EXIT_FAILURE);
    }

    size_t i;
    while(seq_read(f,r) > 0)
    {
      if(change_case == 1) fixcase(r->seq.b,r->seq.end,i,toupper);
      if(change_case == 2) fixcase(r->seq.b,r->seq.end,i,tolower);

      #if defined(FASTQ)
        seq_print_fastq(r, stdout, linewrap);
      #elif defined(FASTA)
        seq_print_fasta(r, stdout, linewrap);
      #elif defined(PLAIN)
        printf("%s\n", r->seq.b); // sequence only, one per line, no wrap option
      #endif
    }

    seq_close(f);
  }
  seq_read_destroy(r);

  return EXIT_SUCCESS;
}
