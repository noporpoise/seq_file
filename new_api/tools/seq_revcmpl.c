#include "seq_file.h"

#define CMDSTR "seqrevcmpl"

void print_usage(const char *err)
{
  if(err != NULL) fprintf(stderr, "%s: %s\n", CMDSTR, err);
  else
  {
    fprintf(stderr, "Usage: %s [OPTIONS] <file1> [file2] ..\n", CMDSTR);
    fprintf(stderr, "  Reverse complement sequences\n");

    fprintf(stderr, "\n  OPTIONS:\n");
    fprintf(stderr, "   -w <n>  wrap lines by <n> characters (FASTA/FASTQ only)\n");
    fprintf(stderr, "   -uc     convert sequence to uppercase\n");
    fprintf(stderr, "   -lc     convert sequence to lowercase\n");
    fprintf(stderr, "   -h      show this help text\n");
  }
  exit(EXIT_FAILURE);
}

char parse_entire_uint(char *str, size_t *result)
{
  char *tmp_str = str;
  long tmp = strtol(str, &tmp_str, 10);

  if(tmp > UINT_MAX || tmp < 0 || tmp_str != str+strlen(str)) return 0;

  *result = (int)tmp;
  return 1;
}

#define fixcase(str,len,i,func) do { \
    for((i)=0; (i)<(len); (i)++) (str)[(i)] = func((str)[(i)]); \
  } while(0)

int main(int argc, char **argv)
{
  // linewrap: 0 => don't
  // change case: 0 => don't, 1 => uppercase, 2 => lowercase
  int argi, change_case = 0;
  size_t linewrap = 80;

  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "-w") == 0)
    {
      if(argi == argc-1) print_usage("-w <n> requires an argument");
      if(!parse_entire_uint(argv[++argi], &linewrap) || linewrap <= 0)
        print_usage("invalid -w argument");
    }
    else if(strcasecmp(argv[argi], "-uc") == 0) change_case = 1;
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
      fprintf(stderr, "%s: Cannot open file %s\n", CMDSTR, argv[argi]);
      exit(EXIT_FAILURE);
    }

    size_t i;
    while(seq_read(f,r) > 0)
    {
      if(change_case == 1) fixcase(r->seq.b,r->seq.end,i,toupper);
      if(change_case == 2) fixcase(r->seq.b,r->seq.end,i,tolower);

      seq_read_reverse_complement(r);

      if(seq_is_fasta(f)) seq_print_fasta(r, stdout, linewrap);
      else if(seq_is_fastq(f)) seq_print_fastq(r, stdout, linewrap);
      else printf("%s\n", r->seq.b);
    }

    seq_close(f);
  }
  seq_read_destroy(r);

  return EXIT_SUCCESS;
}
