#include <limits.h> // INT_MAX
#include <ctype.h> // toupper tolower
#include "seq_file.h"

SETUP_SEQ_FILE();

char parse_entire_int(char *str, int *result)
{
  char *tmp_str = str;
  long tmp = strtol(str, &tmp_str, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || tmp_str != str+strlen(str)) return 0;

  *result = (int)tmp;
  return 1;
}

void print_usage(const char *cmd, const char *err)
{
  if(err != NULL) fprintf(stderr, "%s: %s\n", cmd, err);
  else
  {
    fprintf(stderr, "Usage: %s [OPTIONS] <file1> [file2] ..\n", cmd);
    fprintf(stderr, "  -w   wrap lines by <wrap> characters\n");
    fprintf(stderr, "  -uc  convert sequence to uppercase\n");
    fprintf(stderr, "  -lc  convert sequence to lowercase\n");
    fprintf(stderr, "  -h   show this help text\n");
  }
  exit(EXIT_FAILURE);
}

#define printwrap(str,len,wrap,i,j) do { \
    for(i=0,j=0;i<len;i++,j++) { \
      if(j==wrap) { putc('\n',stdout); j = 0; } \
      putc(str[i],stdout); \
    } putc('\n',stdout); \
  } while(0)

#define fixcase(str,len,i,func) do { \
    for((i)=0; (i)<(len); (i)++) (str)[(i)] = func((str)[(i)]); \
  } while(0)

int main(int argc, char **argv)
{
  // linewrap: 0 => don't
  // change case: 0 => don't, 1 => uppercase, 2 => lowercase
  int argi, linewrap = 0, change_case = 0;

  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "-w") == 0)
    {
      if(argi == argc-1)
        print_usage(argv[0], "-w requires an argument");
      if(!parse_entire_int(argv[++argi], &linewrap) || linewrap <= 0)
        print_usage(argv[0], "invalid -w argument");
    }
    else if(strcasecmp(argv[argi], "-uc") == 0) change_case = 1;
    else if(strcasecmp(argv[argi], "-lc") == 0) change_case = 2;
    else if(strcasecmp(argv[argi], "-h") == 0) print_usage(argv[0], NULL);
    else break;
  }

  if(argi == argc) print_usage(argv[0], NULL);

  read_t *r = seq_read_alloc();
  seq_file_t *f;

  for(; argi < argc; argi++)
  {
    if((f = seq_open(argv[argi])) == NULL)
    {
      fprintf(stderr, "Cannot open file %s\n", argv[argi]);
      exit(EXIT_FAILURE);
    }

    int i,j;
    while(seq_read(f,r) > 0)
    {
      if(change_case == 1) fixcase(r->seq.b,r->seq.end,i,toupper);
      if(change_case == 2) fixcase(r->seq.b,r->seq.end,i,tolower);

      #ifdef FASTQ
        if(r->qual.end == 0)
        {
          buffer_ensure_capacity(&r->qual, r->seq.end);
          r->qual.end = r->seq.end;
          memset(r->qual.b,'.',r->qual.end);
          r->qual.b[r->qual.end] = 0;
        }
        if(linewrap == 0)
        {
          printf("@%s\n%s\n+\n%s\n", r->name.b, r->seq.b, r->qual.b);
        }
        else
        {
          printf("@%s\n", r->name.b);
          printwrap(r->seq.b, r->seq.end, linewrap, i, j);
          printf("+\n");
          printwrap(r->qual.b, r->qual.end, linewrap, i, j);
        }
      #elif defined(FASTA)
        if(linewrap == 0) printf(">%s\n%s\n", r->name.b, r->seq.b);
        else
        {
          printf(">%s\n", r->name.b);
          printwrap(r->seq.b, r->seq.end, linewrap, i, j);
        }
      #else
        // Plain (sequence only, one per line, no wrap option)
        printf("%s\n", r->seq.b);
      #endif
    }

    seq_close(f);
  }
  seq_read_destroy(r);
}
