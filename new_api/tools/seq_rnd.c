
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <stdint.h>

char *cmdstr;

static char parse_entire_uint(char *str, uint32_t *result)
{
  char *tmp_str = str;
  long tmp = strtol(str, &tmp_str, 10);

  if(tmp > UINT_MAX || tmp < 0 || tmp_str != str+strlen(str)) return 0;

  *result = (uint32_t)tmp;
  return 1;
}

static void print_usage()
{
  fprintf(stderr, "Usage: %s [len]\n", cmdstr);
  fprintf(stderr, "  Print random DNA sequence\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  cmdstr = argv[0];

  if(argc > 2) print_usage();

  srand(time(NULL));
  uint32_t len = 0, i;
  if(argc > 1 && !parse_entire_uint(argv[1], &len)) print_usage();

  char alphabet[] = "ACGT";
  for(i = 0; len == 0 || i < len; i++)
    putc(alphabet[rand()&0x3], stdout);
  putc('\n', stdout);

  return EXIT_SUCCESS;
}
