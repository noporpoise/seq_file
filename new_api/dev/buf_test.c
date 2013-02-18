#include "buffered_input.h"

int main(int argc, char **argv)
{
  if(argc != 2) {
    fprintf(stderr, "usage: buf_test <file>\n");
    fprintf(stderr, "  Prints lines from the file\n");
    exit(EXIT_FAILURE);
  }

  buffer_t *in = buffer_new(10);
  buffer_t *buf = buffer_new(10);
  // gzFile f = gzopen(argv[1],"r");
  FILE *f = fopen(argv[1],"r");
  if(f == NULL) exit(EXIT_FAILURE);

  int i;
  for(i = 0; freadline_buf(f,in,&buf->b,&buf->end,&buf->size) > 0; i++)
  {
    printf("line %3i: %s", i, buf->b);
    buf->begin = buf->end = 0;
  }

  fclose(f);
  buffer_free(buf);
}
