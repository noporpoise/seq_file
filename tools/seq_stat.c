/*
https://github.com/noporpoise/seq_file
Isaac Turner <turner.isaac@gmail.com>
Jan 2014, Public Domain
*/

#include "seq_file.h"

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    return -1;
  }

  char *file = argv[1];

  printf("File: %s\n", file);

  read_t r;
  seq_read_alloc(&r);
  seq_file_t *f = seq_open(file);

  if(f == NULL)
  {
    fprintf(stderr, "Cannot open file\n");
    exit(EXIT_FAILURE);
  }

  int minq = -1, maxq = -1, s, fmt;

  fmt = seq_guess_fastq_format(f, &minq, &maxq);
  s = seq_read(f,&r);

  if(s < 0) {
    fprintf(stderr, "Error occurred reading file\n");
    return EXIT_FAILURE;
  }

  if(s == 0) {
    fprintf(stderr, "Cannot get any reads from file\n");
    return EXIT_FAILURE;
  }

  const char zstr[] = " (read with zlib)";

  if(seq_is_sam(f)) printf("Format: SAM\n");
  if(seq_is_bam(f)) printf("Format: BAM\n");
  if(seq_is_fasta(f)) printf("Format: FASTA%s\n", seq_use_gzip(f) ? zstr : "");
  if(seq_is_fastq(f)) printf("Format: FASTQ%s\n", seq_use_gzip(f) ? zstr : "");
  if(seq_is_plain(f)) printf("Format: plain%s\n", seq_use_gzip(f) ? zstr : "");

  char print_qstat = (seq_is_fastq(f) || seq_is_sam(f) || seq_is_bam(f));

  if(print_qstat)
  {
    if(fmt == -1) printf("Couldn't get any quality scores\n");
    else {
      printf("Quality scores: %s, offset: %i, min: %i, max: %i, scores: [%i,%i]\n",
             FASTQ_FORMATS[fmt], FASTQ_OFFSET[fmt], FASTQ_MIN[fmt], FASTQ_MAX[fmt],
             FASTQ_MIN[fmt]-FASTQ_OFFSET[fmt], FASTQ_MAX[fmt]-FASTQ_OFFSET[fmt]);

      printf("Quality ASCII range in first 500bp: [%i,%i]\n", minq, maxq);
    }
  }

  // while(seq_read(f,&r) > 0)
  //   seq_print_fastq(&r, stdout, 0);

  printf("Done.\n");

  seq_close(f);
  seq_read_dealloc(&r);

  return EXIT_SUCCESS;
}
