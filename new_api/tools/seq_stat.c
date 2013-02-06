
#include "seq_file.h"

SETUP_SEQ_FILE();

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    return -1;
  }

  char *file = argv[1];

  printf("File: %s\n", file);

  read_t *r = seq_read_alloc();
  seq_file_t *f = seq_open(file);

  if(f == NULL)
  {
    fprintf(stderr, "Cannot open file\n");
    exit(EXIT_FAILURE);
  }

  int s = seq_read(f,r);

  if(s < 0) {
    printf("Error occurred reading file\n");
    return -1;
  }

  if(s == 0) {
    printf("Cannot get any reads from file\n");
    return -1;
  }

  const char *zstr = " (read with zlib)";

  if(seq_is_sam(f)) printf("Format: SAM\n");
  if(seq_is_bam(f)) printf("Format: BAM\n");
  if(seq_is_fasta(f)) printf("Format: FASTA%s\n", seq_use_gzip(f) ? zstr : "");
  if(seq_is_fastq(f)) printf("Format: FASTQ%s\n", seq_use_gzip(f) ? zstr : "");
  if(seq_is_plain(f)) printf("Format: plain%s\n", seq_use_gzip(f) ? zstr : "");

  char print_qstat = (seq_is_fastq(f) || seq_is_sam(f) || seq_is_bam(f));

  seq_close(f);
  seq_read_destroy(r);

  if(print_qstat)
  {
    int fmt = seq_guess_fastq_format(file);
    printf("Quality scores: %s, offset: %i, min: %i, max: %i, scores: [%i,%i]\n",
           FASTQ_FORMATS[fmt], FASTQ_OFFSET[fmt], FASTQ_MIN[fmt], FASTQ_MAX[fmt],
           FASTQ_MIN[fmt]-FASTQ_OFFSET[fmt], FASTQ_MAX[fmt]-FASTQ_OFFSET[fmt]);
  
    int minq = -1, maxq = -1;
    seq_get_qual_limits(file, 500, &minq, &maxq);
    printf("Quality ASCII range in first 500bp: [%i,%i]\n", minq, maxq);
  }

  return EXIT_SUCCESS;
}
