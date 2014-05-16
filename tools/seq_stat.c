/*
https://github.com/noporpoise/seq_file
Isaac Turner <turner.isaac@gmail.com>
Jan 2014, Public Domain
*/

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include "seq_file.h"
#include <math.h>

static unsigned int num_of_digits(unsigned long num)
{
  unsigned int digits;
  for(digits = 1; num >= 10; digits++) num /= 10;
  return digits;
}

// result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27
// returns pointer to result
char* ulong_to_str(unsigned long num, char* result)
{
  unsigned int digits = num_of_digits(num);
  unsigned int i, num_commas = (digits-1) / 3;
  char *p = result + digits + num_commas;
  *(p--) = '\0';

  for(i = 0; i < digits; i++, num /= 10) {
    if(i > 0 && i % 3 == 0) *(p--) = ',';
    *(p--) = '0' + (num % 10);
  }

  return result;
}

// result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes
// strlen('-9,223,372,036,854,775,808') = 27
// strlen('.') = 1
// +1 for \0
char* double_to_str(double num, int decimals, char* str)
{
  if(isnan(num)) return strcpy(str, "NaN");
  else if(isinf(num)) return strcpy(str, "Inf");

  unsigned long whole_units = (unsigned long)num;
  num -= whole_units;

  char decstr[2+decimals+1];
  sprintf(decstr, "%.*lf", decimals, num);
  if(decstr[0] == '1') whole_units++;

  ulong_to_str(whole_units, str);

  if(decimals > 0)
  {
    size_t offset = strlen(str);
    strcpy(str+offset, decstr+1);
  }

  return str;
}

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

  // We've already read one read
  size_t total_len = r.seq.end, max_rlen = r.seq.end, nreads = 1;

  while(seq_read(f,&r) > 0) {
    total_len += r.seq.end;
    max_rlen = r.seq.end > max_rlen ? r.seq.end : max_rlen;
    nreads++;
  }

  double mean_rlen = (double)total_len / nreads;

  char nbasesstr[100], nreadsstr[100], maxrlenstr[100], meanrlenstr[100];
  ulong_to_str(total_len, nbasesstr);
  ulong_to_str(nreads, nreadsstr);
  ulong_to_str(max_rlen, maxrlenstr);
  double_to_str(mean_rlen, 1, meanrlenstr);

  // Trim excess zeros
  size_t len = strlen(meanrlenstr);
  while(len > 0 && meanrlenstr[len-1] == '0') len--;
  if(meanrlenstr[len-1] == '.') len--;
  meanrlenstr[len] = '\0';

  printf(" Total seq (bp):    %s\n", nbasesstr);
  printf(" Number of reads:   %s\n", nreadsstr);
  printf(" Longest read (bp): %s\n", maxrlenstr);
  printf(" Mean length  (bp): %s\n", meanrlenstr);

  printf("Done.\n");

  seq_close(f);
  seq_read_dealloc(&r);

  return EXIT_SUCCESS;
}
