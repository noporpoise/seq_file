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
#include <stdarg.h>
#include <stdbool.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h> // INT_MAX
#include <ctype.h> // toupper tolower
#include "seq_file.h"

#define OPS_UPPERCASE   1
#define OPS_LOWERCASE   2
#define OPS_REVERSE     4
#define OPS_COMPLEMENT  8
#define OPS_MASK       16

#define FORMAT_DEFAULT 0
#define FORMAT_FASTA   1
#define FORMAT_FASTQ   2
#define FORMAT_PLAIN   4

static struct option longopts[] =
{
// General options
  {"help",       no_argument,       NULL, 'h'},
  {"fasta",      no_argument,       NULL, 'f'},
  {"fastq",      no_argument,       NULL, 'q'},
  {"plain",      no_argument,       NULL, 'p'},
  {"wrap",       required_argument, NULL, 'w'},
  {"uppercase",  no_argument,       NULL, 'u'},
  {"lowercase",  no_argument,       NULL, 'l'},
  {"revcmp",     no_argument,       NULL, 'r'},
  {"reverse",    no_argument,       NULL, 'R'},
  {"complement", no_argument,       NULL, 'C'},
  {"interleave", no_argument,       NULL, 'i'},
  {"mask",       no_argument,       NULL, 'm'},
  {NULL, 0, NULL, 0}
};

const char shortopts[] = "hfqpw:ulrRCim";

const char *cmdstr;

char parse_entire_size(const char *str, size_t *result)
{
  char *strtol_last_char_ptr = NULL;
  if(*str < '0' || *str > '9') return 0;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);
  if(tmp > SIZE_MAX) return 0;
  if(strtol_last_char_ptr == NULL || *strtol_last_char_ptr != '\0') return 0;
  *result = (size_t)tmp;
  return 1;
}

void print_usage(const char *err, ...)
__attribute__((noreturn))
__attribute__((format(printf, 1, 2)));


void print_usage(const char *err, ...)
{
  if(err != NULL) {
    va_list argptr;
    fprintf(stderr, "%s Error: ", cmdstr);
    va_start(argptr, err);
    vfprintf(stderr, err, argptr);
    va_end(argptr);
    fputc('\n', stderr);
  }

  fprintf(stderr, "Usage: %s [OPTIONS] <file1> [file2] ..\n", cmdstr);
  fprintf(stderr, "  Read and manipulate dna sequence.\n"
"\n"
"  -h,--help        show this help text\n"
"  -f,--fasta       print in FASTA format\n"
"  -q,--fastq       print in FASTQ format\n"
"  -p,--plain       print in plain format\n"
"  -w,--wrap <n>    wrap lines by <n> characters [default: 0 (off)]\n"
"  -u,--uppercase   convert sequence to uppercase\n"
"  -l,--lowercase   convert sequence to lowercase\n"
"  -r,--revcmp      reverse complement sequence [i.e. -R and -C]\n"
"  -R,--reverse     reverse sequence\n"
"  -C,--complement  complement sequence\n"
"  -i,--interleave  interleave input files\n"
"  -m,--mask        mask lowercase bases\n");

  exit(EXIT_FAILURE);
}

static void read_print(seq_file_t *sf, read_t *r,
                       uint8_t fmt, uint8_t ops, size_t linewrap)
{
  size_t i;
  if(ops & OPS_UPPERCASE)  seq_read_to_uppercase(r);
  if(ops & OPS_LOWERCASE)  seq_read_to_lowercase(r);
  if((ops & OPS_REVERSE) && (ops & OPS_COMPLEMENT)) seq_read_reverse_complement(r);
  else if(ops & OPS_REVERSE)    seq_read_reverse(r);
  else if(ops & OPS_COMPLEMENT) seq_read_complement(r);

  if(ops & OPS_MASK) {
    for(i = 0; i < r->seq.end; i++)
      if(islower(r->seq.b[i]))
        r->seq.b[i] = 'N';
  }

  if(fmt == 0) {
    if(seq_is_plain(sf)) fmt = FORMAT_PLAIN;
    else if(seq_is_fasta(sf)) fmt = FORMAT_FASTA;
    else fmt = FORMAT_FASTQ;
  }

  switch(fmt) {
    case FORMAT_FASTA: seq_print_fasta(r, stdout, linewrap); break;
    case FORMAT_FASTQ: seq_print_fastq(r, stdout, linewrap); break;
    case FORMAT_PLAIN: fputs(r->seq.b, stdout); fputc('\n', stdout); break;
    default: fprintf(stderr, "Got value: %i\n", (int)fmt); abort();
  }
}

int main(int argc, char **argv)
{
  cmdstr = argv[0];

  bool interleave = false;
  uint8_t ops = 0, fmt = FORMAT_DEFAULT;
  size_t linewrap = 0;

  if(argc == 1) print_usage(NULL);

  // Arg parsing
  int c;
  opterr = 0; // silence getopt error messages

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': print_usage(NULL); break;
      case 'f': fmt |= FORMAT_FASTA; break;
      case 'q': fmt |= FORMAT_FASTQ; break;
      case 'p': fmt |= FORMAT_PLAIN; break;
      case 'w':
        if(!parse_entire_size(optarg, &linewrap))
          print_usage("Bad -w argument: %s\n", optarg);
        break;
      case 'u': ops |= OPS_UPPERCASE; break;
      case 'l': ops |= OPS_LOWERCASE; break;
      case 'r': ops |= OPS_REVERSE | OPS_COMPLEMENT; break;
      case 'R': ops |= OPS_REVERSE; break;
      case 'C': ops |= OPS_COMPLEMENT; break;
      case 'm': ops |= OPS_MASK; break;
      case 'i': interleave = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        print_usage("Bad option: %s\n", argv[optind-1]);
      default: abort();
    }
  }

  if(optind >= argc) print_usage("Please specify at least one input file\n");

  if(!!(fmt&FORMAT_FASTA) + !!(fmt&FORMAT_FASTQ) + !!(fmt&FORMAT_PLAIN) > 1)
    print_usage("Please specify only one output format (-f,-q,-p)\n");

  size_t num_inputs = argc - optind;
  char **input_paths = argv + optind;

  read_t r;
  seq_read_alloc(&r);

  seq_file_t *inputs[num_inputs];
  size_t i;

  for(i = 0; i < num_inputs; i++) {
    if((inputs[i] = seq_open(input_paths[i])) == NULL)
      print_usage("Couldn't read file: %s\n", input_paths[i]);
  }

  if(interleave) {
    // read one entry from each file
    size_t waiting_files = num_inputs;
    while(waiting_files) {
      for(i = 0; i < num_inputs; i++) {
        if(inputs[i] != NULL) {
          if(seq_read(inputs[i],&r) > 0) read_print(inputs[i], &r, fmt, ops, linewrap);
          else { seq_close(inputs[i]); inputs[i] = NULL; waiting_files--; }
        }
      }
    }
  }
  else {
    for(i = 0; i < num_inputs; i++) {
      while(seq_read(inputs[i],&r) > 0) read_print(inputs[i], &r, fmt, ops, linewrap);
      seq_close(inputs[i]);
    }
  }

  seq_read_dealloc(&r);

  return EXIT_SUCCESS;
}
