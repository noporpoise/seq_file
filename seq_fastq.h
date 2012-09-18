
#ifndef SEQ_FASTQ_HEADER_SEEN
#define SEQ_FASTQ_HEADER_SEEN

#include "seq_common.h"

char seq_next_read_fastq(SeqFile *sf);

char seq_read_base_fastq(SeqFile *sf, char *c);
char seq_read_qual_fastq(SeqFile *sf, char *c);

char seq_read_all_bases_fastq(SeqFile *sf, StrBuf *sbuf);
char seq_read_all_quals_fastq(SeqFile *sf, StrBuf *sbuf);

size_t seq_file_write_name_fastq(SeqFile *sf, const char *name);
size_t seq_file_write_seq_fastq(SeqFile *sf, const char *seq, size_t str_len);
size_t seq_file_write_qual_fastq(SeqFile *sf, const char *seq);

#endif
