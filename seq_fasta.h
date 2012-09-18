
#ifndef SEQ_FASTA_HEADER_SEEN
#define SEQ_FASTA_HEADER_SEEN

#include "seq_common.h"

char seq_next_read_fasta(SeqFile *sf);
char seq_read_base_fasta(SeqFile *sf, char *c);
char seq_read_all_bases_fasta(SeqFile *sf, StrBuf *sbuf);

size_t seq_file_write_name_fasta(SeqFile *sf, const char *name);
size_t seq_file_write_seq_fasta(SeqFile *sf, const char *seq, size_t str_len);

#endif
