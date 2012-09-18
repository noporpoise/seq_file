
#ifndef SEQ_SAM_HEADER_SEEN
#define SEQ_SAM_HEADER_SEEN

#include "seq_common.h"

char seq_next_read_sam(SeqFile *sf);

char seq_read_base_sam(SeqFile *sf, char *c);
char seq_read_qual_sam(SeqFile *sf, char *c);

char seq_read_all_bases_sam(SeqFile *sf, StrBuf *sbuf);
char seq_read_all_quals_sam(SeqFile *sf, StrBuf *sbuf);

// Not implemented yet
size_t seq_file_write_name_sam(SeqFile *sf, const char *name);
size_t seq_file_write_seq_sam(SeqFile *sf, const char *seq, size_t str_len);
size_t seq_file_write_qual_sam(SeqFile *sf, const char *seq, size_t str_len);

#endif
