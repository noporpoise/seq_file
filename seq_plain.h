
#ifndef SEQ_PLAIN_HEADER_SEEN
#define SEQ_PLAIN_HEADER_SEEN

#include "seq_common.h"

char seq_next_read_plain(SeqFile *sf);
char seq_read_base_plain(SeqFile *sf, char *c);
char seq_read_all_bases_plain(SeqFile *sf, StrBuf *sbuf);

size_t seq_file_write_seq_plain(SeqFile *sf, const char *seq);

#endif
