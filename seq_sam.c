/*
 seq_sam.c
 project: seq_file
 url: https://github.com/noporpoise/seq_file
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 20-June-2012
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "seq_sam.h"

char seq_next_read_sam(SeqFile *sf)
{
  if(sf->entry_read)
  {
    // Count skipped bases
    sf->total_bases_skipped += (unsigned long)sf->bam->core.l_qseq -
                               sf->entry_offset;
  }

  if(samread(sf->sam_file, sf->bam) < 0)
    return 0;

  // Get name
  strbuf_append_str(sf->entry_name, bam1_qname(sf->bam));

  return 1;
}

char seq_read_base_sam(SeqFile *sf, char *c)
{
  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;
  unsigned long query_len = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= query_len)
    return 0;

  uint8_t *seq = bam1_seq(sf->bam);

  if(is_reversed)
  {
    unsigned long index = query_len - sf->entry_offset - 1;
    int8_t b = bam1_seqi(seq, index);
    *c = bam_nt16_rev_table[seq_comp_table[b]];
  }
  else
  {
    int8_t b = bam1_seqi(seq, sf->entry_offset);
    *c = bam_nt16_rev_table[b];
  }

  return 1;
}

char seq_read_qual_sam(SeqFile *sf, char *c)
{
  char is_reversed = sf->bam->core.flag & 16;
  unsigned long query_len = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= query_len)
    return 0;

  uint8_t *seq = bam1_qual(sf->bam);
  unsigned long index;

  if(is_reversed)
    index = query_len - sf->entry_offset_qual - 1;
  else
    index = sf->entry_offset_qual;

  *c = sf->fastq_ascii_offset + seq[index];

  return 1;
}

char seq_read_all_bases_sam(SeqFile *sf, StrBuf *sbuf)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset >= qlen)
    return 0;

  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;

  strbuf_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam1_seq(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    unsigned long index = (is_reversed ? i : qlen - i - 1);
    int8_t b = bam1_seqi(seq, index);
    char c = bam_nt16_rev_table[is_reversed ? seq_comp_table[b] : b];
    strbuf_append_char(sbuf, c);
  }

  return 1;
}

char seq_read_all_quals_sam(SeqFile *sf, StrBuf *sbuf)
{
  unsigned long qlen = (unsigned long)sf->bam->core.l_qseq;

  if(sf->entry_offset_qual >= qlen)
    return 0;

  // Get reverse
  char is_reversed = sf->bam->core.flag & 16;

  strbuf_ensure_capacity(sbuf, qlen - sf->entry_offset);

  uint8_t *seq = bam1_qual(sf->bam);

  // read in and reverse complement (if needed)
  unsigned long i;
  for(i = sf->entry_offset; i < qlen; i++)
  {
    char c = sf->fastq_ascii_offset + seq[is_reversed ? i : qlen - i - 1];
    strbuf_append_char(sbuf, c);
  }

  return 1;
}
