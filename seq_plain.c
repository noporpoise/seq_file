/*
 seq_plain.c
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

#include "seq_plain.h"

char seq_next_read_plain(SeqFile *sf)
{
  if(sf->read_line_start)
  {
    sf->read_line_start = 0;

    int c;

    while((c = gzgetc(sf->gz_file)) != -1 && c != '\r' && c != '\n')
    {
      sf->total_bases_skipped++;
    }

    if(c == -1)
      return 0;
    else
      sf->line_number++;

    return 1;
  }
  else
  {
    // Check if we can read a base
    int c = gzgetc(sf->gz_file);

    if(c == -1)
    {
      return 0;
    }
    else if(c == '\n' || c == '\r')
    {
      sf->line_number++;
    }
    else
    {
      sf->read_line_start = (char)c;
    }

    return 1;
  }
}

char seq_read_base_plain(SeqFile *sf, char *c)
{
  if(sf->read_line_start != 0)
  {
    *c = sf->read_line_start;
    sf->read_line_start = 0;
    return 1;
  }

  int next = gzgetc(sf->gz_file);

  if(next == -1)
  {
    return 0;
  }
  else if(next == '\r' || next == '\n')
  {
    sf->line_number++;
    return 0;
  }
  else
  {
    *c = (char)next;
    return 1;
  }
}

char seq_read_all_bases_plain(SeqFile *sf, StrBuf *sbuf)
{
  t_buf_pos len = 0;

  if(sf->read_line_start)
  {
    strbuf_append_char(sbuf, sf->read_line_start);
    sf->read_line_start = 0;
    len++;
  }

  len += strbuf_gzreadline(sbuf, sf->gz_file);
  strbuf_chomp(sbuf);
  sf->line_number++;

  return (len > 0);
}

size_t seq_file_write_seq_plain(SeqFile *sf, const char *seq)
{
  size_t num_bytes_printed = 0;

  if(sf->write_state != WS_BEGIN)
  {
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }

  num_bytes_printed += seq_puts(sf, seq);

  return num_bytes_printed;
}
