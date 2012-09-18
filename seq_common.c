
#include "seq_common.h"

// Array for complementing bases read from BAM/SAM files
int8_t seq_comp_table[16] = {0,8,4,12,2,10,9,14,1,6,5,13,3,11,7,15};

size_t _write_wrapped(SeqFile *sf, const char *str, size_t str_len)
{
  size_t num_bytes_printed = 0;

  if(sf->curr_line_length == sf->line_wrap)
  {
    sf->curr_line_length = 0;
    num_bytes_printed += seq_puts(sf, "\n");
    sf->line_number++;
  }
  
  if(sf->curr_line_length + str_len <= sf->line_wrap)
  {
    // Doesn't go over a single line
    sf->curr_line_length += str_len;
    num_bytes_printed += seq_puts(sf, str);
    return num_bytes_printed;
  }

  size_t bytes_to_print = sf->line_wrap - sf->curr_line_length;

  num_bytes_printed += seq_write(sf, str, bytes_to_print);
  num_bytes_printed += seq_puts(sf, "\n");
  sf->line_number++;

  size_t offset;

  for(offset = bytes_to_print; offset < str_len; offset += sf->line_wrap)
  {
    bytes_to_print = MIN(str_len - offset, sf->line_wrap);
    num_bytes_printed += (size_t)seq_write(sf, str + offset, bytes_to_print);

    if(bytes_to_print < sf->line_wrap)
    {
      num_bytes_printed += seq_puts(sf, "\n");
      sf->line_number++;
    }
  }

  sf->curr_line_length = bytes_to_print;

  return num_bytes_printed;
}

