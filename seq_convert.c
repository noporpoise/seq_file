/*
 seq_reader.c
 project: seq_reader
 author: Isaac Turner <turner.isaac@gmail.com>
 Copyright (C) 21 June 2012

 To build type 'make'

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

#include <stdlib.h>
#include <stdio.h>

#include "string_buffer.h"
#include "seq_reader.h"

int main(int argc, char** argv)
{
  if(argc != 3)
  {
    printf("usage: seq_convert <in> <out>\n"
"  Output file must end with: fa|fq|fa.gz|fq.gz|sam|bam\n");
    return -1;
  }

  printf("zlib version: %s\n", ZLIB_VERSION);

  char* in_path = argv[1];
  char* out_path = argv[2];

  SequenceFileType out_file_type = SEQ_UNKNOWN;
  char out_zipped = 0;

  seq_file_guess_filetype_from_path(out_path, &out_file_type, &out_zipped);

  if(out_file_type == SEQ_UNKNOWN)
  {
    fprintf(stderr, "Sorry, I cannot identify the output file's format "
                    "from its path\n");
    exit(EXIT_FAILURE);
  }

  SequenceFile* in_file = seq_file_open(in_path);
  SequenceFile* out_file = seq_file_open_write(out_path, out_file_type,
                                               out_zipped, 0);

  if(in_file == NULL)
  {
    fprintf(stderr, "Couldn't open input file: %s\n", in_path);
    exit(EXIT_FAILURE);
  }

  if(out_file == NULL)
  {
    fprintf(stderr, "Couldn't open output file: %s\n", out_path);
    exit(EXIT_FAILURE);
  }

  printf(" In : %s [%s]\n", in_path, seq_file_get_type_str(in_file));
  printf(" Out: %s [%s]\n", out_path,
         seq_file_convert_type_str(seq_file_get_type(out_file), out_zipped));

  // Start converting
  Sequence* seq = seq_create();

  unsigned long bytes_written = 0;

  while(seq_file_read(in_file, seq))
  {
    if(!(bytes_written += seq_file_write(out_file, seq)))
    {
      fprintf(stderr, "Couldn't write to file\n");
      break;
    }
  }

  seq_free(seq);

  seq_file_close(in_file);
  seq_file_close(out_file);

  printf("%lu bytes written\n", bytes_written);

  printf("Done. \n");

  return EXIT_SUCCESS;
}
