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
  if(argc != 2) {
    printf("usage: seq_reader_test <in.fa|fq|sam|bam>\n");
    return -1;
  }

  char* file_path = argv[1];

  printf(" Filepath: %s\n", file_path);

  Sequence* seq = seq_create();
  SequenceFile* file = seq_file_open(file_path);

  //SequenceFileType file_type = seq_file_get_type(file);
  printf(" Filetype: %s\n", seq_file_get_type_str(file));

  int i, j;
  
  for(i = 0; i < 5 && seq_file_read(file, seq); i++)
  {
    printf("%i>%s\n%s\n%s\n", i, seq->name->buff,
           seq->seq->buff, seq->qual->buff);
  }

  seq_file_close(file);

  // Printing kmers
  int num_of_kmers = 2;
  int kmers[] = {1,20};
  SequenceKmerReader* reader;

  for(i = 0; i < num_of_kmers; i++)
  {
    printf("kmer: %i\n", kmers[i]);
    reader = kmer_reader_open(file_path, kmers[i]);

    for(j = 0; j < 10 && kmer_reader_read_kmer(reader, seq); j++)
    {
      printf(" %i) %i:%i %s\n%s\n%s\n", j, (int)kmer_reader_entry_number(reader),
             seq->offset, seq->name->buff, seq->seq->buff, seq->qual->buff);
    }

    kmer_reader_close(reader);
  }

  seq_free(seq);

  return EXIT_SUCCESS;
}
