#!/bin/bash
#
g++ -c -g file_name_sequence_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling file_name_sequence_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ file_name_sequence_prb.o /$HOME/libcpp/$ARCH/file_name_sequence.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading file_name_sequence_prb.o."
  exit
fi
#
rm file_name_sequence_prb.o
#
mv a.out file_name_sequence_prb
./file_name_sequence_prb > file_name_sequence_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running file_name_sequence_prb."
  exit
fi
rm file_name_sequence_prb
#
echo "Program output written to file_name_sequence_prb_output.txt"
