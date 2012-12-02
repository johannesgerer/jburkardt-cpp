#!/bin/bash
#
g++ -c -g -I/$HOME/include crc_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling crc_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ crc_prb.o /$HOME/libcpp/$ARCH/crc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading crc_prb.o."
  exit
fi
#
rm crc_prb.o
#
mv a.out crc_prb
./crc_prb > crc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running crc_prb."
  exit
fi
rm crc_prb
#
echo "Program output written to crc_prb_output.txt"
