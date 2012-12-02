#!/bin/bash
#
g++ -c -g -I/$HOME/include calpak_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling calpak_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ calpak_prb.o /$HOME/libcpp/$ARCH/calpak.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading calpak_prb.o."
  exit
fi
#
rm calpak_prb.o
#
mv a.out calpak_prb
./calpak_prb > calpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running calpak_prb."
  exit
fi
rm calpak_prb
#
echo "Program output written to calpak_prb_output.txt"
