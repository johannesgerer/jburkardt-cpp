#!/bin/bash
#
g++ -c -g -I/$HOME/include fem1d_pack_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_pack_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ fem1d_pack_prb.o /$HOME/libcpp/$ARCH/fem1d_pack.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_pack_prb.o."
  exit
fi
#
rm fem1d_pack_prb.o
#
mv a.out fem1d_pack_prb
./fem1d_pack_prb > fem1d_pack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_pack_prb."
  exit
fi
rm fem1d_pack_prb
#
echo "Program output written to fem1d_pack_prb_output.txt"
