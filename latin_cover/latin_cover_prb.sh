#!/bin/bash
#
g++ -c -g -I/$HOME/include latin_cover_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_cover_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ latin_cover_prb.o /$HOME/libcpp/$ARCH/latin_cover.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_cover_prb.o"
  exit
fi
#
rm latin_cover_prb.o
#
mv a.out latin_cover_prb
./latin_cover_prb > latin_cover_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latin_cover_prb."
  exit
fi
rm latin_cover_prb
#
echo "Program output written to latin_cover_prb_output.txt"
