#!/bin/bash
#
g++ -c -g -I/$HOME/include dcdflib_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dcdflib_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ dcdflib_prb.o /$HOME/libcpp/$ARCH/dcdflib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dcdflib_prb.o."
  exit
fi
#
rm dcdflib_prb.o
#
mv a.out dcdflib_prb
./dcdflib_prb > dcdflib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dcdflib_prb."
  exit
fi
rm dcdflib_prb
#
echo "Program output written to dcdflib_prb_output.txt"
