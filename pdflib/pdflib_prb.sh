#!/bin/bash
#
g++ -c -g -I/$HOME/include pdflib_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pdflib_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pdflib_prb.o /$HOME/libcpp/$ARCH/pdflib.o /$HOME/libcpp/$ARCH/rnglib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pdflib_prb.o"
  exit
fi
#
rm pdflib_prb.o
#
mv a.out pdflib_prb
./pdflib_prb > pdflib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pdflib_prb."
  exit
fi
rm pdflib_prb
#
echo "Program output written to pdflib_prb_output.txt"
