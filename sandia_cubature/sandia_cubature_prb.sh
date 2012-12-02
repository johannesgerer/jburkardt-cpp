#!/bin/bash
#
g++ -c -g -I/$HOME/include sandia_cubature_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_cubature_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sandia_cubature_prb.o /$HOME/libcpp/$ARCH/sandia_cubature.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_cubature_prb.o."
  exit
fi
#
rm sandia_cubature_prb.o
#
mv a.out sandia_cubature_prb
./sandia_cubature_prb > sandia_cubature_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_cubature_prb."
  exit
fi
rm sandia_cubature_prb
#
echo "Program output written to sandia_cubature_prb_output.txt"
