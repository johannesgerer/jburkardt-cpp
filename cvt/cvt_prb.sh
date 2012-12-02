#!/bin/bash
#
g++ -c -g -I/$HOME/include cvt_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ cvt_prb.o /$HOME/libcpp/$ARCH/cvt.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_prb.o."
  exit
fi
#
rm cvt_prb.o
#
mv a.out cvt_prb
./cvt_prb > cvt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_prb."
  exit
fi
rm cvt_prb
#
echo "Program output written to cvt_prb_output.txt"
