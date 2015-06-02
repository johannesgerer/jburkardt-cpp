#!/bin/bash
#
g++ -c -I/$HOME/include line_cvt_lloyd_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_cvt_lloyd_prb.cpp"
  exit
fi
#
g++ -o line_cvt_lloyd_prb line_cvt_lloyd_prb.o /$HOME/libcpp/$ARCH/line_cvt_lloyd.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_cvt_lloyd_prb.o."
  exit
fi
#
rm line_cvt_lloyd_prb.o
#
./line_cvt_lloyd_prb > line_cvt_lloyd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_cvt_lloyd_prb."
  exit
fi
rm line_cvt_lloyd_prb
#
echo "Program output written to line_cvt_lloyd_prb_output.txt"
