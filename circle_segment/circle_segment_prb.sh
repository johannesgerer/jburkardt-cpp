#!/bin/bash
#
g++ -c -g -I/$HOME/include circle_segment_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_segment_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ circle_segment_prb.o /$HOME/libcpp/$ARCH/circle_segment.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading circle_segment_prb.o."
  exit
fi
#
rm circle_segment_prb.o
#
mv a.out circle_segment_prb
./circle_segment_prb > circle_segment_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running circle_segment_prb."
  exit
fi
rm circle_segment_prb
#
echo "Program output written to circle_segment_prb_output.txt"
