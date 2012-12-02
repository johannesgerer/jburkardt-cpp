#!/bin/bash
#
g++ -c -g -I/$HOME/include point_merge_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling point_merge_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ point_merge_prb.o /$HOME/libcpp/$ARCH/point_merge.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading point_merge_prb.o."
  exit
fi
#
rm point_merge_prb.o
#
mv a.out point_merge_prb
./point_merge_prb > point_merge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running point_merge_prb."
  exit
fi
rm point_merge_prb
#
echo "Program output written to point_merge_prb_output.txt"
