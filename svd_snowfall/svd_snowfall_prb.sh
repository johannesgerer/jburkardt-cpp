#!/bin/bash
#
g++ -c -g -I/$HOME/include svd_snowfall_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling svd_snowfall_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ svd_snowfall_prb.o /$HOME/libcpp/$ARCH/svd_snowfall.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading svd_snowfall_prb.o"
  exit
fi
#
rm svd_snowfall_prb.o
#
mv a.out svd_snowfall_prb
./svd_snowfall_prb > svd_snowfall_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running svd_snowfall_prb."
  exit
fi
rm svd_snowfall_prb
#
echo "Program output written to svd_snowfall_prb_output.txt"
