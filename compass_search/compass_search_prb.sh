#!/bin/bash
#
g++ -c -g -I/$HOME/include compass_search_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling compass_search_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ compass_search_prb.o /$HOME/libcpp/$ARCH/compass_search.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading compass_search_prb.o."
  exit
fi
#
rm compass_search_prb.o
#
mv a.out compass_search_prb
./compass_search_prb > compass_search_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running compass_search_prb."
  exit
fi
rm compass_search_prb
#
echo "Program output written to compass_search_prb_output.txt"
