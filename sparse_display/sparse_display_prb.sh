#!/bin/bash
#
g++ -c -I/$HOME/include sparse_display_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_display_prb.cpp"
  exit
fi
#
g++ -o sparse_display_prb sparse_display_prb.o /$HOME/libcpp/$ARCH/sparse_display.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_display_prb.o."
  exit
fi
#
rm sparse_display_prb.o
#
./sparse_display_prb > sparse_display_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_display_prb."
  exit
fi
rm sparse_display_prb
#
echo "Program output written to sparse_display_prb_output.txt"
