#!/bin/bash
#
g++ -c -I/$HOME/include s2de_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling s2de_prb.cpp"
  exit
fi
#
g++ -o s2de_prb s2de_prb.o /$HOME/libcpp/$ARCH/s2de.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading s2de_prb.o."
  exit
fi
#
rm s2de_prb.o
#
./s2de_prb > s2de_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running s2de_prb."
  exit
fi
rm s2de_prb
#
echo "Program output written to s2de_prb_output.txt"
