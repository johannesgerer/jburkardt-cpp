#!/bin/bash
#
g++ -c -I/$HOME/include ns2de_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ns2de_prb.cpp"
  exit
fi
#
g++ -o ns2de_prb ns2de_prb.o /$HOME/libcpp/$ARCH/ns2de.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ns2de_prb.o."
  exit
fi
#
rm ns2de_prb.o
#
./ns2de_prb > ns2de_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ns2de_prb."
  exit
fi
rm ns2de_prb
#
echo "Program output written to ns2de_prb_output.txt"
