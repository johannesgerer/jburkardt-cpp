#!/bin/bash
#
g++ -c -I/$HOME/include bellman_ford_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling bellman_ford_prb.cpp"
  exit
fi
#
g++ bellman_ford_prb.o /$HOME/libcpp/$ARCH/bellman_ford.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bellman_ford_prb.o"
  exit
fi
#
rm bellman_ford_prb.o
#
mv a.out bellman_ford_prb
./bellman_ford_prb > bellman_ford_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bellman_ford_prb."
  exit
fi
rm bellman_ford_prb
#
echo "Program output written to bellman_ford_prb_output.txt"
