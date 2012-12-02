#!/bin/bash
#
g++ -c -g -I/$HOME/include latin_edge_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_edge_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ latin_edge_prb.o /$HOME/libcpp/$ARCH/latin_edge.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_edge_prb.o."
  exit
fi
#
rm latin_edge_prb.o
#
mv a.out latin_edge_prb
./latin_edge_prb > latin_edge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latin_edge_prb."
  exit
fi
rm latin_edge_prb
#
echo "Program output written to latin_edge_prb_output.txt"
