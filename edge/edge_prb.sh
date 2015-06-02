#!/bin/bash
#
g++ -c -I/$HOME/include edge_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling edge_prb.cpp"
  exit
fi
#
g++ edge_prb.o /$HOME/libcpp/$ARCH/edge.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading edge_prb.o."
  exit
fi
#
rm edge_prb.o
#
mv a.out edge_prb
./edge_prb > edge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running edge_prb."
  exit
fi
rm edge_prb
#
echo "Program output written to edge_prb_output.txt"
