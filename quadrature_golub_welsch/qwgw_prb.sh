#!/bin/bash
#
g++ -c -g -I/$HOME/include qwgw_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwgw_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ qwgw_prb.o /$HOME/libcpp/$ARCH/qwgw.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qwgw_prb.o."
  exit
fi
#
rm qwgw_prb.o
#
mv a.out qwgw_prb
./qwgw_prb > qwgw_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qwgw_prb."
  exit
fi
rm qwgw_prb
#
echo "Program output written to qwgw_prb_output.txt"
