#!/bin/bash
#
g++ -c -I/$HOME/include hpp_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hpp_prb.cpp"
  exit
fi
#
g++ hpp_prb.o /$HOME/libcpp/$ARCH/hpp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hpp_prb.o"
  exit
fi
#
rm hpp_prb.o
#
mv a.out hpp_prb
./hpp_prb > hpp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hpp_prb."
  exit
fi
rm hpp_prb
#
echo "Program output written to hpp_prb_output.txt"
