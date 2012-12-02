#!/bin/bash
#
g++ -c -g -I/$HOME/include uniform_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ uniform_prb.o /$HOME/libcpp/$ARCH/uniform.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading uniform_prb.o."
  exit
fi
#
rm uniform_prb.o
#
mv a.out uniform_prb
./uniform_prb > uniform_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running uniform_prb."
  exit
fi
rm uniform_prb
#
echo "Program output written to uniform_prb_output.txt"
