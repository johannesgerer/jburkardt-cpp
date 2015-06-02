#!/bin/bash
#
g++ -c -I/$HOME/include machine_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling machine_prb.cpp"
  exit
fi
#
g++ machine_prb.o /$HOME/libcpp/$ARCH/machine.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading machine_prb.o."
  exit
fi
#
rm machine_prb.o
#
mv a.out machine_prb
./machine_prb > machine_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running machine_prb."
  exit
fi
rm machine_prb
#
echo "Program output written to machine_prb_output.txt"
