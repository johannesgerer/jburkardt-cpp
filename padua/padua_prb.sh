#!/bin/bash
#
g++ -c -I/$HOME/include padua_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling padua_prb.cpp"
  exit
fi
#
g++ padua_prb.o /$HOME/libcpp/$ARCH/padua.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading padua_prb.o."
  exit
fi
#
rm padua_prb.o
#
mv a.out padua_prb
./padua_prb > padua_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running padua_prb."
  exit
fi
rm padua_prb
#
echo "Program output written to padua_prb_output.txt"
