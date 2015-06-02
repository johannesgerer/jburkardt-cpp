#!/bin/bash
#
g++ -c -I/$HOME/include ihs_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ihs_prb.cpp"
  exit
fi
#
g++ ihs_prb.o /$HOME/libcpp/$ARCH/ihs.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ihs_prb.o."
  exit
fi
#
rm ihs_prb.o
#
mv a.out ihs_prb
./ihs_prb > ihs_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ihs_prb."
  exit
fi
rm ihs_prb
#
echo "Program output written to ihs_prb_output.txt"
