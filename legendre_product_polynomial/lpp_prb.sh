#!/bin/bash
#
g++ -c -I/$HOME/include lpp_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling lpp_prb.cpp"
  exit
fi
#
g++ lpp_prb.o /$HOME/libcpp/$ARCH/lpp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lpp_prb.o"
  exit
fi
#
rm lpp_prb.o
#
mv a.out lpp_prb
./lpp_prb > lpp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lpp_prb."
  exit
fi
rm lpp_prb
#
echo "Program output written to lpp_prb_output.txt"
