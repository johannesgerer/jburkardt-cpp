#!/bin/bash
#
g++ -c -I/$HOME/include toms443_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling toms443_prb.cpp"
  exit
fi
#
g++ -o toms443_prb toms443_prb.o /$HOME/libcpp/$ARCH/toms443.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms443_prb.o."
  exit
fi
#
rm toms443_prb.o
#
./toms443_prb > toms443_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms443_prb."
  exit
fi
rm toms443_prb
#
echo "Program output written to toms443_prb_output.txt"
