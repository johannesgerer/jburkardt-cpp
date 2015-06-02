#!/bin/bash
#
g++ -c -I/$HOME/include hammersley_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hammersley_prb.cpp"
  exit
fi
#
g++ hammersley_prb.o /$HOME/libcpp/$ARCH/hammersley.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hammersley_prb.o."
  exit
fi
#
rm hammersley_prb.o
#
mv a.out hammersley_prb
./hammersley_prb > hammersley_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hammersley_prb."
  exit
fi
rm hammersley_prb
#
echo "Program output written to hammersley_prb_output.txt"
