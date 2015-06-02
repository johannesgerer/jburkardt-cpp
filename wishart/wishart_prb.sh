#!/bin/bash
#
g++ -c -I/$HOME/include wishart_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wishart_prb.cpp"
  exit
fi
#
g++ wishart_prb.o \
  /$HOME/libcpp/$ARCH/wishart.o \
  /$HOME/libcpp/$ARCH/pdflib.o \
  /$HOME/libcpp/$ARCH/rnglib.o -lm

if [ $? -ne 0 ]; then
  echo "Errors linking and loading wishart_prb.o."
  exit
fi
#
rm wishart_prb.o
#
mv a.out wishart_prb
./wishart_prb > wishart_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wishart_prb."
  exit
fi
rm wishart_prb
#
echo "Program output written to wishart_prb_output.txt"
