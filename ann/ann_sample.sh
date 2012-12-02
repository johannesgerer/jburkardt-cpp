#!/bin/bash
#
g++ -c -g -I/$HOME/include ann_sample.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ann_sample.cpp"
  exit
fi
rm compiler.txt
#
g++ ann_sample.o -L/$HOME/libcpp/$ARCH -lann -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ann_sample.o."
  exit
fi
#
rm ann_sample.o
#
mv a.out ~/bincpp/$ARCH/ann_sample
#
echo "Executable installed as ~/bincpp/$ARCH/ann_sample"
