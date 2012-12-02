#!/bin/bash
#
g++ -c -g -I$HOME/include normal_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ normal_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading normal_dataset.o."
  exit
fi
#
rm normal_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/normal_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/normal_dataset"
