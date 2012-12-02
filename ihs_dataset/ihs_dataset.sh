#!/bin/bash
#
g++ -c -g -I $HOME/include ihs_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ihs_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ ihs_dataset.o $HOME/libcpp/$ARCH/ihs.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ihs_dataset.o."
  exit
fi
#
rm ihs_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ihs_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/ihs_dataset"
