#!/bin/bash
#
g++ -c -g -I$HOME/include niederreiter2_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter2_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ niederreiter2_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading niederreiter2_dataset.o."
  exit
fi
#
rm niederreiter2_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/niederreiter2_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/niederreiter2_dataset"
