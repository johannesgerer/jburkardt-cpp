#!/bin/bash
#
g++ -c -I$HOME/include van_der_corput_dataset.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling van_der_corput_dataset.cpp"
  exit
fi
#
g++ van_der_corput_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading van_der_corput_dataset.o."
  exit
fi
#
rm van_der_corput_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/van_der_corput_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/van_der_corput_dataset"
