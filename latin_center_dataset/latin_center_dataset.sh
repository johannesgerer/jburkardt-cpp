#!/bin/bash
#
g++ -c -g -I$HOME/include latin_center_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_center_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ latin_center_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_center_dataset.o."
  exit
fi
#
rm latin_center_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/latin_center_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/latin_center_dataset"
