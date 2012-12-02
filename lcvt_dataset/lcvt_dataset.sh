#!/bin/bash
#
g++ -c -g -I $HOME/include lcvt_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lcvt_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ lcvt_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lcvt_dataset.o."
  exit
fi
#
rm lcvt_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/lcvt_dataset
#
echo "A new version of lcvt_dataset has been created."
