#!/bin/bash
#
g++ -c -g -I $HOME/include cvt_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ cvt_dataset.o $HOME/libcpp/$ARCH/cvt.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_dataset.o."
  exit
fi
#
rm cvt_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/cvt_dataset
#
echo "A new version of cvt_dataset has been created."
