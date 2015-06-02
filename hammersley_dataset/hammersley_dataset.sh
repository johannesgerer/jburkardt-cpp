#!/bin/bash
#
g++ -c -I$HOME/include hammersley_dataset.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hammersley_dataset.cpp"
  exit
fi
#
g++ hammersley_dataset.o $HOME/libcpp/$ARCH/hammersley.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hammersley_dataset.o."
  exit
fi
#
rm hammersley_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hammersley_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/hammersley_dataset"
