#!/bin/bash
#
g++ -c -I$HOME/include uniform_dataset.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform_dataset.cpp"
  exit
fi
#
g++ uniform_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading uniform_dataset.o."
  exit
fi
#
rm uniform_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/uniform_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/uniform_dataset"
