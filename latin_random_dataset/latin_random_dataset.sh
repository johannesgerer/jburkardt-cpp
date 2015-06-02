#!/bin/bash
#
g++ -c -I$HOME/include latin_random_dataset.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_random_dataset.cpp"
  exit
fi
#
g++ latin_random_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_random_dataset.o."
  exit
fi
#
rm latin_random_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/latin_random_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/latin_random_dataset"
