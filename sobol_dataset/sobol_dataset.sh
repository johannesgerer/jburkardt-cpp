#!/bin/bash
#
g++ -c -g -I$HOME/include sobol_dataset.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sobol_dataset.cpp"
  exit
fi
rm compiler.txt
#
g++ sobol_dataset.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sobol_dataset.o."
  exit
fi
#
rm sobol_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/sobol_dataset
#
echo "Executable installed as ~/bincpp/$ARCH/sobol_dataset"
