#!/bin/bash
#
g++ -c -I$HOME/include table_quality.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_quality.cpp"
  exit
fi
rm compiler.txt
#
g++ table_quality.o $HOME/libcpp/$ARCH/quality.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_quality.o + quality.o."
  exit
fi
#
rm table_quality.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/table_quality
#
echo "Executable installed as ~/bincpp/$ARCH/table_quality"
