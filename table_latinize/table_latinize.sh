#!/bin/bash
#
g++ -c -I$HOME/include table_latinize.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_latinize.cpp"
  exit
fi
rm compiler.txt
#
g++ table_latinize.o $HOME/libcpp/$ARCH/latinize.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_latinize.o + latinize.o."
  exit
fi
#
rm table_latinize.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/table_latinize
#
echo "Executable installed as ~/bincpp/$ARCH/table_latinize"
