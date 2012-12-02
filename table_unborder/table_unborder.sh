#!/bin/bash
#
g++ -c -I$HOME/include table_unborder.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_unborder.cpp"
  exit
fi
rm compiler.txt
#
g++ table_unborder.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_unborder.o + unborder.o."
  exit
fi
#
rm table_unborder.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/table_unborder
#
echo "Executable installed as ~/bincpp/$ARCH/table_unborder"
