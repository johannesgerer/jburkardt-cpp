#!/bin/bash
#
g++ -c -I$HOME/include table_border.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling table_border.cpp"
  exit
fi
rm compiler.txt
#
g++ table_border.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading table_border.o + border.o."
  exit
fi
#
rm table_border.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/table_border
#
echo "Executable installed as ~/bincpp/$ARCH/table_border"
