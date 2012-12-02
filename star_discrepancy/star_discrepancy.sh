#!/bin/bash
#
g++ -c star_discrepancy.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling star_discrepancy.cpp"
  exit
fi
rm compiler.txt
#
g++ star_discrepancy.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading star_discrepancy.o."
  exit
fi
#
rm star_discrepancy.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/star_discrepancy
#
echo "Execustar installed as ~/bincpp/$ARCH/star_discrepancy"
