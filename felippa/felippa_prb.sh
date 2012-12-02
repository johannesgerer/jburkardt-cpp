#!/bin/bash
#
g++ -c -g -I/$HOME/include felippa_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling felippa_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ felippa_prb.o /$HOME/libcpp/$ARCH/felippa.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading felippa_prb.o."
  exit
fi
#
rm felippa_prb.o
#
mv a.out felippa_prb
./felippa_prb > felippa_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running felippa_prb."
  exit
fi
rm felippa_prb
#
echo "Program output written to felippa_prb_output.txt"
