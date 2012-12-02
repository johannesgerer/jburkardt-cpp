#!/bin/bash
#
g++ -c -g -I/$HOME/include fekete_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fekete_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ fekete_prb.o /$HOME/libcpp/$ARCH/fekete.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fekete_prb.o."
  exit
fi
#
rm fekete_prb.o
#
mv a.out fekete_prb
./fekete_prb > fekete_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fekete_prb."
  exit
fi
rm fekete_prb
#
echo "Program output written to fekete_prb_output.txt"
