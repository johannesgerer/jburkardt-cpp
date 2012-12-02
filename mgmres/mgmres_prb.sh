#!/bin/bash
#
g++ -c -g -I/$HOME/include mgmres_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mgmres_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ mgmres_prb.o /$HOME/libcpp/$ARCH/mgmres.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mgmres_prb.o."
  exit
fi
#
rm mgmres_prb.o
#
mv a.out mgmres_prb
./mgmres_prb > mgmres_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mgmres_prb."
  exit
fi
rm mgmres_prb
#
echo "Program output written to mgmres_prb_output.txt"
