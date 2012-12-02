#!/bin/bash
#
g++ -c -g -I/$HOME/include beta_nc_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling beta_nc_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ beta_nc_prb.o /$HOME/libcpp/$ARCH/beta_nc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading beta_nc_prb.o."
  exit
fi
#
rm beta_nc_prb.o
#
mv a.out beta_nc_prb
./beta_nc_prb > beta_nc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running beta_nc_prb."
  exit
fi
rm beta_nc_prb
#
echo "Program output written to beta_nc_prb_output.txt"
