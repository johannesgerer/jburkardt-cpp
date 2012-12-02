#!/bin/bash
#
g++ -c -g -I/$HOME/include random_data_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_data_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ random_data_prb.o /$HOME/libcpp/$ARCH/random_data.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading random_data_prb.o."
  exit
fi
#
rm random_data_prb.o
#
mv a.out random_data_prb
./random_data_prb > random_data_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running random_data_prb."
  exit
fi
rm random_data_prb
#
echo "Program output written to random_data_prb_output.txt"
