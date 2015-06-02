#!/bin/bash
#
g++ -c -I/$HOME/include snakes_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling snakes_prb.cpp"
  exit
fi
#
g++ snakes_prb.o /$HOME/libcpp/$ARCH/snakes.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading snakes_prb.o"
  exit
fi
#
rm snakes_prb.o
#
mv a.out snakes_prb
./snakes_prb > snakes_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running snakes_prb."
  exit
fi
rm snakes_prb
#
echo "Program output written to snakes_prb_output.txt"
