#!/bin/bash
#
g++ -c -g -I/$HOME/include doomsday_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling doomsday_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ doomsday_prb.o /$HOME/libcpp/$ARCH/doomsday.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading doomsday_prb.o."
  exit
fi
#
rm doomsday_prb.o
#
mv a.out doomsday_prb
./doomsday_prb > doomsday_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running doomsday_prb."
  exit
fi
rm doomsday_prb
#
echo "Program output written to doomsday_prb_output.txt"
