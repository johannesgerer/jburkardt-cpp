#!/bin/bash
#
g++ -c mgs_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling mgs_prb.cpp"
  exit
fi
#
g++ mgs_prb.o mgs.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mgs_prb.o"
  exit
fi
rm mgs_prb.o
#
mv a.out mgs_prb
./mgs_prb > mgs_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mgs_prb"
  exit
fi
rm mgs_prb
#
echo "The program output was written to mgs_prb_output.txt"
