#!/bin/bash
#
g++ -c -I/$HOME/include cube_exactness_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_exactness_prb.cpp"
  exit
fi
#
g++ -o cube_exactness_prb cube_exactness_prb.o /$HOME/libcpp/$ARCH/cube_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_exactness_prb.o."
  exit
fi
#
rm cube_exactness_prb.o
#
./cube_exactness_prb > cube_exactness_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_exactness_prb."
  exit
fi
rm cube_exactness_prb
#
echo "Program output written to cube_exactness_prb_output.txt"
