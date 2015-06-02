#!/bin/bash
#
g++ -c fem1d_bvp_quadratic_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_bvp_quadratic_prb.cpp"
  exit
fi
#
g++ fem1d_bvp_quadratic_prb.o ~/libcpp/$ARCH/fem1d_bvp_quadratic.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_bvp_quadratic_prb.o"
  exit
fi
rm fem1d_bvp_quadratic_prb.o
#
mv a.out fem1d_bvp_quadratic_prb
./fem1d_bvp_quadratic_prb > fem1d_bvp_quadratic_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_bvp_quadratic_prb"
  exit
fi
rm fem1d_bvp_quadratic_prb
#
echo "Program output written to fem1d_bvp_quadratic_prb_output.txt"
