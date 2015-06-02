#!/bin/bash
#
g++ -c fem2d_bvp_linear_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_linear_prb.cpp"
  exit
fi
#
g++ fem2d_bvp_linear_prb.o ~/libcpp/$ARCH/fem2d_bvp_linear.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_bvp_linear_prb.o"
  exit
fi
rm fem2d_bvp_linear_prb.o
#
mv a.out fem2d_bvp_linear_prb
./fem2d_bvp_linear_prb > fem2d_bvp_linear_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_bvp_linear_prb"
  exit
fi
rm fem2d_bvp_linear_prb
#
echo "Program output written to fem2d_bvp_linear_prb_output.txt"
