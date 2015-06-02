#!/bin/bash
#
g++ -c -I/$HOME/include zero_rc_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling zero_rc_prb.cpp"
  exit
fi
#
g++ zero_rc_prb.o /$HOME/libcpp/$ARCH/zero_rc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading zero_rc_prb.o."
  exit
fi
#
rm zero_rc_prb.o
#
mv a.out zero_rc_prb
./zero_rc_prb > zero_rc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running zero_rc_prb."
  exit
fi
rm zero_rc_prb
#
echo "Program output written to zero_rc_prb_output.txt"
