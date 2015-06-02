#!/bin/bash
#
g++ -c -I/$HOME/include sort_rc_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling sort_rc_prb.cpp"
  exit
fi
#
g++ -o sort_rc_prb sort_rc_prb.o /$HOME/libcpp/$ARCH/sort_rc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sort_rc_prb.o."
  exit
fi
#
rm sort_rc_prb.o
#
./sort_rc_prb > sort_rc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sort_rc_prb."
  exit
fi
rm sort_rc_prb
#
echo "Program output written to sort_rc_prb_output.txt"
