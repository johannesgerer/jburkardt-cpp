#!/bin/bash
#
g++ -c -g -I/$HOME/include prime_serial_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling prime_serial_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ prime_serial_prb.o /$HOME/libcpp/$ARCH/prime_serial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading prime_serial_prb.o."
  exit
fi
#
rm prime_serial_prb.o
#
mv a.out prime_serial_prb
./prime_serial_prb > prime_serial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running prime_serial_prb."
  exit
fi
rm prime_serial_prb
#
echo "Program output written to prime_serial_prb_output.txt"
