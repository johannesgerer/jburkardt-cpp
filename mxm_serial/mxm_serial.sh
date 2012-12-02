#!/bin/bash
#
g++ -c -g -I/$HOME/include mxm_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxm_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ mxm_serial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxm_serial.o."
  exit
fi
#
rm mxm_serial.o
#
mv a.out mxm_serial
./mxm_serial > mxm_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mxm_serial."
  exit
fi
rm mxm_serial
#
echo "Program output written to mxm_serial_output.txt"
