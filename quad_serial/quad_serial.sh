#!/bin/bash
#
g++ -c quad_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ quad_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking quad_serial.o."
  exit
fi
#
rm quad_serial.o
#
mv a.out quad_serial
./quad_serial > quad_serial_output.txt
rm quad_serial
#
echo "Program output written to quad_serial_output.txt"
