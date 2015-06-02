#!/bin/bash
#
g++ -c -g fire_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fire_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ fire_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fire_serial.o"
  exit
fi
rm fire_serial.o
#
mv a.out fire_serial
./fire_serial > fire_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fire_serial"
  exit
fi
rm fire_serial
#
echo "Program output written to fire_serial_output.txt"
