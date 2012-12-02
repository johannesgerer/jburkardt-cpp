#!/bin/bash
#
g++ -c poisson_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ poisson_serial.o
if [ $? -ne 0 ]; then
  echo "Errors loading poisson_serial.o"
  exit
fi
rm poisson_serial.o
#
mv a.out poisson_serial
./poisson_serial > poisson_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson_serial"
  exit
fi
rm poisson_serial
#
echo "The poisson_serial test problem has been executed."
