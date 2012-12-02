#!/bin/bash
#
g++ -c -g search_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling search_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ search_serial.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading search_serial.o"
  exit
fi
rm search_serial.o
#
mv a.out search_serial
./search_serial > search_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running search_serial"
  exit
fi
rm search_serial
#
echo "Program output written to search_serial_output.txt"
