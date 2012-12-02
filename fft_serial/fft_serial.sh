#!/bin/bash
#
g++ -c -g -I/$HOME/include fft_serial.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fft_serial.cpp"
  exit
fi
rm compiler.txt
#
g++ fft_serial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fft_serial.o."
  exit
fi
#
rm fft_serial.o
#
mv a.out fft_serial
./fft_serial > fft_serial_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fft_serial."
  exit
fi
rm fft_serial
#
echo "Program output written to fft_serial_output.txt"
