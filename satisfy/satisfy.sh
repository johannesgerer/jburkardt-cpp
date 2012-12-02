#!/bin/bash
#
g++ -c -g satisfy.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling satisfy.cpp"
  exit
fi
rm compiler.txt
#
g++ satisfy.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading satisfy.o."
  exit
fi
#
rm satisfy.o
#
mv a.out satisfy
./satisfy > satisfy_output.txt
rm satisfy
#
echo "Program output written to satisfy_output.txt"
