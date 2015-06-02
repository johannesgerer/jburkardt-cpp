#!/bin/bash
#
g++ -c array_return.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling array_return.cpp"
  exit
fi
#
g++ array_return.o
if [ $? -ne 0 ]; then
  echo "Errors loading array_return.o"
  exit
fi
rm array_return.o
#
mv a.out array_return
./array_return > array_return_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running array_return"
  exit
fi
rm array_return
#
echo "Program output written to array_return_output.txt"
