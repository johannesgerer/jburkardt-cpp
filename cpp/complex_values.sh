#!/bin/bash
#
g++ -c complex_values.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling complex_values.cpp"
  exit
fi
#
g++ complex_values.o
if [ $? -ne 0 ]; then
  echo "Errors linking complex_values.o."
  exit
fi
#
rm complex_values.o
#
mv a.out complex_values
./complex_values > complex_values_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running complex_values."
  exit
fi
rm complex_values
#
echo "Program output written to complex_values_output.txt"
