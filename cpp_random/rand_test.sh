#!/bin/bash
#
g++ -c rand_test.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rand_test.cpp"
  exit
fi
rm compiler.txt
#
g++ rand_test.o
if [ $? -ne 0 ]; then
  echo "Errors loading rand_test.o"
  exit
fi
rm rand_test.o
#
mv a.out rand_test
./rand_test > rand_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rand_test"
  exit
fi
rm rand_test
#
echo "Program output written to rand_test_output.txt"
