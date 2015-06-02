#!/bin/bash
#
g++ -c hello.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hello.cpp"
  exit
fi
#
g++ hello.o
if [ $? -ne 0 ]; then
  echo "Errors linking hello.o."
  exit
fi
#
rm hello.o
#
mv a.out hello
./hello > hello_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hello."
  exit
fi
rm hello
#
echo "Program output written to hello_output.txt"
