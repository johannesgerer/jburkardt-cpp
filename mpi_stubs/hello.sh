#!/bin/bash
#
g++ -c hello.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling hello.cpp"
  exit
fi
rm compiler.txt
#
g++ hello.o $HOME/libcpp/$ARCH/mpi_stubs.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading hello.o"
  exit
fi
rm hello.o
#
mv a.out hello
./hello > hello_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running hello"
  exit
fi
rm hello
#
echo "Program output written to hello_output.txt"
