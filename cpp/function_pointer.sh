#!/bin/bash
#
g++ -c function_pointer.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling function_pointer.cpp"
  exit
fi
#
g++ function_pointer.o
if [ $? -ne 0 ]; then
  echo "Errors linking function_pointer.o."
  exit
fi
#
rm function_pointer.o
#
mv a.out function_pointer
./function_pointer > function_pointer_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running function_pointer."
  exit
fi
rm function_pointer
#
echo "Program output written to function_pointer_output.txt"
