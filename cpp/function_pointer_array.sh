#!/bin/bash
#
g++ -c function_pointer_array.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling function_pointer_array.cpp"
  exit
fi
rm compiler.txt
#
g++ function_pointer_array.o
if [ $? -ne 0 ]; then
  echo "Errors linking function_pointer_array.o."
  exit
fi
#
rm function_pointer_array.o
#
mv a.out function_pointer_array
./function_pointer_array > function_pointer_array_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running function_pointer_array."
  exit
fi
rm function_pointer_array
#
echo "Program output written to function_pointer_array_output.txt"
