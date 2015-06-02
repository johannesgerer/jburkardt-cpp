#!/bin/bash
#
g++ -c function_pointer_array_new.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling function_pointer_array_new.cpp"
  exit
fi
#
g++ function_pointer_array_new.o
if [ $? -ne 0 ]; then
  echo "Errors linking function_pointer_array_new.o."
  exit
fi
#
rm function_pointer_array_new.o
#
mv a.out function_pointer_array_new
./function_pointer_array_new > function_pointer_array_new_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running function_pointer_array_new."
  exit
fi
rm function_pointer_array_new
#
echo "Program output written to function_pointer_array_new_output.txt"
