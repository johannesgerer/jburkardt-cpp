#!/bin/bash
#
g++ -c dynamic_array_2d.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling dynamic_array_2d.cpp"
  exit
fi
#
g++ dynamic_array_2d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dynamic_array_2d.o."
  exit
fi
#
rm dynamic_array_2d.o
#
chmod ugo+x a.out
mv a.out dynamic_array_2d
dynamic_array_2d > dynamic_array_2d_output.txt
rm dynamic_array_2d
#
echo "Program output written to dynamic_array_2d_output.txtvvvvv"
