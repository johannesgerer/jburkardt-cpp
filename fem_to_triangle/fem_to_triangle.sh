#!/bin/bash
#
g++ -c fem_to_triangle.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_triangle.cpp"
  exit
fi
#
g++ fem_to_triangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_triangle.o."
  exit
fi
#
rm fem_to_triangle.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem_to_triangle
#
echo "Executable installed as ~/bincpp/$ARCH/fem_to_triangle"
