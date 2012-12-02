#!/bin/bash
#
g++ -c -g fsu_halton.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_halton.cpp"
  exit
fi
rm compiler.txt
#
mv fsu_halton.o ~/libcpp/$ARCH/fsu_halton.o
#
echo "A new version of fsu_halton has been created."
