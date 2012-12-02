#!/bin/bash
#
cp fsu.hpp /$HOME/include
#
g++ -c -g fsu_sub.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_sub.cpp"
  exit
fi
rm compiler.txt
#
mv fsu_sub.o ~/libcpp/$ARCH/fsu_sub.o
#
echo "A new version of fsu_sub has been created."
