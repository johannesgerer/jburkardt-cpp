#!/bin/bash
#
g++ -c -g fsu_latinize.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_latinize.cpp"
  exit
fi
rm compiler.txt
#
mv fsu_latinize.o ~/libcpp/$ARCH/fsu_latinize.o
#
echo "A new version of fsu_latinize has been created."
