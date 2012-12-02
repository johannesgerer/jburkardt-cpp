#!/bin/bash
#
g++ -c -g fsu_quality.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_quality.cpp"
  exit
fi
rm compiler.txt
#
mv fsu_quality.o ~/libcpp/$ARCH/fsu_quality.o
#
echo "A new version of fsu_quality has been created."
