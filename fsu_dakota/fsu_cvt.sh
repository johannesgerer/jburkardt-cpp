#!/bin/bash
#
g++ -c -g fsu_cvt.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_cvt.cpp"
  exit
fi
rm compiler.txt
#
mv fsu_cvt.o ~/libcpp/$ARCH/fsu_cvt.o
#
echo "A new version of fsu_cvt has been created."
