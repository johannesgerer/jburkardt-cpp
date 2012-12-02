#!/bin/bash
#
cp test_mat.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_mat.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_mat.cpp"
  exit
fi
rm compiler.txt
#
mv test_mat.o ~/libcpp/$ARCH/test_mat.o
#
echo "Library installed as ~/libcpp/$ARCH/test_mat.o"
