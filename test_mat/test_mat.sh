#!/bin/bash
#
cp test_mat.hpp /$HOME/include
#
g++ -c -I /$HOME/include test_mat.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling test_mat.cpp"
  exit
fi
#
mv test_mat.o ~/libcpp/$ARCH/test_mat.o
#
echo "Library installed as ~/libcpp/$ARCH/test_mat.o"
