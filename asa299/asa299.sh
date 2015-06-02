#!/bin/bash
#
cp asa299.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa299.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa299.cpp"
  exit
fi
#
mv asa299.o ~/libcpp/$ARCH/asa299.o
#
echo "Library installed as ~/libcpp/$ARCH/asa299.o"
