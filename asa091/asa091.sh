#!/bin/bash
#
cp asa091.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa091.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa091.cpp"
  exit
fi
#
mv asa091.o ~/libcpp/$ARCH/asa091.o
#
echo "Library installed as ~/libcpp/$ARCH/asa091.o"
