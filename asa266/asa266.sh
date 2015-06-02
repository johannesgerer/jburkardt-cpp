#!/bin/bash
#
cp asa266.hpp /$HOME/include
#
g++ -c -I /$HOME/include asa266.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa266.cpp"
  exit
fi
#
mv asa266.o ~/libcpp/$ARCH/asa266.o
#
echo "Library installed as ~/libcpp/$ARCH/asa266.o"
