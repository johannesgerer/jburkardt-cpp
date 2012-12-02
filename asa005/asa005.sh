#!/bin/bash
#
cp asa005.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include asa005.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa005.cpp"
  exit
fi
rm compiler.txt
#
mv asa005.o ~/libcpp/$ARCH/asa005.o
#
echo "Library installed as ~/libcpp/$ARCH/asa005.o"
