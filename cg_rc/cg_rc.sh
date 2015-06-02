#!/bin/bash
#
cp cg_rc.hpp /$HOME/include
#
g++ -c -g -I/$HOME/include cg_rc.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_rc.cpp"
  exit
fi
rm compiler.txt
#
mv cg_rc.o ~/libcpp/$ARCH/cg_rc.o
#
echo "Library installed as ~/libcpp/$ARCH/cg_rc.o"
