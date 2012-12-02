#!/bin/bash
#
cp cnf_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include cnf_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cnf_io.cpp"
  exit
fi
rm compiler.txt
#
mv cnf_io.o ~/libcpp/$ARCH/cnf_io.o
#
echo "Library installed as ~/libcpp/$ARCH/cnf_io.o"
