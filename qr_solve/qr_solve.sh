#!/bin/bash
#
cp qr_solve.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include qr_solve.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qr_solve.cpp"
  exit
fi
rm compiler.txt
#
mv qr_solve.o ~/libcpp/$ARCH/qr_solve.o
#
echo "Library installed as ~/libcpp/$ARCH/qr_solve.o"
