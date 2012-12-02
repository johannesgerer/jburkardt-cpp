#!/bin/bash
#
cp test_opt_con.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include test_opt_con.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_opt_con.cpp."
  exit
fi
rm compiler.txt
#
mv test_opt_con.o ~/libcpp/$ARCH/test_opt_con.o
#
echo "Library installed as ~/libcpp/$ARCH/test_opt_con.o"
