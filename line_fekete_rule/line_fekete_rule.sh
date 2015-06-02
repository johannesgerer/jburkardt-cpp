#!/bin/bash
#
cp line_fekete_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include line_fekete_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_fekete_rule.cpp"
  exit
fi
#
mv line_fekete_rule.o ~/libcpp/$ARCH/line_fekete_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/line_fekete_rule.o"
