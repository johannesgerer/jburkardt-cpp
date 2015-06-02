#!/bin/bash
#
cp line_nco_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include line_nco_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_nco_rule.cpp"
  exit
fi
#
mv line_nco_rule.o ~/libcpp/$ARCH/line_nco_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/line_nco_rule.o"
