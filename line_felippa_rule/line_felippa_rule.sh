#!/bin/bash
#
cp line_felippa_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include line_felippa_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling line_felippa_rule.cpp"
  exit
fi
#
mv line_felippa_rule.o ~/libcpp/$ARCH/line_felippa_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/line_felippa_rule.o"
