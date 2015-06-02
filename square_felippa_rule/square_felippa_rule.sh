#!/bin/bash
#
cp square_felippa_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include square_felippa_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling square_felippa_rule.cpp"
  exit
fi
#
mv square_felippa_rule.o ~/libcpp/$ARCH/square_felippa_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/square_felippa_rule.o"
