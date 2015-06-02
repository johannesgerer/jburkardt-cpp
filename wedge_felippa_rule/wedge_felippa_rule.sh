#!/bin/bash
#
cp wedge_felippa_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include wedge_felippa_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_felippa_rule.cpp"
  exit
fi
#
mv wedge_felippa_rule.o ~/libcpp/$ARCH/wedge_felippa_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/wedge_felippa_rule.o"
