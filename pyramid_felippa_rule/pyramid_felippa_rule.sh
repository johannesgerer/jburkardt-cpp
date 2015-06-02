#!/bin/bash
#
cp pyramid_felippa_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include pyramid_felippa_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_felippa_rule.cpp"
  exit
fi
#
mv pyramid_felippa_rule.o ~/libcpp/$ARCH/pyramid_felippa_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/pyramid_felippa_rule.o"
