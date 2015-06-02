#!/bin/bash
#
cp high_card_simulation.hpp /$HOME/include
#
g++ -c -I/$HOME/include high_card_simulation.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling high_card_simulation.cpp"
  exit
fi
#
mv high_card_simulation.o ~/libcpp/$ARCH/high_card_simulation.o
#
echo "Library installed as ~/libcpp/$ARCH/high_card_simulation.o"
