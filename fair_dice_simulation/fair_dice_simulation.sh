#!/bin/bash
#
g++ -c -g fair_dice_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fair_dice_simulation.cpp"
  exit
fi
rm compiler.txt
#
g++ fair_dice_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fair_dice_simulation.o"
  exit
fi
rm fair_dice_simulation.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fair_dice_simulation
#
echo "Executable installed as ~/bincpp/$ARCH/fair_dice_simulation"
