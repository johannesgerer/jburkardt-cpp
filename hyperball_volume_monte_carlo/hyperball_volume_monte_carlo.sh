#!/bin/bash
#
g++ -c hyperball_volume_monte_carlo.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hyperball_volume_monte_carlo.cpp"
  exit
fi
#
g++ hyperball_volume_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hyperball_volume_monte_carlo.o."
  exit
fi
#
rm hyperball_volume_monte_carlo.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hyperball_volume_monte_carlo
#
echo "Executable installed as ~/bincpp/$ARCH/hyperball_volume_monte_carlo"
