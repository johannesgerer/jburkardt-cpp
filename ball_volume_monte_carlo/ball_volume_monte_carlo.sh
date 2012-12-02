#!/bin/bash
#
g++ -c -g ball_volume_monte_carlo.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_volume_monte_carlo.cpp"
  exit
fi
rm compiler.txt
#
g++ ball_volume_monte_carlo.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ball_volume_monte_carlo.o."
  exit
fi
#
rm ball_volume_monte_carlo.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ball_volume_monte_carlo
#
echo "Executable installed as ~/bincpp/$ARCH/ball_volume_monte_carlo"
