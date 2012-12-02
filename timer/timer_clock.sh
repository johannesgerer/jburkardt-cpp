#!/bin/bash
#
g++ -c timer_clock.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_clock.cpp"
  exit
fi
rm compiler.txt
#
g++ timer_clock.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading timer_clock.o."
  exit
fi
#
rm timer_clock.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/timer_clock
#
echo "Executable installed as ~/bincpp/$ARCH/timer_clock"
