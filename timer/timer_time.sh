#!/bin/bash
#
g++ -c timer_time.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_time.cpp"
  exit
fi
rm compiler.txt
#
g++ timer_time.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading timer_time.o."
  exit
fi
#
rm timer_time.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/timer_time
#
echo "Executable installed as ~/bincpp/$ARCH/timer_time"
