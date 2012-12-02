#!/bin/bash
#
g++ -c -g -I/$HOME/include chrpak_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chrpak_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ chrpak_prb.o /$HOME/libcpp/$ARCH/chrpak.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chrpak_prb.o."
  exit
fi
#
rm chrpak_prb.o
#
mv a.out chrpak_prb
./chrpak_prb > chrpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chrpak_prb."
  exit
fi
rm chrpak_prb
#
echo "Program output written to chrpak_prb_output.txt"
