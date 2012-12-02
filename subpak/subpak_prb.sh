#!/bin/bash
#
g++ -c -g -I/$HOME/include subpak_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subpak_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ subpak_prb.o /$HOME/libcpp/$ARCH/subpak.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subpak_prb.o."
  exit
fi
#
rm subpak_prb.o
#
mv a.out subpak_prb
./subpak_prb > subpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subpak_prb."
  exit
fi
rm subpak_prb
#
echo "Program output written to subpak_prb_output.txt"
