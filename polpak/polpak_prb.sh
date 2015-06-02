#!/bin/bash
#
g++ -c -I/$HOME/include polpak_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling polpak_prb.cpp"
  exit
fi
#
g++ polpak_prb.o /$HOME/libcpp/$ARCH/polpak.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading polpak_prb.o."
  exit
fi
#
rm polpak_prb.o
#
mv a.out polpak_prb
./polpak_prb > polpak_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running polpak_prb."
  exit
fi
rm polpak_prb
#
echo "Program output written to polpak_prb_output.txt"
