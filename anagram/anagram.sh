#!/bin/bash
#
g++ -c -g anagram.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling anagram.cpp."
  exit
fi
rm compiler.txt
#
g++ anagram.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading anagram.o."
  exit
fi
#
rm anagram.o
mv a.out ~/bincpp/$ARCH/anagram
#
echo "Executable installed as ~/bincpp/$ARCH/anagram"
