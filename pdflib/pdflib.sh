#!/bin/bash
#
cp pdflib.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include pdflib.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pdflib.cpp"
  exit
fi
rm compiler.txt
#
mv pdflib.o ~/libcpp/$ARCH/pdflib.o
#
echo "Library installed as ~/libcpp/$ARCH/pdflib.o"
