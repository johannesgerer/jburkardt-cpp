#!/bin/bash
#
g++ -c -g fsu_hammersley.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fsu_hammersley.cpp"
  exit
fi
rm compiler.txt
#
mv fsu_hammersley.o ~/libcpp/$ARCH/fsu_hammersley.o
#
echo "A new version of fsu_hammersley has been created."
