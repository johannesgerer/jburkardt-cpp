#!/bin/bash
#
g++ -c -I$HOME/include ppmb_to_bmp.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppmb_to_bmp.cpp"
  exit
fi
rm compiler.txt
#
g++ ppmb_to_bmp.o ~/libcpp/$ARCH/bmp_io.o ~/libcpp/$ARCH/ppmb_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ppmb_to_bmp.o"
  exit
fi
#
rm ppmb_to_bmp.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ppmb_to_bmp
#
echo "A new version of ppmb_to_bmp has been created."
