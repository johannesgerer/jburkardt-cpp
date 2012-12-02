#!/bin/bash
#
g++ -c -I$HOME/include ann_to_fig.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ann_to_fig.cpp"
  exit
fi
rm compiler.txt
#
g++ ann_to_fig.o -L$HOME/libcpp/$ARCH -lann -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ann_to_fig.o + libann"
  exit
fi
#
rm ann_to_fig.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ann_to_fig
#
echo "A new version of ann_to_fig has been created."
