#!/bin/bash
#
g++ -c -I/usr/local/dislin quickplot_bar.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_bar.cpp"
  exit
fi
rm compiler.txt
#
g++ quickplot_bar.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_bar.o."
  exit
fi
#
rm quickplot_bar.o
#
mv a.out quickplot_bar
./quickplot_bar > quickplot_bar_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_bar."
  exit
fi
rm quickplot_bar
#
echo "Program output written to quickplot_bar_output.txt"
