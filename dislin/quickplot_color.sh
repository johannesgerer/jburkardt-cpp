#!/bin/bash
#
g++ -c -I/usr/local/dislin quickplot_color.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_color.cpp"
  exit
fi
rm compiler.txt
#
g++ quickplot_color.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_color.o."
  exit
fi
#
rm quickplot_color.o
#
mv a.out quickplot_color
./quickplot_color > quickplot_color_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_color."
  exit
fi
rm quickplot_color
#
echo "Program output written to quickplot_color_output.txt"
