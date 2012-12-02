#!/bin/bash
#
g++ -c -I/usr/local/dislin quickplot_contour.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quickplot_contour.cpp"
  exit
fi
rm compiler.txt
#
g++ quickplot_contour.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quickplot_contour.o."
  exit
fi
#
rm quickplot_contour.o
#
mv a.out quickplot_contour
./quickplot_contour > quickplot_contour_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quickplot_contour."
  exit
fi
rm quickplot_contour
#
echo "Program output written to quickplot_contour_output.txt"
