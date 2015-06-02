#!/bin/bash
#
g++ -c -g discrete_pdf_sample_2d.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling discrete_pdf_sample_2d.cpp"
  exit
fi
rm compiler.txt
#
g++ discrete_pdf_sample_2d.o
if [ $? -ne 0 ]; then
  echo "Errors while loading discrete_pdf_sample_2d.o"
  exit
fi
rm discrete_pdf_sample_2d.o
#
mv a.out ~/bincpp/$ARCH/discrete_pdf_sample_2d
#
echo "Executable installed as ~/bincpp/$ARCH/discrete_pdf_sample_2d"
