#!/bin/bash
#
g++ -c -g discrete_pdf_sample.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling discrete_pdf_sample.cpp"
  exit
fi
rm compiler.txt
#
g++ discrete_pdf_sample.o
if [ $? -ne 0 ]; then
  echo "Errors while loading discrete_pdf_sample.o"
  exit
fi
rm discrete_pdf_sample.o
#
mv a.out ~/bincpp/$ARCH/discrete_pdf_sample
#
echo "Executable installed as ~/bincpp/$ARCH/discrete_pdf_sample"
