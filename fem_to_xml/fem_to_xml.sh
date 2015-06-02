#!/bin/bash
#
g++ -c fem_to_xml.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_xml.cpp"
  exit
fi
#
g++ fem_to_xml.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_xml.o."
  exit
fi
#
rm fem_to_xml.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fem_to_xml
#
echo "Executable installed as ~/bincpp/$ARCH/fem_to_xml"
