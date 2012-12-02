#!/bin/bash
#
g++ -c mgs.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling mgs.cpp"
  exit
fi
#
echo "The mgs file was compiled."
