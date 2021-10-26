#!/bin/bash

file=write_MuE_MCevents
fileC=${file}.C
fileEXE=${file}.exe

ROOTLIBS="`root-config --cflags --libs`"

ROOTINCDIR=`root-config --incdir`

echo "Compiling"

rootcling -f MuEtreeDict.C -c MuEtree.h MuEtreeLinkDef.h

g++ -I${ROOTINCDIR} ${fileC} MuEtreeDict.C ${ROOTLIBS} -o ${fileEXE}
