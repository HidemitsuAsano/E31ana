#!/bin/sh

if [ -z "$G4SIMHOME" ] ; then
    G4SIMHOME=../geant/knucl4.10/
fi

if [ $# -eq 1 ] ; then
    G4SIMHOME=$1
fi

echo "KnuclRootData & ComCrossSection install from $G4SIMHOME "

cp $G4SIMHOME/src/ComCrossSectionTable.cc src/ComCrossSectionTable.cpp
cp $G4SIMHOME/include/ComCrossSectionTable.hh src/ComCrossSectionTable.h
cp $G4SIMHOME/src/KnuclRootData.cc src/KnuclRootData.cpp
cp $G4SIMHOME/include/KnuclRootData.h src/KnuclRootData.h

sed -i -e "s/ComCrossSectionTable.hh/ComCrossSectionTable.h/g" src/ComCrossSectionTable.cpp
sed -i -e "s/ComCrossSectionTable.hh/ComCrossSectionTable.h/g" src/KnuclRootData.h