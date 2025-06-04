#! /bin/bash
cd ../Include
cp ../Samples/PWR_NEACRP_A1/Include/pwr_neacrp_a1.fh parameters.fh
cp /b parameters.fh+,, 
cd ../Source
scons -j8
# make  -f makefile_ifort
cd ../Input
cp ../Samples/PWR_NEACRP_A1/Input/PWR_NEACRP_A1.DAT PWR_NEACRP_A1.DAT
cp ../Samples/PWR_NEACRP_A1/Input/SKETCH.INI.PWR_A1.ST SKETCH.INI
read
cd ..
./sketch.exe
cd Shell