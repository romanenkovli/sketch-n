#!/bin/bash
cd ../Include
cp ../Samples/PWR_KOEBERG2D/Include/koeberg2d_2x2.fh parameters.fh
cd ../Source
scons
cd ../Input
cp ../Samples/PWR_KOEBERG2D/Input/KOEBERG2D_2x2.DAT KOEBERG2D_2x2.DAT
cp ../Samples/PWR_KOEBERG2D/Input/KOEBERG2D.ref KOEBERG2D.ref
cp ../Samples/PWR_KOEBERG2D/Input/SKETCH.INI.KOEBERG2D_2x2.ST SKETCH.INI
cd ..
./sketch.exe
cd Shell