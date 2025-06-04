#!/bin/bash
cd ../Include
cp ../Samples/PWR_KOEBERG2D/Include/koeberg2d.fh parameters.fh
cd ../Source
scons
cd ../Input
cp ../Samples/PWR_KOEBERG2D/Input/KOEBERG2D.DAT KOEBERG2D.DAT
cp ../Samples/PWR_KOEBERG2D/Input/KOEBERG2D.ref KOEBERG2D.ref
cp ../Samples/PWR_KOEBERG2D/Input/SKETCH.INI.KOEBERG2D.ST SKETCH.INI
cd ..
./sketch.exe
cd Shell