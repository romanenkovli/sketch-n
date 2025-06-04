#! /bin/bash
cd ../Input
cp ../Samples/PWR_NEACRP_A1/Input/PostProc.INI.PWR_A1_KIN PostProc.INI
cp ../Samples/PWR_NEACRP_A1/Input/PP_PWR_A1_KIN.dat PP_PWR_A1_KIN.dat
cd ..
./postproc
cd Shell
