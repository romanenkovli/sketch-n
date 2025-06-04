#! /bin/bash
cd ../Input
cp ../Samples/PWR_NEACRP_A1/Input/PostProc.ini.PWR_A1_ST PostProc.INI
cp ../Samples/PWR_NEACRP_A1/Input/PP_PWR_A1_ST.dat PP_PWR_A1_ST.dat
cd ..
./postproc
cd Shell
