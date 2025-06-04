cd ../Include

cp ../Samples/VVER_IAEA2D_NOR/Include/vver_iaea2d_nor.fh parameters.fh
cd ../Source
  scons
cd ../Input
cp ../Samples/VVER_IAEA2D_NOR/Input/VVER_IAEA2D_NOR_125.DAT VVER_IAEA2D_NOR_125.DAT
cp ../Samples/VVER_IAEA2D_NOR/Input/VVER_IAEA2D_NOR_125.ref VVER_IAEA2D_NOR_125.ref
cp ../Samples/VVER_IAEA2D_NOR/Input/SKETCH.INI.VVER_IAEA2D_NOR_125.ST SKETCH.INI
cd ..
./sketch.exe