#!/bin/bash

cd v2/src
if [ ! -f "kshell.exe" -o ! -f "transit.exe" ]; then
	make clean
	make -j4
fi
cp kshell.exe ../../Ne22_v2
cp transit.exe ../../Ne22_v2
cp kshell.exe ../../V50_v2
cp transit.exe ../../V50_v2
cd ../..

cd v4/src
if [ ! -f "kshell.exe" -o ! -f "transit.exe" ]; then
	make clean
	make -j4
fi
cp kshell.exe ../../Ne22_v4
cp transit.exe ../../Ne22_v4
cp kshell.exe ../../V50_v4
cp transit.exe ../../V50_v4
cd ../..

cd Ne22_v2
rm *.wav
rm log_*
rm summary_*
bash Ne22_usda.sh
cd ..
cd Ne22_v4
rm *.wav
rm log_*
rm summary_*
bash Ne22_usda.sh
cd ..
cd V50_v2
rm *.wav
rm log_*
rm summary_*
bash V50_gxpf1a.sh
cd ..
cd V50_v4
rm *.wav
rm log_*
rm summary_*
bash V50_gxpf1a.sh
cd ..
python compare.py
