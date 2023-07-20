#!/bin/bash

if [ -z "$1" ] || [ "$1" != "long" -a "$1" != "short" ]; then
    echo "Invalid or no argument supplied. Please enter either 'long' or 'short'."
    exit 1
fi

if [ "$1" == "long" ]; then
	cd v2/src
	make clean
	make -j8
	cp kshell.exe ../../Ne22_v2/
	cp transit.exe ../../Ne22_v2/
	cp kshell.exe ../../V50_v2/
	cp transit.exe ../../V50_v2/
	cd ../..

	cd v4/src
	make clean
	make -j8
	cp kshell.exe ../../Ne22_v4/
	cp transit.exe ../../Ne22_v4/
	cp kshell.exe ../../V50_v4/
	cp transit.exe ../../V50_v4/
	cd ../..
fi

cd ../../src
make clean
make -j8
cp kshell.exe ../test/v2_to_v4/Ne22_dev/
cp transit.exe ../test/v2_to_v4/Ne22_dev/
cp kshell.exe ../test/v2_to_v4/V50_dev/
cp transit.exe ../test/v2_to_v4/V50_dev/
cp kshell.exe ../test/v2_to_v4/V51_dev/
cp transit.exe ../test/v2_to_v4/V51_dev/
cd ../test/v2_to_v4

if [ "$1" == "long" ]; then
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
fi

cd Ne22_dev
rm *.wav
rm log_*
rm summary_*
bash Ne22_usda.sh
cd ..

cd V50_dev
rm *.wav
rm log_*
rm summary_*
bash V50_gxpf1a.sh
cd ..

cd V51_dev
rm *.wav
rm log_*
rm summary_*
bash V51_gxpf1a.sh
cd ..

python compare.py "$1"
