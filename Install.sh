#!/bin/bash

echo "\nStart of SMART+ compilation.\n"

current_dir=$(pwd)

if [ ! -d "build" ]
then
	mkdir ${current_dir}/build
	echo "Folder created \n"
else
	echo "Build directory already exists"
	
	while true; do
		read -p "Do you want to erase old compilation files ? \n y for yes or n for no (recommended)" yn
		case $yn in
			[YyOo]* ) rm -r ${current_dir}/build/*; break;;
			[Nn]* ) break;;
			* ) echo "Please answer yes (y) or no (n).";;
		esac
	done
fi

echo ""
cd ${current_dir}/build
cmake ..
echo ""
make 

echo "\n---------------------------"

if [ -f ${current_dir}/build/CMakeFiles/umat.dir/software/umat_single.cpp.o ]
then 
	cp ${current_dir}/build/CMakeFiles/umat.dir/software/umat_single.cpp.o ${current_dir}/build/bin/umat_single.o
	echo "umat_single.o copied in ${current_dir}/build/bin"
fi

if [ -f ${current_dir}/build/CMakeFiles/umatT.dir/software/umat_singleT.cpp.o ]
then 
	cp ${current_dir}/build/CMakeFiles/umatT.dir/software/umat_singleT.cpp.o ${current_dir}/build/bin/umat_singleT.o
	echo "umat_singleT.o copied in ${current_dir}/build/bin"
fi

cp ${current_dir}/build/bin/solver ${current_dir}/exec
echo "Solver copied in ${current_dir}/exec"
echo "libsmartplus.so available in ${current_dir}/lib"
echo "---------------------------"
echo "Compilation of SMART+ done. \n"

