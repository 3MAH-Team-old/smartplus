#!/bin/bash

echo "\nStart the compilation of SMART+\n"

current_dir=$(pwd)

if [ ! -d "build" ]
then
	mkdir ${current_dir}/build
	echo "Folder created"
else
	echo "Build directory already exists"
	
	while true; do
		read -p "Do you want to erase old comilation files ? " yn
		case $yn in
			[YyOo]* ) rm -r ${current_dir}/build/*; break;;
			[Nn]* ) break;;
			* ) echo "Please answer yes or no.";;
		esac
	done
fi

cd ${current_dir}/build
cmake ..
make 

echo "---------------------------"

if [ -d "CMakeFiles/umat.dir" ]
then 
	cp ${current_dir}/build/CMakeFiles/umat.dir/software/umat_single.cpp.o ${current_dir}/build/bin/umat_single.o
	echo "umat_single.o copied in ${current_dir}/build/bin"
fi

if [ -d "CMakeFiles/umatT.dir" ]
then 
	cp ${current_dir}/build/CMakeFiles/umatT.dir/software/umat_singleT.cpp.o ${current_dir}/build/bin/umat_singleT.o
	echo "umat_singleT.o copied in ${current_dir}/build/bin"
fi

cp ${current_dir}/build/bin/solver ${current_dir}/exec
echo "Solver copied in ${current_dir}/exec"
echo "libsmartplus.so available in ${current_dir}/lib"
echo "---------------------------"
echo "Compilation of SMART+ done.\n"

