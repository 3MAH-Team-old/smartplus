#!/bin/bash
echo "\n----------------------------"
echo "Start of SMART+ compilation."
echo "---------------------------\n"

#Find th current directory
current_dir=$(pwd)

#Test if build exist and if it's necessary to erase it
if [ ! -d "build" ]
then
	mkdir ${current_dir}/build
	echo "Folder created.\n"
else
	echo "Build directory already exists."
	
	while true; do
		read -p "Do you want to erase old compilation files (Recommended : No) ? " yn
		case $yn in
			[YyOo]* ) rm -r ${current_dir}/build/*; break;;
			[Nn]* ) break;;
			* ) echo "Please answer yes (y) or no (n).";;
		esac
	done
fi

#Build SMART+
echo ""
cd ${current_dir}/build
cmake ..
echo ""
make 


# Copy all important files
if [ $? -eq 0 ]
then
	echo "\n---------------------------"
	
	#Copy of umat_single.o
	if [ -f ${current_dir}/build/CMakeFiles/umat.dir/software/umat_single.cpp.o ]
	then 
		cp ${current_dir}/build/CMakeFiles/umat.dir/software/umat_single.cpp.o ${current_dir}/build/bin/umat_single.o
		echo "umat_single.o copied in ${current_dir}/build/bin"
	fi
	
	#Copy of umat_singleT.o
	if [ -f ${current_dir}/build/CMakeFiles/umatT.dir/software/umat_singleT.cpp.o ]
	then 
		cp ${current_dir}/build/CMakeFiles/umatT.dir/software/umat_singleT.cpp.o ${current_dir}/build/bin/umat_singleT.o
		echo "umat_singleT.o copied in ${current_dir}/build/bin"
	fi
	
	#if debug existsn, copy of solver from Debug
	if [ -f ${current_dir}/build/bin/Debug/solver ]
	then
		cp ${current_dir}/build/bin/Debug/solver ${current_dir}/build/bin
	fi

    #if debug existsn, copy of solver from Debug
    if [ -f ${current_dir}/build/bin/Debug/identification ]
    then
    cp ${current_dir}/build/bin/Debug/identification ${current_dir}/build/bin
    fi

	cp ${current_dir}/build/bin/solver ${current_dir}/exec
	cp ${current_dir}/build/bin/identification ${current_dir}/exec
	echo "Solver copied in ${current_dir}/exec"
	echo "libsmartplus.so available in ${current_dir}/lib"
	echo "---------------------------"
	echo "SMART+ compilation done.\n"
else
	echo "\n---------------------------"
	echo "SMART+ compilation failed.\n"
fi
