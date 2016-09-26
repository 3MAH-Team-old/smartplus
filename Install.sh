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
	echo "Build folder created.\n"
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
make -j4
echo ""
make test
echo ""
make install && Install_check='OK'

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
	cp ${current_dir}/build/bin/solver ${current_dir}/exec
	echo "solver copied in ${current_dir}/exec"

	#if debug existsn, copy of solver from Debug
	if [ -f ${current_dir}/build/bin/Debug/identification ]
	then
		cp ${current_dir}/build/bin/Debug/identification ${current_dir}/build/bin
	fi
	cp ${current_dir}/build/bin/identification ${current_dir}/exec
	echo "identification copied in ${current_dir}/exec"

	#if debug existsn, copy of solver from Debug
	if [ -f ${current_dir}/build/bin/Debug/L_eff ]
	then
		cp ${current_dir}/build/bin/Debug/L_eff ${current_dir}/build/bin
	fi
	cp ${current_dir}/build/bin/L_eff ${current_dir}/exec
	echo "L_eff copied in ${current_dir}/exec"

	#if debug existsn, copy of solver from Debug
	if [ -f ${current_dir}/build/bin/Debug/ODF ]
	then
		cp ${current_dir}/build/bin/Debug/ODF ${current_dir}/build/bin
	fi
	cp ${current_dir}/build/bin/ODF ${current_dir}/exec
	echo "ODF copied in ${current_dir}/exec"
	
	if [ "${Install_check}" = "OK" ]
	then
		echo "${blue}libsmartplus.so${reset} installed in ${blue}${current_dir}/lib${reset}"
	else
		echo "${blue}libsmartplus.so${reset} not installed. Uncomment 'make install' line in ${0##*/} file if you want to."
	fi
	
	echo "---------------------------"
	echo "SMART+ compilation done.\n"
else
	echo "\n---------------------------"
	echo "SMART+ compilation failed.\n"
fi
