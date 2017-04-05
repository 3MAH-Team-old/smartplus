#!/bin/bash
echo "\n----------------------------"
echo "Start of SMART+ compilation."
echo "---------------------------\n"

blue=`tput setaf 4`
red=`tput setaf 1`
reset=`tput sgr0`

#Find the current directory
current_dir=$(pwd)

#Get the file of the bash file (possibly necessary)
bashFileName=${0##*/}

# Detect the platform (similar to $OSTYPE)
OS="`uname`"
case $OS in
'Linux')
OS='Linux'
alias finder='find'
alias nprocs='nproc'
;;
'Darwin')
OS='Mac'
alias finder='gfind'
alias nprocs='sysctl -n hw.ncpu'
;;
'WindowsNT')
OS='Windows'
echo 'No equivalent function known to be equivalent to "find or gfind" : Please search for one and replace this "echo" command in ${bashFileName} and replace it by an "alias"'
;;
*) ;;
esac

#Number of procs used to compil (+1 because if nproc=1 => 1/2=0)
nproc_used=$(( ($(nprocs)+1)/2 ))

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

#Ask for installation of the SMART+ library
while true; do
	read -p "Do you want to install SMART+ library (necassary to use libsmartplus.so and simmit) ? " yn
	case $yn in
		[YyOo]* ) Install_check='OK'; break;;
		[Nn]* ) Install_check='NO'; break;;
		* ) echo "Please answer yes (y) or no (n).";;
	esac
done

#Build SMART+
echo ""
cd ${current_dir}/build
cmake ..
echo ""
make -j${nproc_used}
echo ""
if [ "${Install_check}" = "OK" ]
then
	make install
fi
Install_OK=$?

if [ $Install_OK -eq 0 ]
then
	make test

	Test_OK=$?


	#Create the list of the file to copy after compilation
	executableToCopy="solver identification L_eff Elastic_props ODF"
	objectToCopy="umat_single umat_singleT"

	# Copy all important files (+ final message)
	if [ $Test_OK -eq 0 ]
	then
		echo "\n---------------------------"
		
		#Treatement of object files
		for object in ${objectToCopy}
		do
			#Copy of the "object".o from build/CMakeFiles/umat.dir/software to build/bin
			if [ -f ${current_dir}/build/CMakeFiles/umat.dir/software/${object}.cpp.o ]
			then 
				cp ${current_dir}/build/CMakeFiles/umat.dir/software/${object}.cpp.o ${current_dir}/build/bin/${object}.o
				echo "${blue}${object}.o${reset} copied in ${blue}${current_dir}/build/bin${reset}"
			fi
		done
		
		#Treatement of executable files
		for file in ${executableToCopy}
		do
			#if debug exists, copy of the file from build/bin/Debug to build/bin
			if [ -f ${current_dir}/build/bin/Debug/${file} ]
			then
				cp ${current_dir}/build/bin/Debug/${file} ${current_dir}/build/bin
			fi
			
			#Copy the file from build/bin to exec
			cp ${current_dir}/build/bin/${file} ${current_dir}/exec
			echo "${blue}${file}${reset} copied in ${blue}${current_dir}/exec${reset}"
		done
		
		if [ "${Install_check}" = "OK" ]
		then
			echo "${blue}libsmartplus.so${reset} installed in ${blue}${current_dir}/lib${reset}"
		else
			echo "${blue}libsmartplus.so${reset} not installed."
		fi
		
		echo "---------------------------"
		echo "SMART+ compilation done.\n"
	else
		echo "\n---------------------------"
		echo "${red} SMART+ tests failed.\n"
	fi

else

	echo "\n---------------------------"
	echo "${red} SMART+ compilation failed"

fi
