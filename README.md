# smartplus

How to install the SMART+ :

copy and paste the content of the library where you want to install SMART+
create a folder /build. So in Linux/Unix systems you can type
Use Cmake to generate source in the folder /build

mkdir build
cd build
ccmake ../Hello 

Open a terminal
Navigate to the folder /build

type the command:
make

After all the compilation of the files, go to /build/bin/Release and copy-paste the solver exec to a folder where you want to run the solver executable file
You have to copy the content of the folder "data" in the same folder.

Copy the content of build/lib/Release to a folder /lib of your project (Or you can manually paste it with appropriate right to /usr/lib you you want to use it from anywhere)

Have fun :)
