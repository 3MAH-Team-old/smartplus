# smartplus

How to install SMART+ :

You should have downloaded a zip file called smartplus-master.zip

1 - Unzip the file in the location where you want to use SMART+ and rename it "smartplus"

2 - Go to the folder "smartplus"
        ** cd 'pathtothefile'/smartplus

3 - create a build folder 
        ** mkdir build

4 - Use Cmake to configure and generate the Makefile
        ** cmake ..

5 - Build the libraries and the executable files
        ** make

6 - After the compilation and linking is done, go to /build/bin/Release

    Copy-paste the solver exec to a folder where you want to run the solver executable file
    for example in a folder exec
    Copy the content of the folder "data" in the same folder.

7 - If you want to build projects with the smartplus lib, you can copy paste the library in the build/lib/Release to the folder /usr/lib, or to the folder /lib of your project 
    In the last case you may need to specify the path of the smartplus dynamic library (with the -I option)

8 - To link smartplus with your application compile it with using -lsmartplus
    To use the solver you can go to the 'exec' folder and use the command ./solver in a terminal (Linux and MacOS)

9 - Go to http://www.lem3.fr/chemisky/smartplus/ . Click on the  the menu 'examples' to find examples on how to use SMART+

Have fun :)
