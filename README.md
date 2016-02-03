# smartplus

How to install SMART+ :

You should have downloaded a zip file called smartplus-master.zip.

1 - Unzip the file in the location where you want to use SMART+ and rename it "smartplus" (or as you wish).

2 - Go to the folder "smartplus"

        ** cd 'pathtothefile'/smartplus

3 - Execute the installation bash file :

        ** sh Install.sh

4 - After the compilation and linking is done, you can enjoy SMART+.

Many possibilities to use it are available :

A ) Use the SMART+ solver.

Copy all necessary files from the folder "exec" to your a work folder. 
Edit configuration files to define your simulation ("path.txt", "material.dat" and others if necessary ).
Solve the problem executing "./solver" in a terminal (Linux and MacOS).

B ) Use SMART+ Umat for Abaqus. 

Copy-paste "umat_single.o" or "umat_singleT.o" from 'pathtothefile'/smartplus/build/bin to your Abaqus work directory and use it like a classical Umat.
    
	Example : abaqus job=mymodel.inp user=umat_single.o

C ) Build your own projects using the SMART+ library called "libsmartplus.so".

For this, you need to specify the path of the SMART+ dynamic library (with the -I option) and to link it with your application.
    
	Example : gcc myproject.cpp -I/'pathtothefile'/smartplus/lib -lsmartplus

For more explanation, go to http://www.lem3.fr/chemisky/smartplus/ and click on the  the menu 'examples' to find examples on how to use SMART+.

Have fun :)
