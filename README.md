# smartplus

How to install SMART+ :

You should have downloaded a zip file called smartplus-master.zip

1 - Unzip the file in the location where you want to use SMART+ and rename it "smartplus" (or as you choice).

2 - Go to the folder "smartplus"
        ** cd 'pathtothefile'/smartplus

3 - Execute the installation bash file :
        ** sh Install.sh

4 - After the compilation and linking is done, you can enjoy SMART+.

Many uses are available :

	a - Use the SMART+ solver.
	In a work folder, copy all necessary files from the folder 'exec'. 
	Edit configurations file to define your simulation (path.txt, material.dat...)
	Solve the problem executing "./solver" in a terminal (Linux and MacOS)

	b - Use SMART+ Umat for Abaqus. 
	Copy-past "umat_single.o" or "umat_singleT.o" from 'pathtothefile'/smartplus/build/bin to your abaqus work directory and use it like a classical Umat.
		Example : abaqus job=mymodel.inp user=umat_single.o

	c - Build your own projects using the smartplus lib so called "libsmartplus.so".
	For this, you need to specify the path of the smartplus dynamic library (with the -I option) and to link smartplus with your application
		Example : gcc myproject.cpp -I/'pathtothefile'/smartplus/lib -lsmartplus

For more explanation, go to http://www.lem3.fr/chemisky/smartplus/ and click on the  the menu 'examples' to find examples on how to use SMART+.

Have fun :)
