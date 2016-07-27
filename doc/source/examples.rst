Examples
========

Run a simulation
----------------

Here is a simple example to run your first simulation using SMART+.
Open the file Path.txt in the folder Control
You should find something that look like that

.. code-block:: none

    #Initial_temperature
    293.5
    #Number_of_blocks
    1

    #Block
    1
    #Loading_type
    1
    #Repeat
    1
    #Steps
    1

    #Mode
    1
    #Number_of_increments
    100
    #time
    30.
    #mechanical_state
    E 0.3 
    S 0 S 0
    S 0 S 0 S 0
    #temperature_state
    T 293.5

The first part of the file describe the initial temperature conditions, under the tag #Initial_temperature.

Just below, under the tag #Number_of_blocks, you define the number of blocks. Here we will start with a single block, so this value is set to 1.
The next part is to define the first block:
#Block defines the block number
#Loading_type defines if the loading is : 1 – linear; 2 – sinusoidal; 3 – tabular (from a file)
#Repeat is the number of time the block is repeated
#Steps is the number of steps of the block

The next part of the file defines the steps of the first block. It always starts with the mode of the step (#mode), which is:
1 – mechanical; 2 – thermomechanical

In this example we will consider that the step mode is mechanical. We therefore need to set up the following

#. The mode of the step (under #mode)
#. The number of the increments (under #number_of_increments)
#. The time of the step (under #time)
#. The mechanical loading stage at the end of the step (#mechanical_state)
#. The thermal loading stage at the end of the step (#temperature_state)

Set up a micro mechanical model
-------------------------------

The first thing you want to do when setting up a micro mechanical model is to define the microstructure. At a certain scale, you should inform the model about the phases, their volume fraction, geometry and their properties.

First, in the file data/material.dat, you need to enter the material properties corresponding to the micro mechanical model you selected:

For Mori-Tanaka and Self-Consistent: 4 material parameters (and a consequent number of state_variables)



#. props(0) : Number of phases
#. props(1) : File number that stores the microstructure properties
#. props(2) : Number of integration points in the 1 direction
#. props(3) : Number of integration points in the 2 direction

For Periodic layers: 2 material parameters (and a consequent number of state_variables)


#. props(0) : Number of phases
#. props(1) : File number that stores the microstructure properties

The file data/material.dat should look like this for a 2-phase material using a Mor-Tanaka model:

.. code-block:: none

	Material
	Name    MIMTN
	Number_of_material_parameters   4
	Number_of_internal_variables    10000

	#Thermal
	density 1.12
	c_p   1.64

	#Mechancial
	nphases 2
	file_number 0
	nItg1 20
	nItg2 20

The density and specific heat capacity c_p are utilized only if you want to solve a thermomechanical boundary-value problem.

The file number represents the number of the Nphases[i].dat file, where [i] is replaced by the number value. In this case we should fill the file Nphases0.dat, which looks like this:

.. code-block:: none

    Number  Coatingof  umat   c    phi_mat  theta_mat  psi_mat  a1  a2  a3  phi_geom  theta_geom  psi_geom  nprops  nstatev  props
    0       0          ELISO  0.8  0        0          0        1   1   1   0.        0.          0.        3       1        3000    0.4   1.E-5
    1       0          ELISO  0.2  0        0          0        1   1   1   0.        0.          0.        3       1        70000   0.4   1.E-5

Note that for Mori-Tanaka the first phase in the file should always be the matrix.
The characteristics of the phases are described below:

#. Number : The number of the phase
#. Coatingof : If the model is a coating of an other phase. 0 if the phase is not a coating
#. umat : Constitutive model considered
#. c : Volume fraction of the phase
#. phi_mat: First Euler angle corresponding to the material orientation
#. theta_mat: Second Euler angle corresponding to the material orientation
#. psi_mat: Third Euler angle corresponding to the material orientation
#. a1:
#. a2:
#. a3:
#. phi_geom: First Euler angle corresponding to the ellipsoid orientation
#. theta_geom: Second Euler angle corresponding to the ellipsoid orientation
#. psi_geom: Third Euler angle corresponding to the ellipsoid orientation
#. npros: Number of material properties
#. nstatev: Number of scalar internal variables
#. props: The list of material properties

For a wide majority of composites, the orientation of the material coincides with the orientation of the reinforcement (For instance transversely isotropic carbon fibers).
However, for metallic polycristals, the two materials systems have to be considered to separate the orientation of the lattice with the orientation of the ellipsoid that represent a grain.
This version of SMART+ currently does not support coated inclusions, but the files Nphase[i].dat is prepared so that you can easily add this to a custom micromechancial model.

Note that the Euler system reference utilised (3-1-3 for the most common) is defined in the parameter.hpp file. For instance this system is defined by default in the parameter.hpp:

.. code-block:: none

    #ifndef axis_psi
    #define axis_psi 3
    #endif

    #ifndef axis_theta
    #define axis_theta 1
    #endif

    #ifndef axis_phi
    #define axis_phi 3
    #endif

In the example here we are defining a 2-phase composite, with spherical reinforcements, considering two phases:

#. An epoxy matrix, 80% volume, with E=3000MPa and nu=0.4, and alpha=1.E-5
#. Aluminium reinforcements: 20% volume, with E=70000MPa and nu=0.3, and alpha=5.E-5

Once these files have been set up, you can run a simulation using the classical solver.
