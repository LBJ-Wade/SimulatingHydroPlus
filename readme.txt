Causal viscous hydro code for central heavy-ion collisions
by Paul Romatschke,	v0.0: 1/2007
			v0.1: 2/2007
hacked by Gregory Ridgway, v1.0 11/2018
Hydro Code

-------------------------------------------------------
File			Explanation

VH1+1.cpp		Main Code source file
crit.cpp    Critical EoS/Mode ammendments to Paul's code
GJE.cpp			Gauss-Jordan elimination code for matrices
paramreader.cpp		Reads actual run parameters from params.txt file
params.txt		Hydro Code parameters (documented)
diags.cpp		Diagnostic routines
vh			this is the executeable code
-----------------------------------------------------
Compilation:

On a standard Linux machine, you should be able to type 'make'
and get an executeable file vh. Don't know about other OS, but
who wants to use other OS anyways?

-------------------------------------------------------

Running the code

Having compiled the code, it should run when you type in

./vh

To change lattice size/hydrodynamic parameters, see params.txt

However, you should create the directories

../data

and 

../data/snapshot

in order for the code to be able to output something useful.
See below

----------------------------------------------------------

Data output will created a file

freezeout.dat 

in ../data/

which contains the information on the freeze-out surface. It has
the following format:

#radius_in_fm	#time_in_fm/c	#u^t  #u^r  #Pi^r_r/(e+p) #Pi^eta_eta/(e+p) #p 
and 4 more numbers that probably are not important to you (if you're 
interested, see outputMeasurements routine to find out what they are

in snapshot/
you will find snapshots of the temperature profile, energy density profile,...
at given intervals in proper time

