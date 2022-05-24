# Calculating the clamped-ion flexo coefficients using the phonon technique

# Preparing the DFPT calculations:

- The first step is to create a file called `Flexo.dat` which has information about the calculation. An example is given here:

```
# BN muXX33

name iBN.in

nqpt 3

qpts
0.0  0.0 0.0
0.01 0.0 0.0
0.02 0.0 0.0

rfatpol 1 4

rfdir 1 1 1

split 1 1 #qpoints and atoms

#Common input variables################################################

#Definition of the k-point grids
 nshiftk 1
 ngkpt   8 8 6
 shiftk 0.0 0.0 0.5
 kptopt 3

 istwfk 384*1

#Definition of the unit cell
acell      4.6921218563E+00  4.6921218563E+00  1.2261816468E+01
rprim      1.0000000000E+00  0.0000000000E+00  0.0000000000E+00
          -5.0000000000E-01  8.6602540378E-01  0.0000000000E+00
           0.0000000000E+00  0.0000000000E+00  1.0000000000E+00

#Definition of the atom types
ntypat 2
znucl 5 7

#Definition of the atoms
natom 4          
typat 1 1 2 2   
xred
0.333333333333 0.666666666666 0.25000000000 #B
0.666666666666 0.333333333333 0.75000000000 #B
0.333333333333 0.666666666666 0.75000000000 #N
0.666666666666 0.333333333333 0.25000000000 #N

#Definition of the planewave basis set
ecut 60         # Maximal kinetic energy cut-off, in Hartree

#Definition of the SCF procedure
nstep 100          # Maximal number of SCF cycles
ixc 7            

nband 8

# add to conserve old < 6.7.2 behavior for calculating forces at each SCF step
 optforces 1

diemac 13

pseudos "B_LDA_noNLCC.psp8, N_LDA_noNLCC.psp8"


```

- Under the "Common input variables" line, you should replace the inputs with those from your calculations. 

- The `name` flag gives the name of the input file that will be created, `nqpt` is the number of $q$ points, `qpts` lists the $q$ points, `rfatpol` and `rfdir` have the same meaning as in abinit, and `split` governs how you want to split the work into multiple separate calculations (first number splits by $q$ points, second number by atoms, 1 indicates no splitting). NOTE: This is NOT an input file for abinit, it will be run with a script, see next bullet.

- You should have just `Flexo.dat` and your pseudopotentials in the folder. 

- Run the script `flexo_abin_gen_ab9.sh` in that folder (available: https://github.com/Dreyer-Research-Group/flexo_scripts). This will generate one or more folders with names like `q1_at1`, `q1_at2`, `q2_at1`, etc. Each folder will have an input file for an abinit DFPT calculation for the given $q$ point(s) and atom(s). If you set `split 1 1` then only one folder wil be produced to run all of the necessary perturbations; for many atoms, this can be quite long, which is why I added the functionality to automatically split them up.

## Run the abinit DFPT calculation

- No go into each `q*_at*` folder and run abinit. Again, you will have to use the version available here: https://github.com/Dreyer-Research-Group/abinit-grp-develop

## Collect the results

- Once the calculations have finished, go to the home folder (i.e., the one where you originally ran `flexo_abin_gen_ab9.sh` to get the `q*_at*` folders). Run the script `dist_col_ab9.sh` which is available here: https://github.com/Dreyer-Research-Group/flexo_scripts . This will generate a file called `totg_col.dat` containing all of the information about the flexo coefficients from your `q*_at*` folders. It will also create a file called `sum_totg_col.dat` which is sumed over atoms.

## Take the q derivatives for CI flexo and piezo coefficients

- Compile the fortran program `flexo_derivatives.f90` using, e.g., the command `gfortran -o flexo_derivatives.x flexo_derivatives.f90`. 

- When you run the program in the same folder as the above step, it will first prompt you whether or not wyou want your piezo coefficients in units of C/m or e/Bohr (flexo is always in C/m^2). Then it will generate three files, `Z_akb.dat`: Contains BECs, `derivatives_flexo_ab9.dat`: Contains CI sublattice-resolved flexo coefficients, `derivatives_piezo_2.dat`: Contains CI sublattice-resolved piezoelectric coefficients.

- NOTE: The above steps will only work for the longitudinal coefficients. Also you must have 3 $q$ points, evenly spaced, in order, with 0 as the first one. For transverse, you need to also calculate the CRG contributions [see PRB 98, 075153 (2018)]. I have not made a script to do this specifically, but will do so soon. Same holds for type I shear coefficients, and in addition, you have to use a different fortran program to take the derivatives, `flexo_shear_deriv_ab9.f90`.
