# Lattice-mediated flexo contributions from numerical derivatives of the dynamical matrix

Here I will give instructions on how to obtain the LM flexoelectric coefficients via taking the second derivative of the dynamical matrix.

## DFT calculation:

- You must use the version of abinit 9 available here: https://github.com/Dreyer-Research-Group/abinit-grp-develop .  Note that in the below example, the variable `nogzero` is set, which is not available in the production version of abinit.

- You must calculate the Born effective charges, as well as the dynamical matrix for three specific q points with magnitude q=0.0, q=0.01, and q=-0.01 (in that order). The direction of q is your choice, but it should be along one of the axes in reduced coordinates, i.e., (0.0,0.01,0.0) or (0.0,0.0,0.01) or (0.01,0.0,0.0). 

- An example of the begining of the input file is give here (make sure to set `rfatpol` to the total number of atoms in your cell):

```
 ndtset 7

#Ground state calculation ##############################################
  kptopt1   1
  tolvrs1   1.0d-18
    iscf1   7
  prtden1   1
########################################################################

#Non-self consistent ground-state calculation##########################

#with q=(0.0 0.01 0.0)
    nqpt2   1
     qpt2   0.0 0.0 0.01
  getwfk2   1
  getden2   1
  kptopt2   3
  tolwfr2   1.0d-18
    iscf2   -2

#with q=(0.0 0.0 -0.01)
    nqpt3   1
     qpt3   0.0 0.0 -0.01
  getwfk3   1
  getden3   1
  kptopt3   3
  tolwfr3   1.0d-18
    iscf3   -2


#########################################################################

#Static Response Function calculation####################################

#Set 4 : Response function calculation of d/dk wave function

    iscf4   -3         # Need this non-self-consistent option for d/dk
  kptopt4   2          # Modify default to use time-reversal symmetry
  rfphon4   0          # Cancel default
  rfelfd4   2          # Calculate d/dk wave function only
  tolwfr4   1.0d-22    # Use wave function residual criterion instead
  getwfk4   1
 rfatpol4   1 4        # Treat displacements of all atoms
   rfdir4   1 1 1      # Do all directions (symmetry will be used)

#Set 5 : Response function calculation of Q=0 phonons and electric field pert.

  getddk5   4          # d/dk wave functions from last dataset
  kptopt5   2          # Modify default to use time-reversal symmetry
  rfelfd5   3          # Electric-field perturbation response only
  tolvrs5   1.0d-8
  rfphon5   1
 rfatpol5   1 4
   rfdir5   1 1 1
  getwfk5   1
 nogzero5   1
     qpt5   0.0 0.0 0.0


# q= 0 0.01 0
  rfphon6   1
 rfatpol6   1 4
   rfdir6   1 1 1
    nqpt6   1
     qpt6   0.0 0.0 0.01
  getwfk6   1
  getwfq6   2
  kptopt6   3
  tolvrs6   1.0d-8
    iscf6   7
 nogzero6 1

#q = 0 -0.01 0
  rfphon7   1
  rfatpol7   1 4
  rfdir7   1 1 1
  nqpt7  1
  qpt7  0.0 0.0 -0.01
  getwfk7   1
  getwfq7   3
  kptopt7   3
  tolvrs7   1.0d-8
    iscf7   7
  nogzero7  1

#######################################################################

```

## Generating dyn.dat and bec.dat:

- Once the above calculation has completed, run the bash script `2dphi_4.sh` available here: https://github.com/Dreyer-Research-Group/flexo_scripts . This will generate two files, dyn.dat and bec.dat from information in your `*.abo` file. `dyn.dat` has the dynamical matrix (not including the masses), as well as some information about the structure, etc. `bec.dat` just has a list of the Born effective charges.

## Taking the derivatives 

- Running `2dphi_4.sh` will also run the executable `twoddyn_pi.x` which is a compiled version of the fotran 90 code: `twoddyn_psdoinv.f90` also available here: https://github.com/Dreyer-Research-Group/flexo_scripts. You can compile that code with the command: `gfortran -lblas -llapack -o twoddyn_pi.x twoddyn_psdoinv.f90`. Make sure you have, e.g., the Intel stack of modules loaded so it can find the blas and lapack libraries.

- The output file will be called `flexoforce.dat` and will list the square brackets and round brackets for each sublattice [defined in Eq. 38 and 39 of Stengel PRB 88, 174106 (2013)], the type I flexo internal strain tensor N (Eq. 42), and type I LM flexo coefficient (last line of Eq. 49).

- Using Eqs. 59 and 61 of Stengel PRB 88, 174106 (2013), you can sum over atoms the round/square brackets to get the elastic constant to compare with a strain perturbation run.