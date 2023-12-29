This is the main code to generate the trajectories.
Before compiling the code, please change the following appropriately.

NPart = Number of cells in the simulation
Delta_t = time_step for the simulations
PHI = packing fraction  
RMEAN = mean value of radii <R>
DelR = dispersion in R
Eval = value of E_i
eta = value of eta
fad = value of fad
EquiLength = length of equilibration run
ProdLength = length of production run


To compile the code, use mpicc -w -O3 mainCode.c -lm
To run mpirun -np $N1 ./a.out $N2. $N1 is the number of independent runs. $N2 is the input argument to create folders for these independent runs. If $N2 is 0 and $N1 is 8 it 
will create folders 000 to 007 for 8 independent runs
