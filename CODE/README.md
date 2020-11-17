This is Giorgio Amati's BGK CODE: it is his personal playground..... :)

├── RUN: working directory

├── SRC: source code

The gpu version needs pgi compiler. 
This version is tested using "Marconi M100" HPC machine (pwr9 + V100), 
ranked #9 among the most powerful HPC machines in the world, in May 2020. 
For any questions, please contact g.amati@cineca.it

To generate the executable for full sponge simulation, the correct commamand line is 
make mpi+openacc "PREPROC= -DORIGINAL -DOBSTACLES -DSPONGES -DOPENACC -DPGI" 

(please, refer to the README.txt file in the RUN folder to choose the right input; contact giacomo.falcucci@uniroma2.it for the full sponge geometry file [approximately, 1.5 Gb dat file]).

For the fluid dynamic "evolutionary" geometries, generate the executable running the following, correct commamand line:
make mpi+openacc "PREPROC= -DOPENACC -DPGI -DCYLINDER_INF -DOBSTACLES -DORIGINAL"

(please, refer to the README.txt file in the RUN folder to choose the right input; contact giacomo.falcucci@uniroma2.it for the various geometry files [approximately, 1 Gb dat file, each]).

