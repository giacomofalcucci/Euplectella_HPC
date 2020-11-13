# Euplectella_HPC

This is the repository for Nature Publishing Policy.
The repository includes:
- the zipped .dat files for the considered geometeries (the Cylinder geometry is directly created within the code, no stl and dat files are prvided);
- the raw data and files to re-generate the figures;
- the CPU version of the LBM solver for the extreme simulations on the fluid dynamic evolutionary geometries;
- the GPU version of the LBM code for the extreme simulations on the entire Euplectella Aspergillum geometry.
For too large datasets, the complete code with the initial and boundary conditions is provided, with the directions and Paraview Files to re-generate the Figures.

The CPU version of the code requires gnu Fortran compiler with openMPI flag; the GPU version runs on nVidia accelerators and requires CUDA compilers. 
