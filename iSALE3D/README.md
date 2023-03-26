For the development of iSALE3D code, it is carried out on Ububtu 20.04 using Fortran 90. The designed regridding method can be found in /iSALE3D/F90/write_dump.F90 from IRP repository. Relevant files that need to be edited in iSALE3D are also in the directory /iSALE3D/F90/. The preliminary test script and fields plotting scripts used for visualization can be found in the directory /iSALE3D/Python/. The input files used for simulations are in the directory /iSALE3D/inp/. 
For the whole iSALE3D code, please visit the branch ms_regird_3d of Github repository /isale-dev  
[https://github.com/isale-code/isale-dev/tree/ms_regrid_3d](url).  
The README.md file provides the instructions for pre-requisites. After the configuration completed, the iSALE3D code can be run under the dictionary /iSALE/example/demo3D/ with command similar to the following command:  
mpirun -n 4 ./iSALE3D -i asteriod.inp &  
Then to regrid the model, the command should be run under the same dictionary:  
mpirun -n 4 ./iSALE3D --ignore -i regrid.inp &                    
The visualisation scripts are written with python 3.8.
