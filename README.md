MemoScip is a MILP solver including a conflict free learning mechanism. 
This project is still in development and currently can only read wcsp files or very specific .lp and .opb files. 

The solver use a modified version of SCIP, you can find in SCIP_Patches the modified sources. I disabled the process tightening the bounds of slack variables and removed a mechanisme relocalising the root node.

To run the project do :

mkdir build
cd build
cmake .. -DSCIP_DIR=/path/to/SCIP"
make

Then to run the learning version :
./memoscipl myinstance.wcsp -l=1
or without learning:
./memoscipl myinstance.wcsp -l=0 
