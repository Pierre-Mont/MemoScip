MemoScip is a MILP solver including a conflict free learning mechanism. 
This project is still in development and currently can only read wcsp files or very specific .lp and .opb files. 

The solver use a modified version of SCIP, you can find in SCIP_Patches the modified sources. 

pricestore.c -> hack to always add slack as priced variables. 
lp.c -> force solving even when farkas proof is incorrect
prop_redcost.c -> correctly detect reduced cost filtering 
prop_rootredcost.c -> correctly detect reduced cost filtering
cons_linear.c -> Do not tighten slack variables
tree.c -> disable root node relocalisation



To run the project do :

mkdir build
cd build
cmake .. -DSCIP_DIR=/path/to/SCIP"
make

Then to run the learning version :
./memoscipl myinstance.wcsp -l=1
or without learning:
./memoscipl myinstance.wcsp -l=0 
