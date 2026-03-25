#!/bin/bash
#neus=(100)
for ((j=1; j<=100; j+=1)); do
echo '#!/bin/bash' >> grid_750_$j.sh
echo '#SBATCH -c 1' >> grid_750_$j.sh
echo '#SBATCH --time=500:00:00' >> grid_750_$j.sh
echo '#SBATCH --mail-user=fabian.hernandez.13@sansano.usm.cl' >> grid_750_$j.sh
echo "#SBATCH --mail-type=BEGIN,END,FAIL" >> grid_750_$j.sh
echo 'source /user/f/fhpinto/projects/bin/activate' >> grid_750_$j.sh
echo 'deltams=(0.142 0.15 0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4)' >> grid_750_$j.sh
echo "xis=(5)" >> grid_750_$j.sh
echo "lams=(0.1)" >> grid_750_$j.sh
echo 'export PYTHIA8DATA=/user/f/fhpinto/RPVLLP/HEPsoftware/pythia8310/share/Pythia8/xmldoc/' >> grid_750_$j.sh
echo 'export LD_LIBRARY_PATH=/user/f/fhpinto/RPVLLP/HEPsoftware/fastjet-install/lib/:$LD_LIBRARY_PATH' >> grid_750_$j.sh
echo 'export LD_LIBRARY_PATH=/user/f/fhpinto/RPVLLP/HEPsoftware/pythia8310/lib/:$LD_LIBRARY_PATH' >> grid_750_$j.sh
echo 'for deltam in "${deltams[@]}"; do' >> grid_750_$j.sh
echo 'for lam in "${lams[@]}"; do' >> grid_750_$j.sh
echo "v0=$(echo "10 * $j" | bc)" >> grid_750_$j.sh
echo "folder=$(echo "$j" | bc)" >> grid_750_$j.sh
echo "cd /user/f/fhpinto/projects/MG5_aMC_v2_9_21" >> grid_750_$j.sh
echo './bin/mg5_aMC scripts/vquint/vquint-$v0-$deltam-$lam' >> grid_750_$j.sh
#echo './disappearingT-ATLAS $v0 $deltam $lam $folder' >> grid_750_$j.sh
echo 'done' >> grid_750_$j.sh
echo 'done' >> grid_750_$j.sh
done
