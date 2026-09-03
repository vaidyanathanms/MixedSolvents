#!/bin/bash
# Optimized code from Subil/Omar -  OLCFHELP-22622 [cades-help] Unable to run Gromacs utilities (Check email: June-03-2025)
# Instructions to build GROMACS are in this email
# Date: July-02-2025

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 06:30:00
#SBATCH -N 1                    # 1 node
#SBATCH --ntasks-per-node=1     # 1 task/node
#SBATCH -c 8                    # 8 cores/task
#SBATCH --mem=0
#SBATCH -J py_jobname
#SBATCH -o rdfout.%J
#SBATCH -e rdferr.%J

# Load modules
module reset
module load intel mkl fftw hwloc cmake

# Export gcc path and num_threads
export LD_LIBRARY_PATH=/sw/cades-open/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Initializing jobs
cd $SLURM_SUBMIT_DIR
echo "begin job @start time: ${date}"
echo $PWD

# Run command
gmx="${HOME}/gromacs-2024.5/install/bin/gmx_mpi "

# Make output directory
inpdir='anainps'
outdir='py_jobname'
mkdir -p ${outdir}

# Set python inputs
salt_cation=py_saltcation
salt_anion=py_saltanion
org_cation=py_orgcation
org_anion=py_organion
diluent=py_diluent
begtime=py_begtime
endtime=py_endtime

# Set flags
rdfflag=py_rdfcalc
difflag=py_difcalc

#----------------------------------------RDF calculations------------------------------------------------------------------------------

if (( rdfflag == 1 )); then
	# Compute RDF
	echo "Run rdf-${salt_cation}-all"
	srun ${gmx} rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_${salt_cation}_all.xvg -cn nrdf_${salt_cation}_all.xvg -b ${begtime} -e ${endtime} -rmpbc yes -ref -sf ${inpdir}/ref${salt_cation}.txt -sel -sf ${inpdir}/ref_all.txt
	wait

	echo "Run rdf-${org_cation}-all"
	srun ${gmx} rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_${org_cation}_all.xvg -cn nrdf_${org_cation}_all.xvg -b ${begtime} -e ${endtime} -rmpbc yes -ref -sf ${inpdir}/ref${org_cation}.txt -sel -sf ${inpdir}/ref_all.txt
	wait

	echo "All RDF calculations completed"
	echo "move files to ${outdir}"

	mv rdf_*xvg ${outdir} 
	mv nrdf_*xvg ${outdir}
	mv rdferr* ${outdir}
	mv rdfout* ${outdir}

fi;
#----------------------------------------Diffusivity calculations----------------------------------------------------------------------
if (( difflag == 1 )); then
	# Compute diffusivity
	echo "Run diffusivity-all"
	srun ${gmx} msd -f traj_npt_main.trr -s npt_main.tpr -o diff_all.xvg -b ${begtime} -e ${endtime} -rmpbc yes -sf ${inpdir}/ref_all.txt
	wait

	echo "All Diffusivity calculations completed"
	echo "move files to ${outdir}"

	mv diff_*xvg ${outdir}
fi; 
