#!/bin/bash
# Optimized code from Subil/Omar -  OLCFHELP-22622 [cades-help] Unable to run Gromacs utilities (Check email: June-03-2025)
# Instructions to build GROMACS are in this email
# Date: July-02-2025

#SBATCH py_qname
#SBATCH -p py_pname
#SBATCH -t 23:00:00
#SBATCH -N 1                    # 1 nodes
#SBATCH --ntasks-per-node=1     # 1 tasks/node
#SBATCH -c 48                   # 48 cores/task
#SBATCH --mem=0
#SBATCH -J py_jobname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

hpcname=py_hpcname
if [[ ${hpcname} == "pf" ]]; then
    # Load modules
    module reset
    module load mkl fftw hwloc cmake gromacs
    gmx="gmx_mpi "
else
    # Load modules
    module reset
    module load intel mkl fftw hwloc cmake
    # Export gcc path and num_threads
    export LD_LIBRARY_PATH=/sw/cades-open/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
    # Run command
    gmx="${HOME}/gromacs-2024.5/install/bin/gmx_mpi "
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Initializing jobs
cd $SLURM_SUBMIT_DIR
echo "begin job @start time: ${date}"
echo $PWD

mkdir -p initdir
mkdir -p outdir
mkdir -p trajfiles

#---------------------------------------------------------Generate initial run----------------------------------
finit_inp=./enermin.tpr
if ! test -f "$finit_inp"; then
	echo "begin generating enermin.tpr.."
	# generate enermin files
	srun ${gmx} grompp -f minim.mdp -c main.pdb -p topol.top -o enermin.tpr
	
fi
wait

#---------------------------------------------------------Minimization -----------------------------------------
fmin_inp=./nvt.tpr 
if ! test -f "$fmin_inp"; then
	
	echo "begin running enermin.tpr.."
	# run enermin.tpr
	srun ${gmx} mdrun -s enermin.tpr -cpo state_min.cpt -cpi state_min.cpt -cpt 2 -g md_min.log -o traj_min.trr -e ener_min.edr -c confout_min.gro  
	wait

	echo "begin generating nvt.tpr.."
	# generate nvt files
	srun ${gmx} grompp -f nvt.mdp -c confout_min.gro -p topol.top -o nvt.tpr 
	wait

        cp md_min.log trajfiles/md_min.log
	cp ener_min.edr trajfiles/ener_min.edr
	cp confout_min.gro trajfiles/confout_min.gro
	
fi
wait

#----------------------------------------------------------NVT Equilibration--------------------------------------------
fnvt_inp=./npt_c_rescale.tpr
if ! test -f "$fnvt_inp"; then

	echo "begin running nvt.tpr.."
	# run nvt.tpr
	srun ${gmx} mdrun -s nvt.tpr -cpo state_nvt.cpt -cpi state_nvt.cpt -cpt 5 -g md_nvt.log -o traj_nvt.trr -e ener_nvt.edr -c confout_nvt.gro  -notunepme
	wait

	echo "begin generating npt_c_rescale.tpr.."
	# generate npt_crescale files
	srun ${gmx} grompp -f npt_crescale.mdp -c confout_nvt.gro -p topol.top -o npt_c_rescale.tpr
	wait

	cp md_nvt.log trajfiles/md_nvt.log
	cp ener_nvt.edr trajfiles/ener_nvt.edr
	cp confout_nvt.gro trajfiles/confout_nvt.gro
fi
wait

#----------------------------------------------------------NPT Equilibration/Production--------------------------------------------
fnpt_inp=./npt_main.tpr
if ! test -f "$fnpt_inp"; then

	echo "begin running npt_c_rescale.tpr.."
	# run npt_c_rescale.tpr
	srun ${gmx} mdrun -s npt_c_rescale.tpr -cpo state_npt_c_rescale.cpt -cpi state_npt_c_rescale.cpt -cpt 5 -g md_npt_c_rescale.log -o traj_npt_c_rescale.trr -e ener_npt_c_rescale.edr -c confout_npt_c_rescale.gro -notunepme 
	wait

	echo "begin generating npt_main.tpr.."
	# generate npt_main files
	srun ${gmx} grompp -f npt_main.mdp -c confout_npt_c_rescale.gro -p topol.top -o npt_main.tpr 
	wait

	cp md_npt_c_rescale.log trajfiles/md_npt_c_rescale.log
	cp ener_npt_c_rescale.edr trajfiles/ener_npt_c_rescale.edr
	cp confout_npt_c_rescale.gro trajfiles/confout_npt_c_rescale.gro
	wait

        echo "begin running npt_main.tpr.."
        # run npt_main.tpr
        srun ${gmx} mdrun -s npt_main.tpr -cpo state_npt_main.cpt -cpi state_npt_main.cpt -cpt 5 -g md_npt_main.log -o traj_npt_main.trr -e ener_npt_main.edr -c confout_npt_main.gro -notunepme       
	wait

        cp md_npt_main.log trajfiles/md_npt_main.log
        cp ener_npt_main.edr trajfiles/ener_npt_main.edr
        cp confout_npt_main.gro trajfiles/confout_npt_main.gro

else

        echo "begin running npt_main.tpr.."
        # run npt_main.tpr
        srun ${gmx} mdrun -s npt_main.tpr -cpo state_npt_main.cpt -cpi state_npt_main.cpt -cpt 5 -g md_npt_main.log -o traj_npt_main.trr -e ener_npt_main.edr -c confout_npt_main.gro -notunepme 
	wait


        cp md_npt_main.log trajfiles/md_npt_main.log
        cp ener_npt_main.edr trajfiles/ener_npt_main.edr
        cp confout_npt_main.gro trajfiles/confout_npt_main.gro
fi
wait
#------------------------------------------------------------------------------------------------------------------------------------

echo "End of run.."
