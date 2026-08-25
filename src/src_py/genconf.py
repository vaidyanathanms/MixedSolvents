# Generate configurations for mixed-solvent simulations
# Author: Vaidyanathan M. S

# Version: Aug-24-2026
#------------------------------------------------------------------
#Import modules 
import os
import sys
import numpy as np
import re
import shutil
import glob
import math
import subprocess
import extra_functions as ef
import directories
#------------------------------------------------------------------ 
# Version Info
print("Generating GROMACS run-time inputs")
print("Version: Aug-18-2026")

#------------------------------------------------------------------
# Input data
# Note: 100% water is kept at 99.99% to avoid NaN errors
run_all  = 0 # 1-copy files and run, 0-NO run (copies files)

cation_types      = ['Fe','Na']
anion_types       = ['SO4','Cl']
electrolyte_1_arr = ['FeSO4','FeCl']
valence_1_arr     = [(1,1),(1,2)]
electrolyte_2_arr = ['NaSO4','NaCl']
valence_2_arr     = [(2,1),(1,1)]

solvent_1_arr     = ['WAT']
solvent_2_arr     = ['EGL','MET','TGL','THF','ACN']
sol1_volfracs_arr = np.array([100,90,75,50,25,0]) # in %
sol1_conc_arr     = np.array([0.1, 1]) # in mol/l
sol2_conc_arr     = np.array([0, 0.05]) # in mol/l

#------------------------------------------------------------------
# Required GMX/sh and default gro/top/itp files
attye_fname = 'ffnonbonded.itp'
mdp_fyles   = ['minim_pyinp.mdp','nvt_pyinp.mdp',\
               'npt_crescale_pyinp.mdp','npt_main_pyinp.mdp']
sh_md_fyle  = 'run_md_pyinp.sh'
sh_pp_fyle  = 'run_preprocess_pyinp.sh'
sh_rep_fyl  = ['repeat_all.sh','repeat_md.sh']
def_inicon  = 'initconf.gro'

#------------------------------------------------------------------
# Simulation details
init_density      = 0.8 # in g/cc

#------------------------------------------------------------------
# Check directories
curr_dir = os.getcwd()
if not os.path.isdir(main_dir):
    raise RuntimeError('Check path to src files and update main_dir')
if not os.path.isdir(scr_dir):
    raise RuntimeError('Check path to working dir and update scr_dir')

#------------------------------------------------------------------
# Create head scratch directory
scr_maindir = scr_dir + '/mix_solv_systems'
if not os.path.isdir(scr_maindir):
    os.mkdir(scr_maindir)
    
#------------------------------------------------------------------
#------------------Main-loop starts here---------------------------
# Loop over solvent-1
for s1id,s1name in enumerate(solvent_1_arr):
                        
    # Loop over solvent-2
    for s2id,s2name in enumerate(solvent_2_arr):
        
        # Loop over volume fractions of solvent-1
        for v1id,v1frac in enumerate(sol1_volfracs_arr):
            # Compute required number of molecules of each solvent
            s1nmol,s2nmol,volsol = ef.compute_solv_molecules(s1name,\
                                                             s2name,\
                                                             v1frac)

            # Loop over concentrations of solvent-1
            for c1id,c1val in enumerate(sol1_conc_arr):

                # Loop over concentrations of solvent-2
                for c2id,c2val in enumerate(sol2_conc_arr):

                    # Loop over electrolyte-1
                    for e1id,e1name in enumerate(electrolyte_1_arr):
                        # Split cations and anions in Elec-1
                        catmol,catname,anmol,anname = \
                            ef.split_ename(e1name,e1nmol,\
                                           cation_types,\
                                           anion_types,\
                                           valence_1_arr[e1id])

                        # Compute number of molecules of electrolyte-1
                        e1nmol = ef.compute_elec_molecules(volsol,\
                                                           c1val,\
                                                           e1name)
                        # Loop over electrolyte-2
                        for e2id,e2name in enumerate(electrolyte_2_arr):
                            # Split cations and anions in Elec-1
                            catmol,catname,anmol,anname = \
                                ef.split_ename(e2name,e2nmol,\
                                               cation_types,\
                                               anion_types,\
                                               valence_2_arr[e2id])

                            # Compute number of molecules of electrolyte-2
                            e2nmol = ef.compute_elec_molecules(volsol,\
                                                               c2val,\
                                                               e2name)
                            elec1dir = scr_dir + '/' + e1name + \
                                '_conc_' + str(c1val)
                            if not os.path.isdir(elec1dir):
                                os.mkdir(elecdir)

                            elec2dir = elec1dir + '/' + e2name + \
                                '_conc_' + str(c2val)
                            if not os.path.isdir(elec2dir):
                                os.mkdir(elec22dir)

                            solvdir = elec2dir + '/' + s1name + '_' s2name
                            if not os.path.isdir(solvdir):
                                os.mkdir(solvdir)

                            workdir = solvdir + '/volfrac_' + str(v1frac)
                            if not os.path.isdir(workdir):
                                os.mkdir(solvdir)
                            
                            print(f'Setting up simulations for
                            Solvent-1: {s1name}; Solvent-2: {s2name};
                            Elec-1: {e1name}; Elec-2: {e2name}')  
                                
                            # Set-up box dimensions 
                            totnatoms,totmass,boxvol,boxlen = \
                                ef.compute_sim_dims(s1name,s1nmol,\
                                                    s2name,s2nmol,\
                                                    e1name,e1nmol,\
                                                    e2name,e2nmol,\
                                                    init_dens = init_density)
                            

                            # Set-up gmx inputs
                            itp_arr,cfg_arr,resname_arr,molname_arr,molval_arr = \
                                ef.generate_gmx_arrs(s1name,s1nmol,s2name,s2nmol,\
                                                     e1name,e1nmol,e2name,e2nmol)

                            sysname = f"{e1name}_{c1val:.1f}_{e2name}_{c2val:.1f}_{s1name}_{volfrac.2f}_{s2name}"

                            
                            # Copy and edit mdp files
                            print('Copying and editing mdp files ...')
                            ef.check_cpy_mdp_files(mdp_dir,workdir,mdp_fyles,\
                                                   resname_arr,Tetau_vr=0.1,\
                                                   Tetau_high_vr=0.1,\
                                                   Tetau_berend=0.1,\
                                                   Tetau_parrah = 0.2,\
                                                   Prtau_berend=4.0,\
                                                   Prtau_crescale=4.0,\
                                                   Prtau_parrah=8.0,\
                                                   ref_temp=300,hi_ref_temp=600,\
                                                   ref_pres=1.0)

                            # Copy coordinate files
                            print('Copying coordinate files ...')
                            coord_fnames = \
                                ef.cpy_coord_files(cfg_dir, cfg_arr, workdir,\
                                                   destdirname="all_coords")


                            # Copy and edit itp files
                            print('Copying itp or prm files ...')
                            itp_fnames = ef.cpy_itp_files(itp_dir,\
                                                          itp_arr,workdir,\
                                                          destdirname="itp_files")
                            itp_fnames = ef.reorder_itp_fnames(itp_fnames)
                            
                            # Copy atype files
                            ef.gencpy(itp_dir,workdir+'/itp_files',attype_fname) 
                            itp_fnames.insert(0,attype_fname)

                            # Generate top files
                            print('Generating top file ...')
                            ef.gen_top_file(itp_dir,itp_fnames,molname_arr,\
                                            molval_arr,workdir,\
                                            destdirname="itp_files",\
                                            top_name="topol.top",\
                                            nbfunc=1,comb_rule=3,gen_pairs='yes',\
                                            fLJ=0.5,fQQ=0.5,sysname=sysname)

                            # Generating packmol files
                            print('Generating packmol file ...')
                            ef.setup_packmol(coord_fnames,molval_arr,workdir,\
                                             box_arr_ang[iarr],\
                                             destdirname="all_coords",\
                                             outname='mixture.pdb',\
                                             packname = 'make_mixture.inp')
                            
                            # Run files
                            print('Generating shell script files ...')
                            ef.cpy_sh_files(sh_dir,workdir,sh_pp_fyle,\
                                            sh_md_fyle,box_arr_nm[iarr],\
                                            runall=run_all,\
                                            outname='mixture.pdb',\
                                            packname='make_mixture.inp',\
                                            top_fyle= "topol.top",\
                                            jname=sysname,qname=qname,\
                                            pname=pname,hpcname=hpcname)

                            
                            # Cleaning up
                            print('Cleaning up directory ..')
                            ef.clean_up(workdir)
