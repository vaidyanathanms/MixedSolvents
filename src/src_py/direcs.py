# All directory data

import os

#HPC details
hpcname   = 'pf' # pf or cades or kestrel or ANL hpcname
qname     = '-q normal' # name of the queue with prefix
pname     = 'serial' # sub-queue name

# Main paths
home_dir  = os.environ["HOME"]
scr_dir   = os.environ["PFSCRATCH"]
head_dir  = home_dir + '/all_codes/files_mixsolv' # Super directory

# Gromacs directories - No need to change this
gmx_dir   = head_dir + '/src/src_gmx' # gmx file super directory
top_dir   = gmx_dir  + '/top_files' # topology dir
cfg_dir   = gmx_dir  + '/gropdb_files' # configuration dir
itp_dir   = gmx_dir  + '/itp_files' # prm/itp file dir
mdp_dir   = gmx_dir  + '/mdp_files' # mdp file dir

excel_file = head_dir + '/all_props.xlsx'

# PACKMOL directory
pack_exec = home_dir + '/tools/packmol/packmol' #packmol executable

# Run-directories
sh_dir    = head_dir + '/src/src_sh'  # sh file dir
scr_dir   = '/scratch/hpcl-phy191/vm5' # scratch (run) directory


