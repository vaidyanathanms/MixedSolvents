# All directory data

hpcname   = 'pf' # pf or cades or kestrel or ANL hpcname
home_dir  = os.environ["HOME"]
head_dir  = home_dir + '/all_codes/files_mixsolv' # Super directory

# Gromacs directories - No need to change this
gmx_dir   = head_dir + '/src_gmx' # gmx file super directory
top_dir   = gmx_dir  + '/top_files' # topology dir
cfg_dir   = gmx_dir  + '/gropdb_files' # configuration dir
itp_dir   = gmx_dir  + '/itp_files' # prm/itp file dir
mdp_dir   = gmx_dir  + '/mdp_files' # mdp file dir

# PACKMOL directory
pack_exec = home_dir + '/tools/packmol/packmol' #packmol executable

# Run-directories
sh_dir    = head_dir + '/src_sh'  # sh file dir
scr_dir   = '/scratch/hpcl-phy191/vm5' # scratch (run) directory


