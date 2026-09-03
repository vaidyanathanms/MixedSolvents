# Supporting files for genconf.py
# Version: Feb-17-2026
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import fileinput
import subprocess
from datetime import datetime 
from pathlib import Path
import fnmatch
import pandas as pd
from types import SimpleNamespace

#------------------------------------------------------------------
# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
    
#------------------------------------------------------------------

# Check for mdp files and copy/edit if not present
def check_cpy_mdp_files(srcdir,destdir,mdp_fyles,resname_arr,Tetau_vr=0.1\
                        ,Tetau_high_vr=0.1,Tetau_berend=0.1,Tetau_parrah = 0.2,\
                        Prtau_berend=4.0,Prtau_crescale=4.0,Prtau_parrah=8.0,ref_temp=300,\
                        hi_ref_temp=600,ref_pres=1.0):

    # Only temperatures need to have multiple coupling groups
    tc_grpdata      = " ".join(map(str, resname_arr))
    lenresname      = len(resname_arr)
    str_Tetau_vr    = " ".join([str(Tetau_vr)]*lenresname)
    str_Tetau_hi_vr = " ".join([str(Tetau_high_vr)]*lenresname)
    str_Tetau_ber   = " ".join([str(Tetau_berend)]*lenresname)
    str_Tetau_par   = " ".join([str(Tetau_parrah)]*lenresname)
    str_temp        = " ".join([str(ref_temp)]*lenresname)
    str_hi_temp     = " ".join([str(hi_ref_temp)]*lenresname)


    for mdp_fname in mdp_fyles:
        if not os.path.exists(srcdir + '/' + mdp_fname):
            raise RuntimeError(mdp_fname," not found in ",srcdir)
        gencpy(srcdir,destdir,mdp_fname)
        py_fname = mdp_fname
        rev_fname = mdp_fname.replace('_pyinp','')
        fr  = open(destdir + '/' + py_fname,'r')
        fw  = open(destdir + '/' + rev_fname,'w')
        fid = fr.read().replace("py_tcgrps",tc_grpdata).\
            replace("py_Temptau_vr",str_Tetau_vr).\
            replace("py_HighTemptau_vr",str_Tetau_hi_vr).\
            replace("py_Temptau_Berend", str_Tetau_ber).\
            replace("py_Temptau_ParRah",str_Tetau_par).\
            replace("py_Prestau_Berend",str(Prtau_berend)).\
            replace("py_Prestau_Crescale",str(Prtau_crescale)).\
            replace("py_Prestau_ParRah",str(Prtau_parrah)).\
            replace("py_ref_t",str_temp).\
            replace("py_Highref_t",str_hi_temp).\
            replace("py_ref_p",str(ref_pres))
        fw.write(fid)
        fw.close()
        fr.close()
        
#------------------------------------------------------------------

# Convert number from d.f to dpf
def convert_number(x):
    x = float(x)
    s = format(x,'g')

    if '.' not in s:
        s += '.0'

    return s.replace('.', 'p')

#------------------------------------------------------------------

# Match files without case sensitivity
def glob_ci(root, wildcard_pattern):
    wp = wildcard_pattern.casefold()
    root = Path(root)
    return [p for p in root.glob("*")
            if p.is_file() and fnmatch.fnmatch(p.name.casefold(), wp)]

#------------------------------------------------------------------

# Compute total number of electrolyte atoms
def compute_solv_molecules(s1name,s2name,v1frac,\
                           excel_fname='all_props.xlsx',\
                           sheet_name='PropertyData'):
    chemicals = load_components(excel_fname,sheet_name)
    s1nmol = 100*v1frac
    s1var  = chemicals[s1name]; s2var = chemicals[s2name]
    s1dens = s1var.density; s2dens = s2var.density
    s1MW   = s1var.MW; s2MW = s2var.MW
    s2nmol = s2var.PureSolvNMols if v1frac == 0 \
        else  round((s1nmol)*((100-v1frac)/(v1frac))*(s2dens/s1dens)*(s1MW/s2MW))
    volsol = s1nmol*(s1MW/s1dens) + s2nmol*(s2MW/s2dens)
    return s1nmol, s2nmol, volsol

#------------------------------------------------------------------

# Compute total number of electrolyte atoms
def compute_elec_molecules(solvol,conc):
    return round(solvol*conc*0.001) #0.001 for solvol in cc to l

#------------------------------------------------------------------

# Split cation and anion names and their valencies
def split_ename(ename,enmol,cation_types,anion_types,valence_tuple):
    
    catname = next((word for word in cation_types if word in ename), None)
    if catname is None:
        raise ValueError(f"No matching cation type found in {ename}.")
    anname  = next((word for word in anion_types if word in ename), None)
    if catname is None:
        raise ValueError(f"No matching cation type found in {ename}.")
    catmol = enmol*int(valence_tuple[0])
    anmol  = enmol*int(valence_tuple[1])
    
    return catmol,catname,anmol,anname

#------------------------------------------------------------------

# Compute total number of atoms and simulation box dimensions
def compute_sim_dims(s1name,s1nmol,s2name,s2nmol,cat1name,cat1mol,\
                     an1name,an1mol,v1tup,cat2name,cat2mol,an2name,\
                     an2mol,v2tup,e1nmol,e2nmol,init_dens = 0.8,\
                     excel_fname='all_props.xlsx',\
                     sheet_name='PropertyData'):

    chemicals = load_components(excel_fname,sheet_name)

    e1cat = cat1name if v1tup[0] == 1 else cat1name+str(v1tup[0])
    e1an  = an1name  if v1tup[1] == 1 else an1name+str(v1tup[1])
    e2cat = cat2name if v2tup[0] == 1 else cat2name+str(v2tup[0])
    e2an  = an2name  if v2tup[1] == 1 else an2name+str(v2tup[1])
    e1name = e1cat+e1an; e2name = e2cat+e2an

    s1var = chemicals[s1name]; s2var = chemicals[s2name]
    e1var = chemicals[e1name]; e2var = chemicals[e2name]
    
    totatoms = s1nmol*s1var.Natoms + s2nmol*s2var.Natoms + \
        e1nmol*e1var.Natoms + e2nmol*e2var.Natoms
    print(s1nmol,s2nmol,e1nmol,e2nmol)
    
    totmass = s1nmol*s1var.MW + s2nmol*s2var.MW + \
        e1nmol*e1var.MW + e2nmol*e2var.MW
    boxvol = totmass/(init_dens*0.6023) # conversion to Ang^3
    return totatoms,totmass,boxvol,boxvol**(1.0/3.0)

#------------------------------------------------------------------

# Generate arrays for writing GMX files
def generate_gmx_arrs(s1name,s1nmol,s2name,s2nmol,e1nmol,cat1name,cat1mol,\
                      an1name,an1mol,e2nmol,cat2name,cat2mol,an2name,an2mol):


    itp_arr = []; cfg_arr = []; resname_arr = []
    molname_arr = []; molval_arr = []
    
    if s1nmol != 0:
        itp_arr.extend([s1name])
        cfg_arr.extend([s1name])
        resname_arr.extend([s1name])
        molname_arr.extend([s1name])
        molval_arr.extend([s1nmol])
    if s2nmol != 0:
        itp_arr.extend([s2name])
        cfg_arr.extend([s2name])
        resname_arr.extend([s2name])
        molname_arr.extend([s2name])
        molval_arr.extend([s2nmol])
    if e1nmol != 0:
        itp_arr.extend([cat1name,an1name])
        cfg_arr.extend([cat1name,an1name])
        resname_arr.extend([cat1name,an1name])
        molname_arr.extend([cat1name,an1name])
        molval_arr.extend([cat1mol,an1mol])      
    if e2nmol != 0:
        itp_arr.extend([cat2name,an2name])
        cfg_arr.extend([cat2name,an2name])
        resname_arr.extend([cat2name,an2name])
        molname_arr.extend([cat2name,an2name])
        molval_arr.extend([cat2mol,an2mol])      

    return itp_arr,cfg_arr,resname_arr,molname_arr,molval_arr
    
#------------------------------------------------------------------

# Copy coordinate files
def cpy_coord_files(cfg_dir,inp_arr,workdir,destdirname="all_coords"):
    destdir = workdir + '/' + destdirname
    if not os.path.isdir(destdir):
        os.mkdir(destdir)


    coordfyle_list = [] 
    for fname in inp_arr:
        allfyles     = glob_ci(cfg_dir,fname+'.pdb')
        if allfyles == []:
            raise RuntimeError(f"No coordinate file for {fname} in {cfg_dir}")
        else:
            cname = max(allfyles,key=os.path.getmtime)
            shutil.copy2(cname,destdir)
            coordfyle_list.append(Path(cname).name)
    return coordfyle_list

#------------------------------------------------------------------

# Copy itp files
def cpy_itp_files(itp_dir,inp_arr,workdir,destdirname="itp_files"):
    destdir = workdir + "/" + destdirname
    extensions = [".itp",".btp",".prm"]
    itp_filelist = []
    if not os.path.isdir(destdir):
        os.mkdir(destdir)
    for fname in inp_arr:
        allfyles = []
        for ext in extensions:
            allfyles.extend(glob_ci(itp_dir,fname+ext))
        if allfyles == []:
            raise RuntimeError(f"No itp/prm for {fname}")
        for fyles in allfyles:
            shutil.copy2(Path(fyles),destdir)
            itp_filelist.append(Path(fyles).name)
    return itp_filelist

#------------------------------------------------------------------

# reorder files with atom types followed by btp and itp files
def reorder_itp_fnames(files):
    # define desired order
    order = {'.atp': 0, '.btp': 1, '.itp': 2, '.prm':3}
    sorted_files = sorted(files, key=lambda x:\
                          order.get(os.path.splitext(x)[1], 999))
    return sorted_files

#------------------------------------------------------------------   

# Generate topology files            
def gen_top_file(ff_indir,itp_fnames,molname_arr,molval_arr,workdir,\
                 destdirname = "itp_files",top_name="topol.top",\
                 nbfunc=1,comb_rule=2,gen_pairs='yes',\
                 fLJ=0.5,fQQ=0.5,sysname="system"):
    
    with open(workdir + "/" + top_name,'w') as ftop:
        now = datetime.now()
        ftop.write('; topology file generated by VMS\n')
        ftop.write('; ' + now.strftime("%B-%d-%Y, %I:%M:%S %p\n"))
        ftop.write('\n')

        ftop.write('[ defaults ]\n')
        ftop.write('; nbfunc\tcomb-rule\tgen-pairs\tfudgeLJ\tfudgeQQ\n')
        ftop.write(f'{nbfunc}\t{comb_rule}\t{gen_pairs}\t{fLJ}\t{fQQ}\n')

        ftop.write('\n')
        # add atomtypes and parameters before [ moleculetype ]
        inc_top = ''
        inc_pre = '#include '
        for i in range(len(itp_fnames)):
            ftop.write(inc_top + inc_pre +  '"./' + destdirname +\
                       "/" + itp_fnames[i] + '"\n')

        # Add system name
        ftop.write('\n')
        ftop.write('[ system ]\n')
        ftop.write(f'{sysname}\n\n')

        # Add molecule names and values
        ftop.write('[ molecules ]\n')
        for i in range(len(molname_arr)):
            ftop.write(f'{molname_arr[i]}\t{int(molval_arr[i])}\n')

#------------------------------------------------------------------

# Set up packmol files
def setup_packmol(inpfnames,molval_arr,workdir,boxl,\
                  destdirname="all_coords",outname='mixture.pdb',\
                  packname='make_mixture.inp'):

    with open(workdir + '/' + packname,'w') as fpack:
        fpack.write('# \n # A mixture of ionic liquids & diluents or pure ionic liquids\n')
        fpack.write('# Min separation: 2 Ang\n')
        fpack.write('# input filetype: pdb\n')

        fpack.write(f'tolerance\t2.0\n')
        fpack.write(f'filetype\tpdb\n')
        fpack.write(f'output\t{outname}\n\n')

        # Add each molecule with their type
        for i,fname in enumerate(inpfnames):
            molname = fname.split(".")[0]
            fpack.write(f'# Add {molname}\n')
            fpack.write(f'structure\t ./{destdirname}/{fname}\n')
            fpack.write(f'number\t{str(int(molval_arr[i]))}\n')
            fpack.write(f'inside box 0. 0. 0. {boxl} {boxl} {boxl}\n')
            fpack.write(f'end structure\n\n')

#------------------------------------------------------------------

def make_anainps(workdir,saltcatname,orgcatname,resnames,\
                 rdftuple=list(),rdfcalc=1,diffcalc=1,\
                 autocorrcalc=1):
    if rdfcalc:
        if not os.path.isdir(workdir + '/anainps'):
            os.mkdir(workdir + '/anainps')


        with open(workdir + '/anainps/ref_all.txt','w') as frdf:
            for res in range(len(resnames)-1):
                frdf.write(f'resname {resnames[res]};\n')
            frdf.write(f'resname {resnames[-1]}')
            
            if rdftuple:
                for res,typ in rdftuple:
                    frdf.write(f';\nresname {res} and type {typ}')

            rdfsaltcat = '/anainps/ref'+saltcatname+'.txt'
            with open(workdir + rdfsaltcat,'w') as fs:
                fs.write(f'resname {saltcatname}')

            rdforgcat = '/anainps/ref'+orgcatname+'.txt'
            with open(workdir + rdforgcat,'w') as fs:
                fs.write(f'resname {orgcatname}')

#------------------------------------------------------------------

# Copy and run shell script files
def cpy_anash_files(srcdir,destdir,sh_ana_fyle,saltcat,\
                    saltan,orgcat,organ,diluent,begtime,endtime,\
                    rdfcalc=1,diffcalc=1,runall=0,jname="default"):

    gencpy(srcdir,destdir,sh_ana_fyle)
    py_fname = destdir + '/' + sh_ana_fyle
    rev_fname = py_fname.replace('_pyinp','')
    fr  = open(py_fname,'r')
    fw  = open(rev_fname,'w')
    fid = fr.read().replace("py_jobname","ana_"+jname).\
        replace("py_saltcation",saltcat).\
        replace("py_saltanion",saltan).\
        replace("py_orgcation",orgcat).\
        replace("py_organion",organ).\
        replace("py_diluent",diluent).\
        replace("py_begtime",str(begtime)).\
        replace("py_endtime",str(endtime)).\
        replace("py_rdfcalc",str(rdfcalc)).\
        replace("py_difcalc",str(diffcalc))
    fw.write(fid)
    fr.close(); fw.close()

    if runall:
        os.chdir(destdir)
        print(f"Submitting {rev_fname} ..")
        subprocess.call(["chmod", "777", rev_fname])
        subprocess.call(["sbatch", rev_fname])
#------------------------------------------------------------------

# Copy and run shell script files
def cpy_sh_files(srcdir,destdir,sh_pp_fyle,sh_md_fyle,boxl,\
                 runall=0,outname='mixture.pdb',\
                 packname='make_mixture.inp',top_fyle="topol.top",\
                 jname="default",qname='-A chem',pname='burst',\
                 hpcname='cades'):

    # Check for tpr files to provide run conditions
    if glob.glob(destdir+'/*.tpr') == []:
        print('No tpr files found. Beginning new runs..')
        if not os.path.exists(srcdir+'/' + sh_pp_fyle):
            raise RuntimeError(sh_pp_fyle+" not found in "+\
                               srcdir)
        gencpy(srcdir,destdir,sh_pp_fyle)

        py_fname = destdir + '/' + sh_pp_fyle
        rev_fname = py_fname.replace('_pyinp','')
        fr  = open(py_fname,'r')
        fw  = open(rev_fname,'w')
        fid = fr.read().replace("py_qname",qname).\
            replace("py_hpcname",hpcname).\
            replace("py_pname",pname).\
            replace("py_jobname","pp_"+jname).\
            replace("py_packname",packname).\
            replace("py_confname",outname).\
            replace("py_boxl",str(boxl+0.2)).\
            replace("py_topname",top_fyle)
        fw.write(fid)
        fr.close(); fw.close()
        
    else:
        print('tpr files found. Continuing runs..')
        if not os.path.exists(srcdir+'/' + sh_md_fyle):
            raise RuntimeError(sh_md_fyle+" not found in "+\
                               srcdir)
        gencpy(srcdir,destdir,sh_md_fyle)

        py_fname = destdir + '/' + sh_md_fyle
        rev_fname = py_fname.replace('_pyinp','')
        fr  = open(py_fname,'r')
        fw  = open(rev_fname,'w')
        fid = fr.read().replace("py_qname",qname).\
            replace("py_hpcname",hpcname).\
            replace("py_pname",pname).\
            replace("py_jobname",'md_'+jname).\
            replace("py_topol",top_fyle)
        fw.write(fid)
        fr.close(); fw.close()

        
    if runall:
        os.chdir(destdir)
        print(f"Submitting {rev_fname} ..")
        subprocess.call(["chmod", "777", rev_fname])
        subprocess.call(["sbatch", rev_fname])
#------------------------------------------------------------------

def clean_up(workdir):
    all_files = glob.glob(workdir + '/*_pyinp*')
    if not all_files == []:
        for fyle in all_files:
            os.remove(fyle)
#------------------------------------------------------------------

def load_components(excel_file, sheet_name='PropertyData'):
    df = pd.read_excel(excel_file, sheet_name=sheet_name)

    df = df.rename(columns={
        "Molecule name": "molname",
        "Abbreviation": "abbreviation",
        "# of atoms": "Natoms",
        "MW (g/mol)": "MW",
        "density (g/cc)": "density",
        "Pure_Solv_#_molecules": "PureSolvNMols",
    })

    df = df.dropna(subset=["abbreviation"])

    components = {}

    for _, row in df.iterrows():
        abbreviation = str(row["abbreviation"]).strip()

        components[abbreviation] = SimpleNamespace(
            Natoms=int(row["Natoms"]),
            MW=float(row["MW"]),
            density=float(row["density"]),
            PureNMols=float(row["PureSolvNMols"]),
        )

    return components
#------------------------------------------------------------------
