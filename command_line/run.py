# LIBTBX_SET_DISPATCHER_NAME phenix.cryo_fit
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT 

# Steps of cryo_fit:
# 1_Make_gro
# 2_Clean_gro
# 3_Prepare_to_Minimize
# 4_Minimize
# 5_Make_restraints
# 6_Make_0_charge
# 7_Make_tpr_for_EM_map_fitting
# 8_EM_map_fitting_itself
# 9_Arrange_outputs (including draw_a_figure_of_cc)

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from iotbx import file_reader
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out
import mmtbx.utils
import platform
import shutil # for rmdir

# this is needed to import all common functions
path = subprocess.check_output(["which", "phenix.cryo_fit"])
splited = path.split("/")
command_path = ''
for i in range(len(splited)-3):
  command_path = command_path + splited[i] + "/"
command_path = command_path + "modules/cryo_fit/"
common_functions_path = command_path + "common_functions/"
sys.path.insert(0, common_functions_path)
from common_functions import *

legend = """\
Goal
    - Changes an input biomolecule structure to fit into the cryo-EM map
    
How to use
    - phenix.cryo_fit <input cif/pdb file> <input map file>
    - Don't run at a phenix folder such as /Users/<user>/bin/phenix-dev-2906/modules/cryo_fit

Input:
    - A .cif or .pdb file
         A template/starting structure that is aligned to a target cryo EM 
         density map structurally (for example by USCF Chimera)
    - A .ccp4 (MRC) or .map (MRC) or .sit (Situs) file, a cryo EM density map 
    
Output folder:
    - cryo_fitted.x: Fitted biomolecule structure to a target cryo-EM map
    - Correlation coefficient record: Record of correlation coefficient between cryo-EM map
    and current biomolecule structure
      
Usage example with minimum input requirements (all other options will run with default values):
    - phenix.cryo_fit GTPase_activation_center.pdb GTPase_activation_center.map

Usage example with step 7~8 only
    - phenix.cryo_fit GTPase_activation_center.pdb GTPase_activation_center.map step_1=False step_2=False step_3=False step_4=False step_5=False step_6=False
    
Most useful options (GUI has more explanation about these):
    - number_of_cores_to_use
    - number_of_steps_for_minimization
    - emweight_multiply_by
    - number_of_steps_for_cryo_fit
"""

#master_params_str are used for default values of options in GUI
master_params_str = """
cryo_fit {
include scope libtbx.phil.interface.tracking_params
Input{
  model_file_name = None
    .type = path
    .short_caption = Starting model file 
    .multiple = False
    .help = Either a homology model or a model from different organism/experimental method (.cif/.pdb)
    .style = file_type:pdb bold input_file
  map_file_name = None
    .type = path
    .short_caption = Target map file 
    .help = Cryo-EM map file (.ccp4/.map/.sit)
    .style = bold input_file
  cryo_fit_path = None
    .type = path
    .short_caption = gromacs_cryo_fit executable path
    .help = Path that has gromacs_cryo_fit executables (such as mdrun). For example, /Users/doonam/bin/cryo_fit/bin
    .style = bold directory
}
Steps
{
  step_1 = True
    .type = bool
    .short_caption = 1. Make gro file
    .help = Make gro file from user input pdb file
  step_2 = True
    .type = bool
    .short_caption = 2. Clean gro file
    .help = Clean gro file to be compatible for Amber03 forcefield
  step_3 = True
    .type = bool
    .short_caption = 3. Prepare to minimize a starting structure
    .help = Make a tpr file for mimization
  step_4 = True
    .type = bool
    .short_caption = 3. Minimize a starting structure
    .help = Minimize a starting structure to avoid blow-up later
  step_5 = True
    .type = bool
    .short_caption = 4. Make contact potential
    .help = Make contact potential to better keep secondary structure
  step_6 = True
    .type = bool
    .short_caption = 5. Neutralize
    .help = Neutralize charges
  step_7 = True
    .type = bool
    .short_caption = 6. Make a tpr file
    .help = Make tpr file for cryo_fit
  step_8 = True
    .type = bool
    .short_caption = 7. Fit to a cryo-EM map
    .help = Fit a user given structure into a user given cryo-EM map
}
Options
{
  restraint_algorithm_minimization = *default none none_default
    .type = choice
    .help = default will be "default" \
            If a user sees this error during minimization Too many lincs warnings try either none or none_default. \
            none_default will minimize twice. \
            e.g. first with restraint_algorithm_minimization = none \
            then with restraint_algorithm_minimization = default
  time_step_for_minimization = 0.001
    .type = float
    .help = Default value is 0.001. Try 0.0005 if a user see this error during minimization
    "Fatal error: A charge group moved too far between two domain decomposition steps \
    This usually means that your system is not well equilibrated"
  number_of_steps_for_minimization = None
    .type = int
    .short_caption = Number of steps for minimization
    .help = Specify number of steps for minimization. \
           If this is left blank, cryo_fit will estimate it depending on molecule size. \
           Enough minimization will prevent "blow-up" during MD simulation later.
  emsteps          = None
    .type          = int
    .short_caption = EM steps
    .help = emsteps is the number of integration steps between re-evaluation of the simulated map and forces. \
            The longer the emsteps be, the faster overall cryo_fit running time. \
            If it is left blank, the cryo_fit will automatically determine the emsteps
  emweight_multiply_by = 8
    .type          = int
    .short_caption = Multiply EM weight by this number
    .help = Multiply by this number to the number of atoms for weight for cryo-EM map bias. \
            Default value is 8 (Paul/Serdal recommended 2, at 2~8, tRNA lost base-pairs according to phenix.secondary_structure_restraints. However, visual inspection shows all good base-pairs.). \
            For example, emweight = (number of atoms in gro file) x (emweight_multiply_by) \
            The higher the weight, the stronger bias toward EM map rather than MD force field and stereochemistry preserving restraints. \
            If user's map has a better resolution, higher value of emweight_multiply_by is recommended since map has much information. \
            If user's map has have a worse resolution, lower value of emweight_multiply_by is recommended for more likely geometry. \
            If CC (correlation coefficient) needs to be improved faster, higher number of emweight_multiply_by is recommended.
  emwritefrequency = None
    .type          = int
    .short_caption = EM write frequency
    .help          = Frequency with which the simulated maps are written to file. \
                     If this frequency is too small, it can cause extremely large amounts of data to be written.\
                     If it is left blank, the cryo_fit will use default value of 1,000,000
  no_rerun         = False
    .type          = bool
    .short_caption = No automatic re-run of cryo_fit
    .help          = If checked/true, cryo_fit does not re-run even its cc values kept increasing.
  nstxtcout        = 100
    .type          = int
    .short_caption = Frequency for trajectory
    .help          = A frequency to write coordinates to xtc trajectory. \
                     By default, this is 100 (both commandline and GUI).
  number_of_steps_for_cryo_fit = None
    .type                      = int
    .short_caption             = Number of steps for the 1st iteration of cryo_fit
    .help = This is the initial number of steps for the 1st iteration of cryo_fit. \
            Eventually, cryo_fit will increase it iteratively until it reaches cc plateau. \
            This value should be > 100. \
            (Just for tutorial files, this will be a fixed value, e.g. 70,000, unless a user specifies it)
  time_step_for_cryo_fit = 0.002
    .type                = float
    .short_caption       = Time step for MD simulation during cryo_fit
    .help = Default value is 0.002. Try 0.001 if a user see this error during cryo_fit \
      "Fatal error: A charge group moved too far between two domain decomposition steps \
      This usually means that your system is not well equilibrated"
}

many_step_____n__dot_pdb = False
  .type = bool
  .short_caption = If true, emweight_multiply_by=1, lincs_order=1, annealing_gen_temp=40
devel = False
  .type = bool
  .short_caption = If true, just quick check for sanity
force_field = *amber03 gromos96 
    .type = choice
    .short_caption = Force field
    .help = We plan to support amber03 more to include more ligands 
ignh = True
  .type = bool
  .short_caption = If true, ignore hydrogen atoms that are in the coordinate file
  .help = http://manual.gromacs.org/programs/gmx-pdb2gmx.html
initial_cc_wo_min = False
    .type = bool
    .short_caption = It is useful when cc_cryo_fit only is required.
initial_cc_w_min = False
    .type = bool
    .short_caption = It is useful when cc_cryo_fit only is required.
kill_mdrun_mpirun_in_linux = False
  .type = bool
  .short_caption = If true, kill any existing md run and mpirun.
  .help = This command works only for Linux
lincs_order = None
  .type = int
  .short_caption = LINear Constraint Solver
  .help = The accuracy in set with lincs-order, which sets the number of matrices in the expansion for the matrix inversion. \
          If it is not specified, the cryo_fit will use 4.
max_emweight_multiply_by = 7
                   .type = int
                   .short_caption = Up to this value, emweight_multiply_by keeps increasing by 2 times
missing = True
  .type = bool
  .short_caption = If true, Continue when atoms are missing, dangerous
  .help = http://manual.gromacs.org/programs/gmx-pdb2gmx.html
ns_type = *grid simple
  .type = choice
  .short_caption = Method to determine neighbor list (simple, grid) during minimization
  .help = "Grid" is needed for domain decomposition (dd) for faster execution and ran well with GTPase_activation_center, \
          beta-galactosidase, and nucleosome, but "simple" was needed for trouble-shooting of 80S ribosome.
number_of_cores_to_use = 2 4 8 12 16 24 32 *max
  .type = choice
  .short_caption = Number of cores to use for minimization and cryo_fit
  .help = Specify number of cores for minimization and cryo_fit. \
          If it is not specified, or max is chosen, the cryo_fit will try to use most cores automatically (up to 16)
perturb_xyz_by = 0.05
  .type = float
  .short_caption = perturb xyz coordinates of 0,0,0 atoms by this much after gromacs' pdb2gmx
  .help = This exists for troubleshooting
remove_metals    = True
  .type          = bool
  .short_caption = If true, remove MG and ZN during cleaning before pdb2gmx
debug = False
  .type = bool
  .expert_level=3
  .short_caption = debug output
  .help = Debug output
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
}
"""
master_params = master_params_str
master_phil = phil.parse(master_params_str, process_includes=True)
# This sentence works before main function

'''
def end_regression(starting_dir,write_this):
  print "end regression (both for each step and all steps)"
  os.chdir (starting_dir)
  print write_this
  
  if (os.path.isfile("cc_record_full") == True):
    rm_command_string = "rm cryo_fit* cc_record_full"
  else:
    rm_command_string = "rm cryo_fit*"
  libtbx.easy_run.fully_buffered(rm_command_string)
  
  exit(1)
########## end of end_regression
'''

def get_release_tag():
  release_tag = os.environ.get("PHENIX_RELEASE_TAG", None)
  return release_tag

def get_version():
    version = os.environ.get("PHENIX_VERSION", None)
    if (version is None):
      tag_file = libtbx.env.under_dist("libtbx", "../TAG")
      if (os.path.isfile(tag_file)):
        try: version = open(tag_file).read().strip()
        except KeyboardInterrupt: raise
        except: pass
    return version

def print_author():
  version = get_version()
  release_tag = get_release_tag()
  print """\
 %s
  cryo_fit %s 
    - Doo Nam Kim, Nigel Moriarty, Serdal Kirmizialtin, Billy Poon, Karissa Sanbonmatsu
 %s""" % ("-"*78, version, "-"*78)
# end of print_author()

def show_header(title):
  print "\n"
  print '#'*95
  number_of_remaining_sharp = 95 - len(title)
  put_this_number_of_sharp = int(int(number_of_remaining_sharp)/2)
  print '#'*(put_this_number_of_sharp-1) + " " + title + " " + '#'*(put_this_number_of_sharp-1)
  print '#'*95
# end of show_header function

def validate_params(params): # validation for GUI
  # check if file type is OK
  #if (params.cryo_fit.Input.cryo_fit_path is None):
  #  raise Sorry("cryo_fit_path should be given, please install gromacs_cryo_fit")
  if (params.cryo_fit.Input.map_file_name is None):
    raise Sorry("Map file should be given")
  if (params.cryo_fit.Input.model_file_name is None):
    raise Sorry("Model file should be given")
  
  file_reader.any_file(
    file_name = params.cryo_fit.Input.model_file_name).check_file_type(expected_type = 'pdb')

  #file_reader.any_file(
  #  file_name = params.cryo_fit.map_file_name).check_file_type(expected_type = 'map')
  # Doonam commented this for now since it resulted in "AttributeError: 'scope_extract' object has no attribute 'map_file_name'"

  print "\tvalidate_params pass"
  return True
############### end of validate_params function

  
def step_1(logfile, command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, \
           force_field, ignh, missing, remove_metals, cryo_fit_path, *args):
  show_header("Step 1: Make gro and topology file by regular gromacs")
  remake_and_move_to_this_folder(starting_dir, "steps/1_make_gro")
  
  cw_dir = os.getcwd()
  print "\tCurrent working directory: %s" % cw_dir
  
  # can't copy "3f2q-FMN riboswitch-fit.pdb"
  #cp_command_string = "cp " + model_file_with_pathways + " ."
  #libtbx.easy_run.fully_buffered(cp_command_string)
  
  # copied "3f2q-FMN riboswitch-fit.pdb" well
  shutil.copy(model_file_with_pathways, ".")

  start = time.time()
  cp_command_string = "cp " + command_path + "files_for_steps/1_make_gro/1_before_pdb2gmx_prepare_pdb.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  if (model_file_without_pathways.find("_cleaned_for_gromacs") == -1):
    run_this = "python 1_before_pdb2gmx_prepare_pdb.py " + model_file_without_pathways + " 0 0 0 " + \
               str(remove_metals)
    print "\tcommand: ", run_this
    f_out = open('log.step_1_1_before_pdb2gmx_prepare_pdb', 'wt')
    write_this_input_command = run_this + "\n"
    f_out.write(write_this_input_command)
    f_out.close()
    libtbx.easy_run.fully_buffered(run_this)

  pdb_file_is_cleand = False
  for check_this_file in glob.glob("*_cleaned_for_gromacs.pdb"): 
   pdb_file_is_cleand = True
    
  if (pdb_file_is_cleand != True):
    print "pdb file cleaning is not done, exit now"
    print "\nPlease email phenixbb@phenix-online.org or doonam@lanl.gov for any feature request/help."
    exit(1)
    
  cp_command_string = "cp " + command_path + "files_for_steps/1_make_gro/2_runme_make_gro.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  os.chdir (starting_dir)
  
  cw_dir = os.getcwd()
  print "\tCurrent working directory: %s" % cw_dir
  
  number_of_atoms_in_input_pdb = know_number_of_atoms_in_input_pdb(model_file_with_pathways)  
  if (number_of_atoms_in_input_pdb < 7000): # GTPase_activation_center for development.
    print "\tApproximately, for this number of atoms, one 3.1 GHz Intel Core i7 took 7 seconds to make a gro file.\n"
  elif (number_of_atoms_in_input_pdb < 20000): # nucleosome has 14k atoms (pdb), 25k atoms (gro)
    print "\tApproximately, for this number of atoms, one 3.1 GHz Intel Core i7 took 4 minutes to make a gro file.\n"
  elif (number_of_atoms_in_input_pdb < 50000): # beta-galactosidase has 32k atoms (pdb), 64k atoms (gro)
    print "\tApproximately, for this number of atoms, one 3.1 GHz Intel Core i7 took 7 minutes to make a gro file.\n"
  else: # ribosome has 223k atoms (lowres_SPLICE.pdb)
    print "\tApproximately, for this number of atoms, one 3.1 GHz Intel Core i7 took 2 hours to make a gro file.\n"
    
  new_path = starting_dir + "/steps/1_make_gro"
  os.chdir( new_path )
  
  command_script = "python 2_runme_make_gro.py " + str(command_path) + " " + force_field + " " + \
            str(ignh) + " " + str(missing) + " " + str(cryo_fit_path)
  # there is only 1 pdb file in this folder, so it is ok not to provide pdb arguments
  
  print "\tcommand: ", command_script
  libtbx.easy_run.call(command_script)
  end = time.time()
  
  this_step_was_successfully_ran = "failed" # just an initial value
  for check_this_file in glob.glob("*_by_pdb2gmx.gro"): # there will be only one *_by_pdb2gmx.gro file
    this_step_was_successfully_ran = check_whether_the_step_was_successfully_ran("Step 1", check_this_file, logfile)
  if (this_step_was_successfully_ran == "failed"):
    logfile.write("Step 1 didn't run successfully")
    
    print "Step 1 didn't run successfully"
    print "\nUser's command was"
    f_in = open('../../cryo_fit.input_command')
    for line in f_in:
      print line
    #check_whether_mdrun_is_accessible() is ran above, confirmed that there is no reason to suggest install gromacs_cryo_fit here again
    print "\nphenix.cryo_fit alone without any arguments introduces full options."
    
    print_this = '''\nIf a user sees a message like\n\
    "Fatal error: \n\
    Atom xx in residue xx xxx was not found in rtp entry xx with xx atoms \n\
    while sorting atoms."
    
    or \n\
    "Fatal error: \n\
    Residue xxx not found in residue topology database"
    
    Solution if these residue/atoms are important:\n
    \tFix wrong names of atoms/residues. Running real_space_refine via phenix GUI will show which atoms need to be fixed.
    \tIf gromacs amber03 force field doesn't have parameters for these residue/atoms, you may need to add appropriate parameters.
    \tIf you added parameters, please email me (doonam@lanl.gov), I want to recognize your contribution publicly.
    \tMost MD simulation force fields do not support all kinds of rare residue/atoms.
    \tcryo_fit2 is under development to address this issue using phenix.eLBOW
    
    Solution if these residue/atoms are not important:\n
    \tRemove these \"wrong\" atoms/residues from user's input pdb file. Run cryo_fit again.
    
    Solution if user's input pdb file is big:\n
    \tProbably cryo_fit will change conformation just minimally, I would extract out these \"wrong\" atoms/residues from user's input pdb file, then add these extracted lines to cryo_fitted file later.
    '''
    print print_this
    logfile.write(print_this)
    logfile.close()
    
    print "\nEmail phenixbb@phenix-online.org or doonam@lanl.gov for any feature request/help."
    exit(1)
  print "Step 1", (show_time(start, end))
  
  os.chdir (starting_dir)
  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_GAC.pdb") \
   or (model_file_without_pathways == "regression_Adenylate.pdb") \
   or (model_file_without_pathways == "regression_tRNA_EFTU_within_10.pdb")):
    this_is_test_for_each_step = True
  
  if (os.path.isfile("steps/1_make_gro/prefix_of_chain_ID_removed") == True):
      write_this = "The 4th character of residue name (prefix_of_chain ID) is removed. Please see https://www.phenix-online.org/documentation/faqs/cryo_fit_FAQ.html#how-can-i-use-double-digit-character-id-pdb-file or email doonam.kim@gmail.com\n\n"
      print write_this
      logfile.write(write_this)
      
  return this_is_test_for_each_step
############################################## end of step_1 function


def step_2(logfile, command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, \
           force_field, perturb_xyz_by, remove_metals):
  show_header("Step 2: Clean gro file to be compatible for amber03 forcefield")
  os.chdir (starting_dir)
  remake_and_move_to_this_folder(starting_dir, "steps/2_clean_gro")

  cp_command_string = ''

  this_is_test_for_each_step = False # default
  #if ((model_file_without_pathways == "regression_GAC.pdb") or (model_file_without_pathways == "regression_Adenylate.pdb")):
  if ((model_file_without_pathways == "regression_GAC.pdb") \
   or (model_file_without_pathways == "regression_Adenylate.pdb") \
   or (model_file_without_pathways == "regression_tRNA_EFTU_within_10.pdb")):
    this_is_test_for_each_step = True
    cp_command_string = "cp ../../data/input_for_step_2/*_cleaned_for_gromacs_by_pdb2gmx.gro ."
  else: # regular running and emd (both wo_restart and w_restart)
    cp_command_string = "cp ../1_make_gro/*.gro ."
  print "cp_command_string:",cp_command_string
  
  libtbx.easy_run.fully_buffered(cp_command_string) #copy step_1 output

  start_time_renaming = time.time()
  print "\nStep 2: Add C prefix to terminal amino acid/nucleic acid for minimization by gromacs"
  #command_script = "cp " + command_path + "steps/2_clean_gro/*.py ."
  command_script = "cp " + command_path + "files_for_steps/2_clean_gro/*.py ."
  libtbx.easy_run.fully_buffered(command_script)

  command_script = "python 1_rename_term_res_to_Cres.py " # there will be only 1 gro file, so it is ok
  # this runs both 2_rename_term_res_to_Cres_by_resnum.py and 3_rename_term_res_to_Cres_by_oc.py
  print "\tcommand: ", command_script
  libtbx.easy_run.fully_buffered(command_script)

  for this_file in glob.glob("*_c_term_renamed_by_resnum_oc.gro"): # there will be only one file like this
    command_string = "python 4_slightlty_change_xyz_for_no_more_000.py " + this_file + " " + str(perturb_xyz_by)
    print "\tcommand: ", command_string 
    libtbx.easy_run.call(command_string)
  
  if (this_is_test_for_each_step == True):
    return True
    
  this_step_was_successfully_ran = "failed" # just an initial value
  for check_this_file in glob.glob("*.gro"): # there will be only "will_be_minimized_cleaned.gro"
    this_step_was_successfully_ran = check_whether_the_step_was_successfully_ran("Step 2", check_this_file, logfile)

  if (this_step_was_successfully_ran == "failed"):
    color_print (("Step 2 didn't run successfully"), 'red')
    exit(1)
  
  for check_this_file in glob.glob("*.gro"): # there will be only one file like this # to work in Karissa's old MacOS
    command_string = "mv " + check_this_file + " will_be_minimized_cleaned.gro"
    libtbx.easy_run.call(command_string)
    
  end_time_renaming = time.time()
  print "Step 2", (show_time(start_time_renaming, end_time_renaming))
  
  return False # this is not a test for each step
############################### end of step_2 (clean gro) function


def step_3(logfile, command_path, starting_dir, ns_type, restraint_algorithm_minimization, number_of_steps_for_minimization, \
           time_step_for_minimization, model_file_without_pathways, devel, cryo_fit_path):
  show_header("Step 3: Make a tpr file for minimization")
  os.chdir (starting_dir)

  remake_this_folder("steps/3_make_tpr_to_minimize")
  remake_and_move_to_this_folder(starting_dir, "steps/3_make_tpr_to_minimize")

  cp_command_script = "cp " + command_path + "files_for_steps/3_make_tpr_to_minimize/minimization_template.mdp ."
  libtbx.easy_run.fully_buffered(cp_command_script)

  if (("regression_" in model_file_without_pathways) or (devel == True)):
    number_of_steps_for_minimization = 5
  
  print "\tBe number_of_steps_for_minimization as ", number_of_steps_for_minimization
  with open("minimization_template.mdp", "rt") as fin:
    with open("minimization.mdp", "wt") as fout:
      for line in fin:
        splited = line.split()
        if splited[0] == "nsteps":
          new_line = "nsteps          = " + str(number_of_steps_for_minimization) + " ; Maximum number of minimization steps to perform\n"
          fout.write(new_line)
        elif splited[0] == "ns_type":
          new_line = "ns_type         = " + str(ns_type) + " ; Method to determine neighbor list (simple, grid)\n"
          fout.write(new_line)
        else:
          fout.write(line)
      print "\ttime_step_for_minimization:", time_step_for_minimization
      if time_step_for_minimization != 0.001:
        print "\ttime_step_for_minimization != 0.001"
        new_line = "\ndt = " + str(time_step_for_minimization) + "\n"
        fout.write(new_line)
      if str(restraint_algorithm_minimization) == "None" or str(restraint_algorithm_minimization) == "none":
        print "\trestraint_algorithm_minimization = none"
        new_line = "\nrestraint-algorithm: none\n"
        fout.write(new_line)
    fout.close()
  fin.close()
  
  cp_command_script = "cp " + command_path + "files_for_steps/3_make_tpr_to_minimize/runme_make_tpr.py ."
  libtbx.easy_run.fully_buffered(cp_command_script)
  
  cp1_command_string = '' #initialization
  
  this_is_test_for_each_step = False # default
  
  if ("regression_pdb5khe" in model_file_without_pathways):
    cp1_command_string = "cp ../2_clean_gro/*.gro . "  
    cp2_command_string = "cp ../1_make_gro/*.top . "
    libtbx.easy_run.fully_buffered(cp2_command_string)
  elif ("regression_" in model_file_without_pathways):
    this_is_test_for_each_step = True
    cp1_command_string = "cp ../../data/input_for_step_3/* ."
  else: # regular running
    if str(restraint_algorithm_minimization) != "none_default":
      cp1_command_string = "cp ../2_clean_gro/*.gro . "
    else:
      cp1_command_string = "mv ../../minimized_c_term_renamed_by_resnum_oc.gro . "
    cp2_command_string = "cp ../1_make_gro/*.top . "
    libtbx.easy_run.fully_buffered(cp2_command_string)
  libtbx.easy_run.fully_buffered(cp1_command_string)
  
  if (this_is_test_for_each_step == False):
    if str(restraint_algorithm_minimization) != "none_default":
      cp1_command_string = "cp ../2_clean_gro/*.gro . "  
      cp2_command_string = "cp ../1_make_gro/*.top . "
      libtbx.easy_run.fully_buffered(cp2_command_string)
    else:
      cp1_command_string = "mv ../../minimized_c_term_renamed_by_resnum_oc.gro . "
  else: # this_is_test_for_each_step = True
    cp1_command_string = "cp ../../data/input_for_step_3/* ."

  command_string = "python runme_make_tpr.py " + str(cryo_fit_path)
  print "\tcommand: ", command_string
  start = time.time()
  libtbx.easy_run.call(command_string)
  end = time.time()

  check_whether_the_step_3_was_successfully_ran(logfile, "to_minimize.tpr")
  print "Step 3", (show_time(start, end))
  os.chdir( starting_dir )
  if (this_is_test_for_each_step == True):
    return True
  return False # this is not a test for each step
############################ end of step_3 (prepare minimization) function


def step_4(logfile, command_path, starting_dir, ns_type, number_of_available_cores, \
           number_of_cores_to_use, model_file_without_pathways, cryo_fit_path):
  show_header("Step 4: Minimize a gro file (to prevent \"blowup\" during Molecular Dynamics Simulation)")
  os.chdir (starting_dir)
  
  print "\nStep 4-1: Minimization itself"
  remake_and_move_to_this_folder(starting_dir, "steps/4_minimize")

  cp_command_script = "cp " + command_path + "files_for_steps/4_minimize/runme_minimize.py ."
  libtbx.easy_run.fully_buffered(cp_command_script)
  
  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_GAC.pdb") or \
      (model_file_without_pathways == "regression_Adenylate.pdb") or \
      (model_file_without_pathways == "regression_tRNA_EFTU_within_10.pdb")):
    cp_command_string = "cp ../../data/input_for_step_4/* ."
    this_is_test_for_each_step = True
  else:
    cp_command_string = "cp ../3_make_tpr_to_minimize/to_minimize.tpr ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  # when there are both mpi and thread cryo_fit exist, thread cryo_fit was used in commandline mode
  command_string = "python runme_minimize.py to_minimize.tpr " + str(command_path) + " " + \
                str(ns_type) + " " + str(number_of_available_cores) + " " + str(2) + " " + str(cryo_fit_path)
              # set number_of_cores_to_use = 2 to minimize a possibility of having cell size error
  print "\tcommand: ", command_string
  print "\n\tA user can check progress of minimization at ", starting_dir + "/steps/4_minimize\n"
  start = time.time()
  
  libtbx.easy_run.call(command_string)
  
  ''' # seems not working
  f_in = open('log.step_4_1_minimization_real_command')
  
  # progress is shown to both commandline & GUI
  # I thank https://stackoverflow.com/questions/42553481/check-on-the-stdout-of-a-running-subprocess-in-python
  
  for line in f_in:
    splited = line.split('\'')
    from subprocess import Popen, PIPE, STDOUT
    double_splited = splited[1].split()
    
    p = Popen([double_splited[0], double_splited[1], double_splited[2], \
              double_splited[3], double_splited[4], double_splited[5], \
              double_splited[6], double_splited[7], double_splited[8], \
              double_splited[9]], stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    for line in p.stdout:
      print(line)

    #libtbx.easy_run.call(splited[1]) # progress was not shown to GUI
  
  #for line in f_in:
  #  splited = line.split('\'')
  #  os.system(splited[1]) # progress was shown to commandline
  
  f_in.close()
  '''
  end = time.time()
  
  final_gro_file_name = ''
  for gro_file_name in glob.glob("*.gro"): # there will be only one file like this
    final_gro_file_name = gro_file_name
    
  returned = check_whether_the_step_was_successfully_ran("Step 4-1", final_gro_file_name, logfile)
  if returned == "failed":
    #bool_enable_mpi = know_output_bool_enable_mpi_by_ls()
    #if bool_enable_mpi == True:
      color_print ("\n<Case 1> When Doonam encountered this error message", 'red')
      color_print ("\t\"[simplednlanlgov.local:14805] [[58528,0],0] usock_peer_recv_connect_ack: received different version from [[58528,1],0]: 2.1.1 instead of 2.1.0\"", 'red')
      color_print ("\t\"-------------------------------------------------------\"", 'red')
      color_print ("\t\"Primary job terminated normally, but 1 process returned\"", 'red')
      color_print ("\t\"a non-zero exit code.. Per user-direction, the job has been aborted.\"", 'red')
      color_print ("\t\"-------------------------------------------------------\"", 'red')
      color_print ("\t\"[simplednlanlgov.local:14805] [[58528,0],0] usock_peer_recv_connect_ack: received different version from [[58528,1],1]: 2.1.1 instead of 2.1.0\"", 'red')
      color_print ("\t\"--------------------------------------------------------------------------\"", 'red')
      color_print ("\t\"mpirun detected that one or more processes exited with non-zero status, thus causing\"", 'red')
      color_print ("\t\"the job to be terminated. The first process to do so was:\"", 'red')
      color_print ("he used /usr/local/bin/mpirun.", 'red')
      color_print ("Using /Users/doonam/bin/openmpi-2.1.1/bin/mpirun (after reinstalling openmpi) solved the problem.", 'green')
      color_print ("Using /Users/doonam/EMAN2/bin/mpirun solved the problem as well.", 'green')
  
      color_print ("\n<Case 2> One user encountered this error message,", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] Signal: Bus error: 10 (10)\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] Signal code:  (2)\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] Failing at address: 0x104043018\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 0] 2   libsystem_platform.dylib            0x00007fff874045aa _sigtramp + 26\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 1] 3   ???                                 0x00007f0300000000 0x0 + 139650861629440\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 2] 4   libopen-pal.0.dylib                 0x00000001040b9c42 mca_base_components_open + 1698\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 3] 5   libopen-pal.0.dylib                 0x00000001040c97f6 opal_memory_base_open + 38\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 4] 6   libopen-pal.0.dylib                 0x00000001040aaf7c opal_init + 76\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 5] 7   libopen-rte.0.dylib                 0x00000001040487d0 orte_init + 32\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 6] 8   mpirun                              0x0000000104036650 orterun + 432\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 7] 9   mpirun                              0x0000000104036402 main + 34\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] [ 8] 10  libdyld.dylib                       0x00007fff8dddb5fd start + 1\"", 'red')
      color_print ("\t\"[Karissas-MacBook-Pro-2:11236] *** End of error message ***\"", 'red')
      color_print ("when she used an old openmpi (version 1.2.8)", 'red')
      color_print ("To solve a problem with mpirun, reinstalling cryo_fit WITHOUT mpi mode is recommended.", 'green')
      color_print ("Otherwise, a user may reinstall openmpi by python \
                   <user_phenix>/modules/cryo_fit/command_line/install_openmpi.py openmpi-2.1.1.tar.gz", 'green')
      color_print ("and use newly installed mpirun by setting PATH.", 'green')
      exit(1)
  print "Step 4-1", (show_time(start, end))
  
  print "\nStep 4-2: Add C prefix to terminal amino acids to minimized.gro for grompp by gromacs"
  cp_command_string = "cp " + command_path + "files_for_steps/2_clean_gro/*_rename_term_res_to_Cres*.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  command_string = "python 1_rename_term_res_to_Cres.py "
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)

  check_whether_the_step_was_successfully_ran("Step 4-2", "minimized_c_term_renamed_by_resnum_oc.gro", logfile)
  os.chdir( starting_dir )
  if (this_is_test_for_each_step == True):
    return True
  return False # this is not a test for each step
############################ end of step_4 (minimization) function
    

def step_5(logfile, command_path, starting_dir, model_file_without_pathways, cryo_fit_path):
  show_header("Step 5: Make contact potential (restraints) and topology file with it")
  remake_and_move_to_this_folder(starting_dir, "steps/5_make_restraints")
  
  start = time.time()
  cp_command_string = "cp " + command_path + "files_for_steps/5_make_restraints/runme_make_contact_potential.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_Adenylate.pdb") or
      (model_file_without_pathways == "regression_GAC.pdb") or
      (model_file_without_pathways == "regression_tRNA_EFTU_within_10.pdb")):
    cp_command_string = "cp ../../data/input_for_step_5/* ."
    this_is_test_for_each_step = True
  else:
    cp_command_string = "cp ../4_minimize/*.gro ."
  libtbx.easy_run.fully_buffered(cp_command_string)
    
  command_string = "python runme_make_contact_potential.py *.gro " + str(command_path) + " " + str(cryo_fit_path)
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)

  check_whether_the_step_was_successfully_ran("Step 5", "minimized_c_term_renamed_by_resnum_oc_including_disre2_itp.top", logfile)
  end = time.time()
  print "Step 5", (show_time(start, end))
  
  #color_print ((show_time("Step 5", start, end)), 'green')
  # [keep] looks as "[32mStep 5 finished in 10.66 seconds (wallclock).[0m" in GUI
  
  os.chdir( starting_dir )
  if (this_is_test_for_each_step == True):
    return True
  return False # this is not a test for each step
########################### end of step_5 (make restraints) function


def step_6(logfile, command_path, starting_dir, model_file_without_pathways):
  show_header("Step 6: Make all charges of atoms be 0")
  remake_and_move_to_this_folder(starting_dir, "steps/6_make_0_charge")

  cp_command_string = "cp " + command_path + "files_for_steps/6_make_0_charge/changetop.awk ."
  libtbx.easy_run.fully_buffered(cp_command_string)
      
  cp_command_string = "cp " + command_path + "files_for_steps/6_make_0_charge/runme_make_0_charge.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_Adenylate.pdb") or
      (model_file_without_pathways == "regression_GAC.pdb") or
      (model_file_without_pathways == "regression_tRNA_EFTU_within_10.pdb")):
    cp_command_string = "cp ../../data/input_for_step_6/* ."
    this_is_test_for_each_step = True
  else:
    cp_command_string = "cp ../5_make_restraints/*including_disre2_itp.top ." ## In normal case, there will be minimized_c_term_renamed_by_resnum_oc_including_disre2_itp.top
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  command_string = "python runme_make_0_charge.py *.top"
  print "\tcommand: ", command_string
  start = time.time()
  libtbx.easy_run.fully_buffered(command_string)
  end = time.time()

  for check_this_file in glob.glob("*_0_charge.top"): # there will be only one file like this
    check_whether_the_step_was_successfully_ran("Step 6", check_this_file, logfile)
    
  print "Step 6", (show_time(start, end))
  os.chdir( starting_dir )
  if (this_is_test_for_each_step == True):
    return True
  return False # this is not a test for each step
######################################## end of step_6 (neutralize charge) function


def step_7(logfile, command_path, starting_dir, number_of_steps_for_cryo_fit, emweight_multiply_by, \
           emsteps, emwritefrequency, lincs_order, nstxtcout, time_step_for_cryo_fit, \
           model_file_without_pathways, cryo_fit_path, many_step_____n__dot_pdb):
  show_header("Step 7 : Make a tpr file for cryo_fit")
  remake_and_move_to_this_folder(starting_dir, "steps/7_make_tpr_with_disre2")
  
  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_Adenylate.pdb") or
      (model_file_without_pathways == "regression_GAC.pdb") or
      (model_file_without_pathways == "regression_tRNA_EFTU_within_10.pdb")):
    cp1_command_string = "cp ../../data/input_for_step_7/* ."
    libtbx.easy_run.fully_buffered(cp1_command_string)
    this_is_test_for_each_step = True
  else: # regular running or emd regression
    cp1_command_string = "cp ../5_make_restraints/*.gro ." # there will be minimized_c_term_renamed_by_resnum_oc.gro only
    libtbx.easy_run.fully_buffered(cp1_command_string)
    cp2_command_string = "cp ../5_make_restraints/disre2.itp ."
    libtbx.easy_run.fully_buffered(cp2_command_string)
    cp3_command_string = "cp ../6_make_0_charge/*0_charge.top ." # there is only one *0_charge.top file
    libtbx.easy_run.fully_buffered(cp3_command_string)
  
  print "\tBe number_of_steps_for_cryo_fit as ", number_of_steps_for_cryo_fit
    
  fout = open("for_cryo_fit.mdp", "wt")
  fin = ''
  if (many_step_____n__dot_pdb == False):
    cp_command_string = "cp " + command_path + "files_for_steps/7_make_tpr_with_disre2/template_for_cryo_fit.mdp ."
    libtbx.easy_run.fully_buffered(cp_command_string)
    fin = open("template_for_cryo_fit.mdp", "rt")
  else:
    cp_command_string = "cp " + command_path + "files_for_steps/7_make_tpr_with_disre2/template_for_cryo_fit_many_step_____n__dot_pdb.mdp ."
    libtbx.easy_run.fully_buffered(cp_command_string)
    fin = open("template_for_cryo_fit_many_step_____n__dot_pdb.mdp", "rt")
  write_for_cryo_fit_mdp(fout, fin, emsteps, time_step_for_cryo_fit, number_of_steps_for_cryo_fit, \
                         emweight_multiply_by, emwritefrequency, lincs_order, nstxtcout)
  
  cp_command_string = "cp " + command_path + "files_for_steps/7_make_tpr_with_disre2/runme_make_tpr_with_disre2.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  start_make_tpr = time.time()
  command_string = "python runme_make_tpr_with_disre2.py " + str(command_path) + " " + str(cryo_fit_path)
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string) 
  end_make_tpr = time.time()
  
  
  if (this_is_test_for_each_step == True):
    return True
  returned = check_whether_the_step_was_successfully_ran("Step 7", "for_cryo_fit.tpr", logfile)
  print "Step 7", (show_time(start_make_tpr, end_make_tpr))
  if (returned == "failed"):
    exit(1)
  os.chdir( starting_dir )
  return False # this is not a test for each step
########################### end of step_7 (make tpr for cryo_fit) function        
   
           
def step_8(logfile, command_path, starting_dir, number_of_available_cores, number_of_cores_to_use, \
       map_file_with_pathways, no_rerun, devel, restart_w_longer_steps, re_run_with_higher_map_weight, \
       model_file_without_pathways, cryo_fit_path):
  show_header("Step 8: Run cryo_fit")
  print "\tmap_file_with_pathways:",map_file_with_pathways
  
  remake_and_move_to_this_folder(starting_dir, "steps/8_cryo_fit")
  
  command_string = "cp " + command_path + "files_for_steps/8_cryo_fit/* ."
  libtbx.easy_run.fully_buffered(command_string)
  
  if (str(restart_w_longer_steps) == "True"):
    command_string = "cp ../../state.cpt . "
    libtbx.easy_run.fully_buffered(command_string)
  
  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_GAC.pdb") or (model_file_without_pathways == "regression_Adenylate.pdb")):
    this_is_test_for_each_step = True  
    
  command_string = "python runme_cryo_fit.py " + str(command_path) + " " + str(number_of_available_cores) \
              + " " + number_of_cores_to_use + " " + map_file_with_pathways + " " + str(starting_dir) \
              + " " + str(this_is_test_for_each_step) + " " + str(restart_w_longer_steps) + " " + str(cryo_fit_path)
  print "\n\tcommand: ", command_string
  print "\n\tA user can check progress of step_8 at ", starting_dir + "/steps/8_cryo_fit\n"
  time_start_cryo_fit = time.time()
  libtbx.easy_run.call(command_string)
  time_end_cryo_fit = time.time()
  
  # progress is shown to monitor in GUI (but not super-fast), but after run is done
  '''
  f_in = open('md.log')
  for line in f_in:
    print(line)
  f_in.close()
  '''
  
  '''
  # 04/18/2018, this doesn't show progress to GUI, but keep for now
  f_in = open('log.step_8_cryo_fit_used_command') 
  for line in f_in:
    # progress is shown to monitor in GUI (but not super-fast)
    from subprocess import Popen, PIPE, STDOUT
    splited = line.split()
    p_cryo_fit = Popen([splited[0], splited[1], splited[2], \
              splited[3], splited[4], splited[5], \
              splited[6], splited[7], splited[8], \
              splited[9], splited[10], splited[11], \
              splited[12], splited[13], splited[14]], \
              stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    for line_current in p_cryo_fit.stdout:
      print(line_current)
    #os.system(line) # progress is shown to commandline, not shown to GUI
    #libtbx.easy_run.call(line) # progress is shown to commandline, not shown to GUI
  f_in.close()
  '''
  
  command_string = "cat md.log | grep correlation > cc_record"
  libtbx.easy_run.fully_buffered(command_string)
  
  command_string = "cat md.log | grep correlation >> ../../cc_record_full"
  libtbx.easy_run.fully_buffered(command_string)
  
  if (this_is_test_for_each_step == True):
    return this_is_test_for_each_step
  
  returned = check_whether_the_step_was_successfully_ran("Step 8", "cc_record", logfile)
  
  if (returned != "success"):
    write_this = "Step 8 (Run cryo_fit) didn't run successfully"
    print write_this
    logfile.write(write_this)
    
    write_this = "\nPlease see https://www.phenix-online.org/documentation/faqs/cryo_fit_FAQ.html\n"
    print write_this
    logfile.write(write_this)
    
    if (returned == "failed_with_nan_in_cc"):
      return "failed_with_nan_in_cc"
    
    searched = search_charge_in_md_log()
    if searched == 0: # no "charge group... " message in md.log
      return "failed"
    else:
      return "re_run_w_smaller_MD_time_step"

  '''
  if ((str(restart_w_longer_steps) == "False") and (str(re_run_with_higher_map_weight) == "False")):
    user_s_cc = get_users_cc("cc_record")
    user_s_cc = round(float(user_s_cc),3)
    write_this = "\nA user's provided atomic model had " + str(user_s_cc) + " cc\n\n"
    logfile.write(write_this)
  '''

  if (no_rerun == False):
    this_is_test = False
    if "regression_" in model_file_without_pathways:
      this_is_test = True
    cc_has_been_increased = check_whether_cc_has_been_increased(logfile, "cc_record", this_is_test)
    print "\tVerdict of cc_has_been_increased function with all cc evaluations:", cc_has_been_increased

    if (devel == True):
        no_rerun = True
    if (no_rerun == False):
      if cc_has_been_increased == "increased":
        return "re_run_with_longer_steps"
      elif cc_has_been_increased == "re_run_with_higher_map_weight":
        return "re_run_with_higher_map_weight"
      else: # cc_has_been_increased = "saturated"  #### cc_has_been_increased = False
        write_this = "\ncc has been saturated, so cryo_fit will go to the next step (e.g. final_output_arrange_step)"
        print write_this
        logfile.write(write_this)

  f_out = open('log.step_8', 'at+')
  write_this_time = show_time(time_start_cryo_fit, time_end_cryo_fit)
  write_this_time = "\n\nStep 8" + write_this_time + "\n"
  f_out.write(write_this_time)
  f_out.close()

  renumber_cc_record_full("../../cc_record_full")
  os.remove("../../cc_record_full")

  f_in = open('../../cc_record_full_renumbered')
  cc_record = list()
  for line in f_in:
    splited = line.split()
    step = splited[1]
    cc = splited[4]
    cc_record.append((float(step), float(cc)))
  f_in.close()
  
  results = dict()
  results['cc_record'] = cc_record # results should have ['cc_record'] only, if it has ['cc_has_been_increased'] as well, GUI will error
  
  print "\nStep 8", (show_time(time_start_cryo_fit, time_end_cryo_fit))
  
  if (os.path.isfile("../../state.cpt") == True):
    os.remove("../../state.cpt")
  
  os.chdir( starting_dir )
  
  return results
######################################## end of step_8 (cryo_fit itself) function


def step_final(logfile, command_path, starting_dir, model_file_without_pathways, cryo_fit_path, no_rerun):
  os.chdir( starting_dir )
  time_start = time.time()
  show_header("Step 9 (final): Arrange output")
  remake_and_move_to_this_folder(starting_dir, "output")
  
  this_is_test_for_each_step = False # just initial value assignment
  
  cp_command_string = "cp " + command_path + "files_for_steps/9_after_cryo_fit/*.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  logfile.close() # to write user's cc for now
 
  this_is_test_for_each_step = False # default
  if ((model_file_without_pathways == "regression_GAC.pdb") or (model_file_without_pathways == "regression_Adenylate.pdb")):
    cp_command_string = "cp ../data/input_for_step_final/* ."
    libtbx.easy_run.fully_buffered(cp_command_string)
    this_is_test_for_each_step = True
  else: # regular running or emd regression
    cp_command_string = "mv ../cc_record_full_renumbered ."
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    cp_command_string = "cp ../steps/8_cryo_fit/cc_record ." # needed for extract_3_highest_cc_gro
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    cp_command_string = "cp ../steps/8_cryo_fit/for_cryo_fit.tpr ../steps/8_cryo_fit/traj.xtc ."
    libtbx.easy_run.fully_buffered(cp_command_string)

  print "\n\tExtract .gro files from the 3 highest cc values."
  if os.path.isfile("extract_3_highest_cc_gro.py") == False:
    print "extract_3_highest_cc_gro.py is not found, please email doonam@lanl.gov. Exit now"
    exit(1)

  log_file_name = "../cryo_fit.overall_log"
  logfile = open(log_file_name, "a+") # append
  
  command_string = "python extract_3_highest_cc_gro.py " + str(this_is_test_for_each_step) + " " + str(cryo_fit_path) + " " + str(no_rerun)
  write_this = command_string + "\n\n"
  logfile.write(write_this)
  print "\tcommand: ", write_this
  libtbx.easy_run.call(command_string)

  if (os.path.isfile("extract_gro_failed.txt") == True):
    return "failed"
  
  for extracted_gro in glob.glob("*.gro"):
    returned_file_size = file_size(extracted_gro)
    if (returned_file_size == 0):
        exit(1)
        
  print "\n\tConvert .gro -> .pdb"
  print "\t\t(.gro file is for Chimera/Gromacs/Pymol/VMD)"
  print "\t\t(.pdb file is for Chimera/ChimeraX/Pymol/VMD)"
  
  for extracted_gro in glob.glob("*.gro"): # just deals .gro files in alphabetical order not in cc order
    command_string = cryo_fit_path + "editconf -f " + extracted_gro + " -o " + extracted_gro[:-4] + ".pdb"
    write_this = command_string + "\n"
    logfile.write(write_this)
    print "\t\tcommand: ", write_this
    libtbx.easy_run.fully_buffered(command_string)
  
  # deal with trajectory files
  make_trajectory_gro(cryo_fit_path)
    
  os.mkdir("trajectory")
  run_this = "mv for_cryo_fit.tpr trajectory.gro traj.xtc trajectory"
  print "\n\tMove trajectory files into trajectory folder "
  print "\t\tcommand: ", run_this
  libtbx.easy_run.call(run_this)
  
  pdb_file_with_original_chains = ''
  for pdb_with_original_chains in glob.glob("../steps/1_make_gro/*.pdb"):
    pdb_file_with_original_chains = pdb_with_original_chains
  
  # log_file_name = "../cryo_fit.overall_log"
  # logfile = open(log_file_name, "a+") # append
  
  if (this_is_test_for_each_step == False): # recover chain information
    write_this = "\n\tRecover chain information (since gromacs erased it).\n"
    logfile.write(write_this)
    print write_this
    
    for pdb in glob.glob("*.pdb"):
      command_string = "python recover_chain.py " + pdb_file_with_original_chains + " " + pdb # worked perfectly with GTPase_activation_center and Dieter's molecule
      print "\n\t\tcommand: ", command_string
      libtbx.easy_run.fully_buffered(command_string)
      
      chain_recovered = pdb[:-4] + "_chain_recovered.pdb"
      number_of_atoms_in_pdb_after_cryo_fit = know_number_of_atoms_in_input_pdb(pdb) # I need to use "pdb" not "pdb_file_with_original_chains"
      number_of_atoms_in_pdb_after_cryo_fit_chain_recovered = know_number_of_atoms_in_input_pdb(chain_recovered)
      run_this = ''
      write_this = '' 
      if (number_of_atoms_in_pdb_after_cryo_fit == number_of_atoms_in_pdb_after_cryo_fit_chain_recovered):
        print "\t\tnumber_of_atoms_in_pdb_after_cryo_fit = number_of_atoms_in_pdb_after_cryo_fit_chain_recovered"
        print "\t\tTherefore, chain_recovery is successful"
        write_this = "\tchain_recovery is successful\n"
        run_this = "rm " + pdb
      else:
        print "\t\tnumber_of_atoms_in_pdb_after_cryo_fit != number_of_atoms_in_pdb_after_cryo_fit_chain_recovered"
        print "\t\tTherefore, chain_recovery is not successful"
        write_this = "\tchain_recovery is not successful\n"
        run_this = "rm " + chain_recovered
      logfile.write(write_this)
      print "\t\trm: ", run_this
      libtbx.easy_run.call(run_this)
  
  print "\n\tClean pdb for realspace_refine and molprobity"
  for pdb in glob.glob("*.pdb"):
    run_this = "python clean_pdb_for_realspace_refine_molprobity.py " + pdb
    print "\t\tcommand: ", run_this
    libtbx.easy_run.call(run_this)
  
  returned = ''
  if (this_is_test_for_each_step == False):
    for pdb in glob.glob("*.pdb"):
      returned = check_whether_the_step_was_successfully_ran("Step final", pdb, logfile)
      break
    #returned = check_whether_the_step_was_successfully_ran("Step final", "cc_record_full_renumbered", logfile)
    # sometimes misclassified as success run although it failed to convert to pdb file
  
  print "\n\tOutputs are in \"output\" folder"
  print "  \t\tIf cryo_fit fitted better than a user provided atomistic model, a model with the highest cc value is cryo_fitted.pdb"
  print "  \t\tThis best fitted bio-molecule may not necessarily be the \"best\" atomic model depending on user's specific purposes."
  print "  \t\tTherefore, a user may use other slightly less fitted extracted_x_steps_x_ps.gro/pdb as well."
  print "\n\t\tTo draw a figure for cc,"
  print "  \t\t\tpython <phenix_path>/modules/cryo_fit/files_for_steps/9_after_cryo_fit/draw_cc/draw_cc.py cc_record"
  print "  \t\t\t(Phenix GUI shows this figure automatically)."
  print "  \t\t\t(If a user is using ssh linked linux, set DISPLAY to avoid \"Unable to access the X Display, is $DISPLAY set properly?\")"
  print "\n\t\tTo watch/record a trajectory movie, see \"https://www.phenix-online.org/documentation/tutorials/cryo_fit_movie.html\""
  
  if (this_is_test_for_each_step == True): # test for each individual steps like steps 1~final
    return this_is_test_for_each_step
    
  
  if (returned != "success"):
    write_this = "Step final (arrange output) didn't run successfully"
    print write_this
    logfile.write(write_this)
    
    logfile.close()
    return "failed"
    #exit(1) # even when step final failed, whole regression test appear as success
  
  for py in glob.glob("*.py"): # most users will not need *.py
    run_this = "rm " + py
    libtbx.easy_run.call(run_this)
    
  logfile.write("Step final (arrange output) is successfully ran\n")
  time_end = time.time()
  print "\nStep final", (show_time(time_start, time_end))
  
  os.chdir( starting_dir )
  return this_is_test_for_each_step
############################## end of step_final (arrange output) function


''' not used now, but keep
def step_9(command_path, starting_dir):
  show_header("Step 9: Show Correlation Coefficient")
  remake_and_move_to_this_folder(starting_dir, "steps/9_after_cryo_fit/draw_cc")
  
  command_string = "cp " + command_path + "steps/9_after_cryo_fit/draw_cc/draw_cc.py ."
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)
  
  command_string = "cp ../8_cryo_fit/md.log ."
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)
  
  cc_record = model_file_without_pathways[:-4] + "_fitted_to_" + map_file_without_pathways[:-4]
  command_string = "cat md.log | grep correlation > " + cc_record
  print "\n\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)
  
  returned = check_whether_the_step_was_successfully_ran("Step 9", cc_record)
  if returned == "failed":
    exit(1)
    
  command_string = "python draw_cc.py " + cc_record
  print "\n\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)
#end of step_9 (cc draw) function
'''



def run_cryo_fit(logfile, params, inputs):
  mdrun_path = check_whether_mdrun_is_accessible()
  cryo_fit_path = ''
  if (mdrun_path == False):
    if (params.cryo_fit.Input.cryo_fit_path != None):
      cryo_fit_path = params.cryo_fit.Input.cryo_fit_path
      cryo_fit_path = cryo_fit_path + "/" # for later steps
    else:
      long_message  = """
        cryo_fit can't find a gromacs executable (e.g. mdrun)
        
        If gromacs_cryo_fit is not installed, install it according to http://www.phenix-online.org/documentation/reference/cryo_fit.html
        
        If gromacs_cryo_fit is installed, source ~/.bash_profile or ~/.bashrc or open a new terminal so that cryo_fit path is included
        For example, if user's executables are installed at /Users/doonam/bin/cryo_fit/bin,
        add \"export PATH=\"/Users/doonam/bin/cryo_fit/bin\":$PATH" + " to ~/.bash_profile or ~/.bashrc and source it
        """
      print long_message
      print "\texit now"
      exit(1)
  else:
    cryo_fit_path = mdrun_path
  print "\tcryo_fit_path:",cryo_fit_path
  
  show_header("Step 0: Prepare to run cryo_fit")

  starting_dir = os.getcwd()
  print "\tCurrent working directory: %s" % starting_dir

  bool_step_1 = params.cryo_fit.Steps.step_1
  bool_step_2 = params.cryo_fit.Steps.step_2
  bool_step_3 = params.cryo_fit.Steps.step_3
  bool_step_4 = params.cryo_fit.Steps.step_4
  bool_step_5 = params.cryo_fit.Steps.step_5
  bool_step_6 = params.cryo_fit.Steps.step_6
  bool_step_7 = params.cryo_fit.Steps.step_7
  bool_step_8 = params.cryo_fit.Steps.step_8
  #bool_step_9 = params.cryo_fit.Steps.step_9
  
  returned = assign_model_name(params, starting_dir, inputs, params.cryo_fit.Input.model_file_name)
  model_file_with_pathways = returned[0]
  model_file_without_pathways = returned[1]
  
  # Options used for GUI based specification as well
  restraint_algorithm_minimization = params.cryo_fit.Options.restraint_algorithm_minimization
  emsteps = params.cryo_fit.Options.emsteps
  
  if (params.cryo_fit.Options.emweight_multiply_by == None): 
    params.cryo_fit.Options.emweight_multiply_by = 2
  emweight_multiply_by = params.cryo_fit.Options.emweight_multiply_by
  
  
  emwritefrequency = params.cryo_fit.Options.emwritefrequency
  no_rerun = params.cryo_fit.Options.no_rerun
  
  if (params.cryo_fit.Options.nstxtcout == None): # if nstxtcout is empty in GUI, there was an error. That's why we have this if clause.
    params.cryo_fit.Options.nstxtcout = 100
  nstxtcout = params.cryo_fit.Options.nstxtcout
  
  if (params.cryo_fit.Options.time_step_for_cryo_fit == None): 
    params.cryo_fit.Options.time_step_for_cryo_fit = 0.002
  time_step_for_cryo_fit = params.cryo_fit.Options.time_step_for_cryo_fit
  
  if (params.cryo_fit.Options.time_step_for_minimization == None): 
    params.cryo_fit.Options.time_step_for_minimization = 0.001
  time_step_for_minimization = params.cryo_fit.Options.time_step_for_minimization
  
  
  user_entered_number_of_steps_for_cryo_fit = params.cryo_fit.Options.number_of_steps_for_cryo_fit
  user_entered_number_of_steps_for_minimization = params.cryo_fit.Options.number_of_steps_for_minimization
  
  # Development options
  devel = params.cryo_fit.devel
  force_field = params.cryo_fit.force_field
  ignh = params.cryo_fit.ignh
  initial_cc_wo_min = params.cryo_fit.initial_cc_wo_min
  initial_cc_w_min = params.cryo_fit.initial_cc_w_min
  kill_mdrun_mpirun_in_linux = params.cryo_fit.kill_mdrun_mpirun_in_linux
  lincs_order = params.cryo_fit.lincs_order
  missing = params.cryo_fit.missing
  ns_type = params.cryo_fit.ns_type
  number_of_cores_to_use = params.cryo_fit.number_of_cores_to_use
  perturb_xyz_by = params.cryo_fit.perturb_xyz_by
  remove_metals = params.cryo_fit.remove_metals
  many_step_____n__dot_pdb = params.cryo_fit.many_step_____n__dot_pdb
  
  if (many_step_____n__dot_pdb == True):
    emweight_multiply_by = 1
    lincs_order = 1
    time_step_for_cryo_fit = 0.001
  
  number_of_steps_for_minimization = determine_number_of_steps_for_minimization(model_file_without_pathways,\
                                                                            model_file_with_pathways, \
                                                                            user_entered_number_of_steps_for_minimization, devel)
  params.cryo_fit.Options.number_of_steps_for_minimization = number_of_steps_for_minimization
  print "\tparams.cryo_fit.Options.number_of_steps_for_minimization (a real value that will be used eventually): ", \
    params.cryo_fit.Options.number_of_steps_for_minimization
  
  number_of_steps_for_cryo_fit = determine_number_of_steps_for_cryo_fit(model_file_without_pathways,\
                                                                            model_file_with_pathways, \
                                                                            user_entered_number_of_steps_for_cryo_fit, devel)
  if (nstxtcout > number_of_steps_for_cryo_fit):
    write_this = "\nnstxtcout (" + str(nstxtcout) + ") > number_of_steps_for_cryo_fit (" \
              + str(number_of_steps_for_cryo_fit) + ")"
    print write_this
    logfile.write(write_this)
    write_this = "Please reset so that nstxtcout < number_of_steps_for_cryo_fit to extract gro/pdb files properly later. Exit now.\n"
    print write_this
    logfile.write(write_this)
    exit(1)
    
  if (initial_cc_wo_min == True):
    no_rerun = True
    number_of_cores_to_use = str(2) # because of option choice above, it should be assigned as string
    number_of_steps_for_minimization = 0
    number_of_steps_for_cryo_fit = 100
  
  if (initial_cc_w_min == True):
    no_rerun = True
    number_of_cores_to_use = str(2) # because of option choice above, it should be assigned as string
    number_of_steps_for_cryo_fit = 100
    
  
  params.cryo_fit.Options.number_of_steps_for_cryo_fit = number_of_steps_for_cryo_fit
  print "\tparams.cryo_fit.Options.number_of_steps_for_cryo_fit (a real value that will be used eventually): ", \
    params.cryo_fit.Options.number_of_steps_for_cryo_fit
  
  steps_list = [bool_step_1, bool_step_2, bool_step_3, bool_step_4, bool_step_5, bool_step_6, bool_step_7\
                , bool_step_8]
  print "\tsteps_list: ", steps_list # this is shown in GUI
  make_new_steps_folder = True
  for i in range(len(steps_list)-1): # don't care step_8 for now
    if steps_list[i] == False:
      make_new_steps_folder = False
    
  if (make_new_steps_folder == True):
    remake_this_folder("steps")
  else:
    if (os.path.isdir("steps") != True):
      print "\tMake \"steps\" folder"
      os.mkdir("steps")
    else:
      print "\tkeep existing \"steps\" folder"
    
  command_path = locate_Phenix_executable()
  
  number_of_available_cores = know_total_number_of_cores()
  number_of_available_cores = number_of_available_cores[:-1] # to remove "\n" at the end
  this_is_test_for_each_step = False # by default
  
  if ((platform.system() == "Linux") and (kill_mdrun_mpirun_in_linux == True)):
    kill_mdrun_mpirun_in_linux()    
  if (steps_list[0] == True):
    this_is_test_for_each_step = step_1(logfile, command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, force_field, ignh, missing, remove_metals, cryo_fit_path)
    write_this = "Step 1 (Make gro and topology file by regular gromacs) is successfully ran\n"
    if (this_is_test_for_each_step == True):
      end_regression(starting_dir, write_this)
      #return 0 # exit the whole main fn as expected
    else:
      logfile.write(write_this)
    
  if (steps_list[1] == True):
    this_is_test_for_each_step = step_2(logfile, command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, force_field, \
           perturb_xyz_by, remove_metals)
    write_this = "Step 2 (Clean gro file to be compatible for amber03 forcefield) is successfully ran\n"
    if (this_is_test_for_each_step == True):
      end_regression(starting_dir, write_this)
    else:
      logfile.write(write_this)
  
  print "\trestraint_algorithm_minimization:",restraint_algorithm_minimization
  
  if str(restraint_algorithm_minimization) != "none_default": # this is default for "regular running" and "regression running"
    if (steps_list[2] == True):
      this_is_test_for_each_step = step_3(logfile, command_path, starting_dir, ns_type, restraint_algorithm_minimization, number_of_steps_for_minimization, \
             time_step_for_minimization, model_file_without_pathways, devel, cryo_fit_path)
      write_this = "Step 3 (Make a tpr file for minimization) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)
      
    if (steps_list[3] == True):
      this_is_test_for_each_step = step_4(logfile, command_path, starting_dir, ns_type, number_of_available_cores, \
                            number_of_cores_to_use, model_file_without_pathways, cryo_fit_path)
      write_this = "Step 4 (Minimize a gro file (to prevent \"blowup\" during MD Simulation)) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)
  else: #str(restraint_algorithm_minimization) = "none_default"
    if (steps_list[2] == True):
      this_is_test_for_each_step = step_3(logfile, command_path, starting_dir, ns_type, "none", number_of_steps_for_minimization, \
             time_step_for_minimization, model_file_without_pathways, devel, cryo_fit_path)
      write_this = "Step 3 (Make a tpr file for minimization) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)
      
    if (steps_list[3] == True):
      this_is_test_for_each_step = step_4(command_path, starting_dir, ns_type, number_of_available_cores, number_of_cores_to_use, \
             model_file_without_pathways, cryo_fit_path)
      write_this = "Step 4 (Minimize a gro file (to prevent \"blowup\" during MD Simulation)) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)
      
    cp_command_string = "cp steps/4_minimize/minimized_c_term_renamed_by_resnum_oc.gro . "
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    shutil.rmtree("steps/3_make_tpr_to_minimize")
    shutil.rmtree("steps/4_minimize")
    
    if (steps_list[2] == True):
      this_is_test_for_each_step = step_3(logfile, command_path, starting_dir, ns_type, "none_default", number_of_steps_for_minimization, \
             time_step_for_minimization, model_file_without_pathways, devel,cryo_fit_path)
      write_this = "Step 3 (Make a tpr file for minimization) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)
    
    if (steps_list[3] == True):
      this_is_test_for_each_step = step_4(command_path, starting_dir, ns_type, number_of_available_cores, number_of_cores_to_use, \
             model_file_without_pathways, cryo_fit_path)
      write_this = "Step 4 (Minimize a gro file (to prevent \"blowup\" during MD Simulation)) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)
  
  if (steps_list[4] == True):
    this_is_test_for_each_step = step_5(logfile, command_path, starting_dir, model_file_without_pathways, cryo_fit_path)
    write_this = "Step 5 (Make contact potential (restraints) and topology file with it) is successfully ran\n"
    if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
    else:
      logfile.write(write_this)
  
  if (steps_list[5] == True):
    this_is_test_for_each_step = step_6(logfile, command_path, starting_dir, model_file_without_pathways)
    write_this = "Step 6 (Make all charges of atoms be 0) is successfully ran\n"
    if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
    else:
      logfile.write(write_this)
  
    
  cc_has_been_increased = True # just an initial value
  charge_group_moved = True # just an initial value
  re_run_with_higher_map_weight = True # just an initial value
  
  returned = assign_map_name(params, starting_dir, inputs, params.cryo_fit.Input.map_file_name) # if needed, mrc_to_sit will run here
  map_file_with_pathways = returned[0]
  map_file_without_pathways = returned[1]
  
  if (model_file_without_pathways == "GTPase_activation_center_tutorial.pdb"):
    no_rerun = True
    number_of_steps_for_cryo_fit = 70000 # this is the ideal steps for reaching a plateau
  
  if (devel == True):
    no_rerun = True
    number_of_steps_for_cryo_fit = 100
    
  restart_w_longer_steps = False # this is a proper initial assignment
  re_run_with_higher_map_weight = False # this is a proper initial assignment
  
  # Iterate until any condition is met
  iteration_numner = 0
  many_stepxb = False
  while ((cc_has_been_increased == True) or (charge_group_moved == True) \
         or (re_run_with_higher_map_weight == True)):
    iteration_numner += 1
    
    if "regression_" in model_file_without_pathways:
      if (iteration_numner >= 5):
        break
    else:
      if (iteration_numner >= 20):
        break
    
    if (emweight_multiply_by > 100):
      break
    
    if ((this_is_test_for_each_step == True) \
       or (steps_list[0] == False and steps_list[1] == False and steps_list[2] == False \
           and steps_list[3] == False and steps_list[4] == False and steps_list[5] == False \
           and steps_list[6] == False and steps_list[7] == False)):
      break
    
    if (re_run_with_higher_map_weight == True):
      restart_w_longer_steps = False
      # w/o this, cryo_fit will think that restart_w_longer_steps is still True
      
    if (steps_list[6] == True):
      this_is_test_for_each_step = step_7(logfile, command_path, starting_dir, number_of_steps_for_cryo_fit, emweight_multiply_by, emsteps, \
             emwritefrequency, lincs_order, nstxtcout, time_step_for_cryo_fit, model_file_without_pathways, cryo_fit_path, many_step_____n__dot_pdb)
      write_this = "Step 7 (Make a tpr file for cryo_fit) is successfully ran\n"
      if (this_is_test_for_each_step == True):
        end_regression(starting_dir, write_this)
      else:
        logfile.write(write_this)

      
    if (steps_list[7] == True):
      results = step_8(logfile, command_path, starting_dir, number_of_available_cores, number_of_cores_to_use, 
             map_file_with_pathways, no_rerun, devel, restart_w_longer_steps, \
             re_run_with_higher_map_weight, model_file_without_pathways, cryo_fit_path)

      if (results == True): # this is a test for each step
        end_regression(starting_dir, "This is a test for each step, so break early of this step 7 & 8 loop")
      
      
      
      ################### (begin) check user_s_cc sanity
      cwd = os.getcwd()
      
      last_folder_name = os.path.basename(cwd)
      user_s_cc = ''
      if (last_folder_name != "8_cryo_fit"):
        user_s_cc = check_first_cc("steps/8_cryo_fit/cc_record")
      else:
        user_s_cc = check_first_cc("cc_record")
      
      if (user_s_cc == ''):
        print_this = "cryo_fit cannot calculate CC with a user input pdb file and map file. \ncc_record is not found. \n Please contact doonam.kim@gmail.com"
        print print_this
        logfile.write(print_this)
        return "failed" # flatly failed
      
      try:
        user_s_cc_rounded = str(round(float(user_s_cc), 3)) # if user_s_cc is stil '', "ValueError: could not convert string to float:"
      except:
        print_this = "cryo_fit cannot calculate CC with a user input pdb file and map file. Please contact doonam@lanl.gov"
        print print_this
        logfile.write(print_this)
        return "failed" # flatly failed

      if (float(user_s_cc) < 0.0001):
        write_this = "\nA user's provided input pdb file has less than 0.0001 cc\n"
        print write_this
        logfile.write(write_this)
        write_this = "Please read https://www.phenix-online.org/documentation/faqs/cryo_fit_FAQ.html#i-see-user-s-provided-atomic-model-had-0-0-cc-in-my-cryo-fit-overall-log\n"
        print write_this
        logfile.write(write_this)
        write_this = "Sleep 10,000 seconds, so that this error is recognized instantly \n"
        print write_this
        logfile.write(write_this)
        write_this = "Exit cryo_fit by ctrl+C \n"
        print write_this
        logfile.write(write_this)
        time.sleep(10000)
      ################ (end) check user_s_cc sanity

      if (results == "failed_with_nan_in_cc"):
        write_this = "\n\tStep 8 failed with nan error in cc calculation\n"
        print write_this
        logfile.write(write_this)
        write_this = "\nPlease read https://www.phenix-online.org/documentation/faqs/cryo_fit_FAQ.html\n"
        print write_this
        logfile.write(write_this)
        return "failed" # flat fail
      
      if (results == "failed"): # flat fail
        return "failed"
        # For a small peptide, this error was possible "The initial cell size (0.392390) is smaller than the cell size limit (0.421512), change options -dd, -rdd or -rcon, see the log file for details

      elif results == "re_run_w_smaller_MD_time_step":
        charge_group_moved = True
        re_run_with_higher_map_weight = False
        print "\tstep 7 & 8 will re-run with smaller time_step_for_cryo_fit (" + str(time_step_for_cryo_fit*0.5) + ")\n\n"
        logfile.write("\tstep 7 & 8 will re-run with smaller time_step_for_cryo_fit (" + str(time_step_for_cryo_fit*0.5) + ")\n\n")
        if (time_step_for_cryo_fit < 0.0001): # to avoid infinite loop
          write_this = "time_step_for_cryo_fit < 0.0001, exit now"
          print write_this
          logfile.write(write_this)
          break
        os.chdir( starting_dir )
        
      elif results == "re_run_with_longer_steps":
        if (no_rerun == True): # usually for development purpose
          write_this = "re_run_with_longer_steps is recommended, but no_rerun = True\n"
          print write_this
          logfile.write(write_this)
          
          write_this = "Step 8 (cryo_fit itself) is successfully ran\n"
          print write_this
          logfile.write(write_this)
          
          this_is_test_for_each_step = step_final(logfile, command_path, starting_dir, model_file_without_pathways, no_rerun) # just to arrange final output
          
          write_this = "result was re_run_with_longer_steps, but this is a test for each step"
          if (this_is_test_for_each_step == "failed"):
            exit(1)
          elif (this_is_test_for_each_step == True):
            end_regression(starting_dir, write_this)
          else:
            logfile.write(write_this)
        
        restart_w_longer_steps = True
        re_run_with_higher_map_weight = False
        
        
        # Copy for a next restart step
        if (os.path.isfile("state.cpt") == False):
          write_this = 'state.cpt not found, step_8 may be full of stepxb_nx.pdb. \nVisit https://www.phenix-online.org/documentation/faqs/cryo_fit_FAQ.html \n'
          print write_this
          logfile.write(write_this)
          
          if (emweight_multiply_by == 1):
            many_stepxb = True
            break
          
          # this long 1 line is essential for proper writing into log file
          write_this = 'Maybe emweight_multiply_by (' + str(emweight_multiply_by) + ') is too high.\n Cryo_fit will divide emweight_multiply_by by 3 (then round, so that emweight_multiply_by becomes ' + str(int(round(emweight_multiply_by/3,0))) + ') and re-run again.\n'
          
          emweight_multiply_by = int(round((emweight_multiply_by/3), 0))
          if (emweight_multiply_by < 1):
            emweight_multiply_by = 1
          
          print write_this
          logfile.write(write_this)
          
          os.chdir( starting_dir ) # needed for re-running
          continue
        
        if (many_stepxb == True):
          break
        
        shutil.copy("state.cpt", "../..")
  
        charge_group_moved = False # just initial value
        
        #number_of_steps_for_cryo_fit = number_of_steps_for_cryo_fit * 2 # Karissa seems to be concerned over speed
        #number_of_steps_for_cryo_fit = number_of_steps_for_cryo_fit + 5000 # for a unknown reason, this method resulted in "0 step only run" eventually
        number_of_steps_for_cryo_fit = int(number_of_steps_for_cryo_fit * 1.3) 
        
        write_this = "\nStep 8 (cryo_fit itself) is ran well, but correlation coefficient values tend to be increased recently.\n"
        print write_this
        logfile.write(write_this)
        
        write_this = "Therefore, step 7 & 8 will re-run with longer steps (including all previously ran steps, up to " + str(number_of_steps_for_cryo_fit) + ")\n\n"
        print write_this
        logfile.write(write_this)
        
        # aim_this_step_when_restart.txt is essential to extract gro from traj.xtc when restarted WITH LONGER STEPS (not restarted w/ higher map weight)
        restart_record = open("../../aim_this_step_when_restart.txt", "a+")
        write_this = str(int(number_of_steps_for_cryo_fit)) + "\n" 
        restart_record.write(write_this)
        restart_record.close()
        
        if (number_of_steps_for_cryo_fit > 1000000000000000 ): # to avoid infinite loop
          write_this = "number_of_steps_for_cryo_fit > 1000000000000000, exit now"
          print write_this
          logfile.write(write_this)
          break # no need to use exit(1) since this break will break this while loop
        os.chdir( starting_dir ) # needed for re-running

      elif results == "re_run_with_higher_map_weight":
        write_this = "\nStep 8 itself ran well, but correlation coefficient values tend to decrease over the last 30 steps\n"
        print write_this
        logfile.write(write_this)
        
        if (max_emweight_multiply_by <= emweight_multiply_by*2):
          emweight_multiply_by = emweight_multiply_by * 2 # this new emweight_multiply_by will be used at step 7 (tpr file generation)
        else:
          write_this = "\nemweight_multiply_by x 2 exceeds max_emweight_multiply_by. Therefore, cryo_fit will finish now.\n"
          print write_this
          logfile.write(write_this)
          break
        
        # this check is important to avoid infinite loop and "vtot is inf: inf" error
        if (emweight_multiply_by > 1000 ):
          # 1,024 resulted in "vtot is inf: inf" with christl's molecule
          write_this = "emweight_multiply_by > 1000, cryo_fit will exit, because 1,024 emweight_multiply_by once resulted in \"vtot is inf: inf\". This often means that user's initial atomistic model already fit to cryo-EM map well."
          print write_this
          logfile.write(write_this)
            # "; emweight is the energetic weight of the em term.  Since the total stabilizing energy in an all-atom SBM is set to the number of atoms,
            #we typically give a weight that is 1-2 times the overall stabilizing energy.  Here, AKE has 22k atoms, so we will use a weight of 44k."
            
          break # no need to use exit(1) since this break will break this while loop
          
        write_this = "Therefore, step 7 & 8 will re-run with a higher emweight_multiply_by (e.g. " + str(emweight_multiply_by) + ")\n\n"
        print write_this
        logfile.write(write_this)
        re_run_with_higher_map_weight = True
        
        os.chdir( starting_dir ) # needed for re-running

      else: # normal ending of cryo_fit
        print_this = "\nA user's provided input pdb file has " + str(round(float(user_s_cc), 3)) + " cc\n"
        print print_this
        logfile.write(print_this)
        
        charge_group_moved = False
        cc_has_been_increased = False
        re_run_with_higher_map_weight = False
  if (many_stepxb == True):
    logfile.write("\nStep 8 ran, but not successfully.\n")
  else:
    logfile.write("\nStep 8 ran.\n")

  if ((initial_cc_wo_min == False) and (initial_cc_w_min == False) and (many_stepxb == False)):
    this_is_test_for_each_step = step_final(logfile, command_path, starting_dir, model_file_without_pathways, \
                              cryo_fit_path, no_rerun) # just to arrange final output
    if (os.path.isfile("steps/1_make_gro/prefix_of_chain_ID_removed") == True):
      write_this = "The 4th character of residue name (prefix_of_chain ID) is removed. Please see https://www.phenix-online.org/documentation/faqs/cryo_fit_FAQ.html#how-can-i-use-double-digit-character-id-pdb-file or email doonam.kim@gmail.com\n\n"
      print write_this
      logfile.write(write_this)
      
    write_this = "step final is done"
    if (this_is_test_for_each_step == "failed"):
      exit(1)
    elif (this_is_test_for_each_step == True):
      end_regression(starting_dir, write_this)
    else:
      if "regression_" in model_file_without_pathways: # regression for all steps
        return True
      else: # regular running
        return results
            
  # keep for now for this cc draw
  #if (steps_list[8] == True):
  #  step_9(command_path, starting_dir, model_file_without_pathways, map_file_without_pathways)
############################# end of run_cryo_fit function


# parse through command line arguments
def cmd_run(args, validated=False, out=sys.stdout):
  time_total_start = time.time()
  print_author()
  if (len(args) < 2 and validated==False):
    print >> out, "-"*79
    print >> out, "                               cryo_fit"
    print >> out, "-"*79
    print >> out, legend
    print >> out, "-"*79
    #master_params.show(out=out)
    explanation_only = True
    return explanation_only
    
  log = multi_out()
  log.register("stdout", out)
  
  log_file_name = "cryo_fit.overall_log"
  logfile = open(log_file_name, "w")
  logfile.write("Overall log of cryo_fit\n\n")
  log.register("logfile", logfile)
  
  print >> log, "Input parameters:", args

  input_command_file = open("cryo_fit.input_command", "w")
  logfile.write("Input command: phenix.cryo_fit ")
  input_command_file.write("phenix.cryo_fit ")
  for i in range(len(args)):
    input_command_file.write(args[i] + " ")
    logfile.write(args[i] + " ")
  input_command_file.write("\n")
  logfile.write("\n\n")
  input_command_file.close()

  # very simple parsing of model and map
  for i, arg in enumerate(args):
    #if arg.endswith('.cif') or arg.endswith('.ent') or arg.endswith('.pdb'): # EMD-3981 has 6exv.ent instead of .pdb
    if arg.endswith('.cif') or arg.endswith('.pdb'): # .ent brought an error in GUI
      if arg.find('=')==-1:
        args[i]='model=%s' % arg
    elif arg.endswith('.ccp4') or arg.endswith('.map') or arg.endswith('.mrc') or arg.endswith('.sit'):
      if arg.find('=')==-1:
        args[i]='map=%s' % arg
  
  # for mrc_to_sit
  crystal_symmetry=None
  inputs = mmtbx.utils.process_command_line_args(args = args,
        cmd_cs=crystal_symmetry,
        master_params = master_phil)
  
  time_process_command_line_args_start = time.time()
  print "\tReading user provided map started..."
  print "\t(If a user provided a big .sit file like 1.6 GB, this may take more than 6 minutes)"
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="cryo_fit",
  )
  
  pdbs = []
  sits = []
  phils = []
  phil_args = []
  for arg in args:
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
      elif arg.endswith('.map'): # not the smartest
        sits.append(arg)
      else:
        try :
          file_phil = phil.parse(file_name=arg)
        except RuntimeError :
          pass
        else :
          phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  working_phil.show()
  working_params = working_phil.extract()
  
  if (not validated):
    validate_params(working_params)
  time_process_command_line_args_end = time.time()
  
  starting_dir = os.getcwd()
  print "\tCurrent working directory: %s" % starting_dir
  
  results = run_cryo_fit(logfile, working_params, inputs)
  
  time_total_end = time.time()
  time_took = show_time(time_total_start, time_total_end)
  print "\nTotal cryo_fit", time_took
  
  write_this = "\nTotal cryo_fit " + time_took + "\n"
  logfile = open(log_file_name, "a+") # append
  logfile.write(write_this)
  logfile.close()
  
  if (results == "failed") or (results == "re_run_w_smaller_MD_time_step"): # errored
    exit(1)
  
  if (results == True): # regression test for all steps
    end_regression(starting_dir,"end regression for all steps")
  
  return results
################### end of cmd_run function


# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    import os
    from wxGUI2 import utils
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = cmd_run(args=self.args, validated=True, out=sys.stdout)
    return result
# =============================================================================


if (__name__ == "__main__") :
  cmd_run(args = sys.argv[1:])
