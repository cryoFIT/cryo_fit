# LIBTBX_SET_DISPATCHER_NAME phenix.cryo_fit
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT 

# Steps of cryo_fit:
# 1_Make_gro
# 2_Clean_gro
# 3_Prepare_to_Minimize
# 4_Minimize
# 5_Make_constraints
# 6_Make_0_charge
# 7_Make_tpr_for_EM_map_fitting
# 8_EM_map_fitting_itself
# 9_Arrange_outputs
# (not used now) # 9_Draw_a_figure_of_cc

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from iotbx import file_reader
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out
import mmtbx.utils
import platform
import shutil # for rmdir
from subprocess import check_output

# this is needed to import all common functions
path = check_output(["which", "phenix.cryo_fit"])
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
         density map structurally (for example by USCF chimera)
    - A .ccp4 (MRC) or .map (MRC) or .sit (Situs) file, a cryo EM density map 
    
Output folder:
    - cryo_fitted.x: Fitted biomolecule structure to a target cryo-EM map
    - Correlation coefficient record: Record of correlation coefficient between cryo-EM map
    and current biomolecule structure
      
Usage example with minimum input requirements (all other options will run with default values):
    - phenix.cryo_fit tRNA.pdb tRNA.map

Usage example with step 7~8 only
    - phenix.cryo_fit tRNA.pdb tRNA.map step_1=False step_2=False step_3=False step_4=False step_5=False step_6=False
    
Most useful options (GUI has more explanation about these):
    - emweight_multiply_by
    - number_of_cores_to_use
    - number_of_steps_for_cryo_fit
    - number_of_steps_for_minimization
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
    .help = (.cif/.pdb) Either a homology model or a model from different organism/experimental method.
    .style = file_type:pdb bold input_file
  map_file_name = None
    .type = path
    .short_caption = Target map file 
    .help = Cryo-EM map file (.ccp4/.map/.sit)
    .style = bold input_file
  cryo_fit_path = None
    .type = path
    .short_caption = gromacs_cryo_fit installation path
    .help = Path to gromacs_cryo_fit installation.
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
  step_9 = False
    .type = bool
    .short_caption = 8. Show cc (optional)
    .help = Show correlation coefficient change
}
Options
{
  constraint_algorithm_minimization = *default none none_default
    .type = choice
    .help = Try none or none_default if a user see this error during minimization "Too many lincs warnings" \
            none_default will minimize twice. \
            First with constraint_algorithm_minimization = none \
            then with constraint_algorithm_minimization = default
  emsteps = None
    .type = int
    .short_caption = EM steps
    .help = emsteps is the number of integration steps between re-evaluation of the simulated map and forces. \
            The longer the emsteps be, the faster overall cryo_fit running time. \
            If it is left blank, the cryo_fit will automatically determine the emsteps
  emweight_multiply_by = 8
    .type = int
    .short_caption = EM weight multiply by this number
    .help = Multiply by this number to the number of atoms for weight for cryo-EM map bias. \
            For example, emweight = (number of atoms in gro file) x (emweight_multiply_by which is 8) \
            The higher the weight, the stronger bias toward EM map rather than MD force field and stereochemistry preserving constraints. \
            If user's map has a better resolution, higher value of emweight_multiply_by is recommended since map has much information. \
            If user's map has have a worse resolution, lower value of emweight_multiply_by is recommended for more likely geometry. \
            If CC (correlation coefficient) needs to be improved faster, higher number of emweight_multiply_by is recommended.
  emwritefrequency = None
    .type = int
    .short_caption = EM write frequency
    .help = Frequency with which the simulated maps are written to file. \
            If this frequency is too small, it can cause extremely large amounts of data to be written.\
            If it is left blank, the cryo_fit will use default value of 1,000,000
  number_of_steps_for_minimization = None
    .type = int
    .short_caption = Number of steps for minimization
    .help = Specify number of steps for minimization. \
           If it is left blank, cryo_fit will estimate it automatically depending on molecule size.
  number_of_steps_for_cryo_fit = None
    .type = int
    .short_caption = Number of steps for cryo_fit
    .help = Specify number of steps for cryo_fit. \
           If it is left blank, cryo_fit will estimate it automatically depending on molecule size and cc trend.
  time_step_for_cryo_fit = 0.002
    .type = float
    .short_caption = Time step for MD simulation during cryo_fit
    .help = Default value is 0.002. Try 0.001 if a user see this error during cryo_fit \
    "Fatal error: A charge group moved too far between two domain decomposition steps \
    This usually means that your system is not well equilibrated"
  time_step_for_minimization = 0.001
    .type = float
    .help = Default value is 0.001. Try 0.0005 if a user see this error during minimization
    "Fatal error: A charge group moved too far between two domain decomposition steps \
    This usually means that your system is not well equilibrated"
}

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
kill_mdrun_mpirun_in_linux = False
  .type = bool
  .short_caption = If true, kill any existing md run and mpirun.
  .help = This command works only for Linux
lincs_order = None
  .type = int
  .short_caption = LINear Constraint Solver
  .help = The accuracy in set with lincs-order, which sets the number of matrices \
          in the expansion for the matrix inversion. \
          If it is not specified, the cryo_fit will use 4.
no_rerun = False
  .type = bool
  .short_caption = If true, no_rerun
missing = True
  .type = bool
  .short_caption = If true, Continue when atoms are missing, dangerous
  .help = http://manual.gromacs.org/programs/gmx-pdb2gmx.html
ns_type = *grid simple
  .type = choice
  .short_caption = Method to determine neighbor list (simple, grid) during minimization
  .help = "Grid" is needed for domain decomposition (dd) for faster execution and ran well with tRNA, \
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
remove_metals = True
  .type = bool
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
    - Doo Nam Kim, Serdal Kirmizialtin, Nigel Moriarty, Tom Terwilliger, Billy Poon
 %s""" % ("-"*78, version, "-"*78)
# end of print_author()

def show_header(title):
  print "\n"
  print '#'*105
  number_of_remaining_sharp = 105 - len(title)
  put_this_number_of_sharp = int(int(number_of_remaining_sharp)/2)
  print '#'*(put_this_number_of_sharp-1) + " " + title + " " + '#'*(put_this_number_of_sharp-1)
  print '#'*105
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
# end of validate_params function

def step_1(logfile, command_path, starting_dir, model_file_with_pathways, starting_pdb_without_path, \
           force_field, ignh, missing, remove_metals, cryo_fit_path, *args):
  show_header("Step 1: Make gro and topology file by regular gromacs")
  remake_and_move_to_this_folder(starting_dir, "steps/1_make_gro")

  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  cp_command_string = ''
  for i in range(len(splited_starting_dir)):
    if splited_starting_dir[i] == "phenix_regression":
      this_is_test = True
      print "\tthis_is_test = True"
    
  cw_dir = os.getcwd()
  print "\tCurrent working directory: %s" % cw_dir
  
  cp_command_string = "cp " + model_file_with_pathways + " ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  start = time.time()
  cp_command_string = "cp " + command_path + "steps/1_make_gro/1_before_pdb2gmx_prepare_pdb.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  if (starting_pdb_without_path.find("_cleaned_for_gromacs") == -1):
    run_this = "python 1_before_pdb2gmx_prepare_pdb.py " + starting_pdb_without_path + " 0 0 0 " + \
               str(remove_metals)
    print "\tcommand: ", run_this
    f_out = open('log.step_1_1_before_pdb2gmx_prepare_pdb', 'wt')
    write_this_input_command = run_this + "\n"
    f_out.write(write_this_input_command)
    f_out.close()
    libtbx.easy_run.fully_buffered(run_this)

  cp_command_string = "cp " + command_path + "steps/1_make_gro/2_runme_make_gro.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  os.chdir (starting_dir)
  
  cw_dir = os.getcwd()
  print "\tCurrent working directory: %s" % cw_dir
  
  number_of_atoms_in_input_pdb = know_number_of_atoms_in_input_pdb(model_file_with_pathways)  
  if (number_of_atoms_in_input_pdb < 7000): # tRNA for development (transmin1_gro.pdb). This is just for short test purpose only
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
    this_step_was_successfully_ran = check_whether_the_step_was_successfully_ran("Step 1", check_this_file)
  if (this_step_was_successfully_ran == "failed"):
    logfile.write("Step 1 didn't run successfully")
    logfile.close()
    print "Step 1 didn't run successfully"
    print "\nUser's command was"
    f_in = open('../../cryo_fit.input_command')
    for line in f_in:
      print line

    bool_enable_mpi = know_output_bool_enable_mpi_by_ls()
    
    print "\nphenix.cryo_fit alone without any arguments introduces full options."
    print "Please email phenixbb@phenix-online.org or doonam@lanl.gov for any feature request/help."
    exit(1)
  print "Step 1", (show_time(start, end))
  return this_is_test
# end of step_1 function

def step_2(command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, \
           force_field, perturb_xyz_by, remove_metals):
  show_header("Step 2: Clean gro file to be compatible for amber03 forcefield")
  os.chdir (starting_dir)
  remake_and_move_to_this_folder(starting_dir, "steps/2_clean_gro")

  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  cp_command_string = ''
  
  if (model_file_without_pathways != "regression_pdb5khe.pdb"):  
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        print "this_is_test = True"
        cp_command_string = "cp ../../data/input_for_step_2/*_cleaned_for_gromacs_by_pdb2gmx.gro ."

  if (this_is_test == False):
    cp_command_string = "cp ../1_make_gro/*.gro ."
  
  libtbx.easy_run.fully_buffered(cp_command_string) #copy step_1 output

  start_time_renaming = time.time()
  print "\nStep 2: Add C prefix to terminal amino acid/nucleic acid for minimization by gromacs"
  command_script = "cp " + command_path + "steps/2_clean_gro/*.py ."
  libtbx.easy_run.fully_buffered(command_script)

  command_script = "python 1_rename_term_res_to_Cres.py " # there will be only 1 gro file, so it is ok
  print "\tcommand: ", command_script
  libtbx.easy_run.fully_buffered(command_script)

  for this_file in glob.glob("*_c_term_renamed_by_resnum_oc.gro"): # there will be only one file like this
    command_string = "python 4_slightlty_change_xyz_for_no_more_000.py " + this_file + " " + str(perturb_xyz_by)
    print "\tcommand: ", command_string 
    libtbx.easy_run.call(command_string)
  
  this_step_was_successfully_ran = "failed" # just an initial value
  for check_this_file in glob.glob("*.gro"): # there will be only "will_be_minimized_cleaned.gro"
    this_step_was_successfully_ran = check_whether_the_step_was_successfully_ran("Step 2", check_this_file)

  if (this_step_was_successfully_ran == "failed"):
    color_print (("Step 2 didn't run successfully"), 'red')
    exit(1)
  
  for check_this_file in glob.glob("*.gro"): # there will be only one file like this # to work in Karissa's old MacOS
    command_string = "mv " + check_this_file + " will_be_minimized_cleaned.gro"
    libtbx.easy_run.call(command_string)
    
  end_time_renaming = time.time()
  print "Step 2", (show_time(start_time_renaming, end_time_renaming))
  return this_is_test
# end of step_2 (clean gro) function


def step_3(logfile, command_path, starting_dir, ns_type, constraint_algorithm_minimization, number_of_steps_for_minimization, \
           time_step_for_minimization, model_file_without_pathways, devel, cryo_fit_path):
  show_header("Step 3: Make a tpr file for minimization")
  os.chdir (starting_dir)

  remake_this_folder("steps/3_make_tpr_to_minimize")
  remake_and_move_to_this_folder(starting_dir, "steps/3_make_tpr_to_minimize")

  cp_command_script = "cp " + command_path + "steps/3_make_tpr_to_minimize/minimization_template.mdp ."
  libtbx.easy_run.fully_buffered(cp_command_script)
  
  if (devel == True):
    number_of_steps_for_minimization = 10
  
  print "\tBe number_of_steps_for_minimization as ", number_of_steps_for_minimization
  with open("minimization_template.mdp", "rt") as fin:
    with open("minimization.mdp", "wt") as fout:
      for line in fin:
        splited = line.split()
        if splited[0] == "nsteps":
          new_line = "nsteps  = " + str(number_of_steps_for_minimization) + " ; Maximum number of minimization steps to perform\n"
          fout.write(new_line)
        elif splited[0] == "ns_type":
          new_line = "ns_type  = " + str(ns_type) + " ; Method to determine neighbor list (simple, grid)\n"
          fout.write(new_line)
        else:
          fout.write(line)
      print "\ttime_step_for_minimization:", time_step_for_minimization
      if time_step_for_minimization != 0.001:
        print "\ttime_step_for_minimization != 0.001"
        new_line = "\ndt = " + str(time_step_for_minimization) + "\n"
        fout.write(new_line)
      print "\tconstraint_algorithm_minimization:", constraint_algorithm_minimization
      if str(constraint_algorithm_minimization) == "None" or str(constraint_algorithm_minimization) == "none":
        print "\tconstraint_algorithm_minimization = none"
        new_line = "\nconstraint-algorithm: none\n"
        fout.write(new_line)
    fout.close()
  fin.close()
  
  cp_command_script = "cp " + command_path + "steps/3_make_tpr_to_minimize/runme_make_tpr.py ."
  libtbx.easy_run.fully_buffered(cp_command_script)
  
  this_is_test_for_tRNA = False
  splited_starting_dir = starting_dir.split("/")
  cp1_command_string = ''
  cp2_command_string = ''
    
  for i in range(len(splited_starting_dir)):
    if (model_file_without_pathways != "regression_pdb5khe.pdb"):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test_for_tRNA = True
        cp1_command_string = "cp ../../data/input_for_step_3/* ."
        
  if (this_is_test_for_tRNA == False):
    if str(constraint_algorithm_minimization) != "none_default":
      cp1_command_string = "cp ../2_clean_gro/*.gro . "
    else:
      cp1_command_string = "mv ../../minimized_c_term_renamed_by_resnum_oc.gro . "
    cp2_command_string = "cp ../1_make_gro/*.top . "
    libtbx.easy_run.fully_buffered(cp2_command_string)
  
  libtbx.easy_run.fully_buffered(cp1_command_string) #copy step_2 output

  command_string = "python runme_make_tpr.py " + str(cryo_fit_path)
  print "\tcommand: ", command_string
  start = time.time()
  libtbx.easy_run.call(command_string)
  end = time.time()

  check_whether_the_step_3_was_successfully_ran(logfile, "to_minimize.tpr")
  print "Step 3", (show_time(start, end))
  os.chdir( starting_dir )
  return this_is_test_for_tRNA
# end of step_3 (prepare minimization) function


def step_4(command_path, starting_dir, ns_type, number_of_available_cores, \
           number_of_cores_to_use, model_file_without_pathways, cryo_fit_path):
  show_header("Step 4: Minimize a gro file (to prevent \"blowup\" during Molecular Dynamics Simulation)")
  os.chdir (starting_dir)
  
  print "\nStep 4-1: Minimization itself"
  remake_and_move_to_this_folder(starting_dir, "steps/4_minimize")

  cp_command_script = "cp " + command_path + "steps/4_minimize/runme_minimize.py ."
  libtbx.easy_run.fully_buffered(cp_command_script)

  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  cp_command_string = ''
  if (model_file_without_pathways != "regression_pdb5khe.pdb"):
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        cp_command_string = "cp ../../data/input_for_step_4/* ."
  if (this_is_test == False):
    cp_command_string = "cp ../3_make_tpr_to_minimize/to_minimize.tpr ."
  
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  # when there are both mpi and thread cryo_fit exist, thread cryo_fit was used in commandline mode
  command_string = "python runme_minimize.py to_minimize.tpr " + str(command_path) + " " + \
                str(ns_type) + " " + str(number_of_available_cores) + " " + str(2) + " " + str(cryo_fit_path)
              # set number_of_cores_to_use = 2 to minimize a possibility of having cell size error
  print "\tcommand: ", command_string
  print "\n\tA user can check progress at ", starting_dir + "/steps/4_minimize\n"
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
    
  returned = check_whether_the_step_was_successfully_ran("Step 4-1", final_gro_file_name)
  if returned == "failed":
    bool_enable_mpi = know_output_bool_enable_mpi_by_ls()
    if bool_enable_mpi == True:
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
  cp_command_string = "cp " + command_path + "steps/2_clean_gro/*_rename_term_res_to_Cres*.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  command_string = "python 1_rename_term_res_to_Cres.py "
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)

  check_whether_the_step_was_successfully_ran("Step 4-2", "minimized_c_term_renamed_by_resnum_oc.gro")
  os.chdir( starting_dir )
  return this_is_test
# end of step_4 (minimization) function
    
    
def step_5(command_path, starting_dir, model_file_without_pathways, cryo_fit_path):
  show_header("Step 5: Make contact potential (constraints) and topology file with it")
  remake_and_move_to_this_folder(starting_dir, "steps/5_make_constraints")
  
  start = time.time()
  cp_command_string = "cp " + command_path + "steps/5_make_constraints/runme_make_contact_potential.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  cp_command_string = ''
  
  if (model_file_without_pathways != "regression_pdb5khe.pdb"):
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        cp_command_string = "cp ../../data/input_for_step_5/* ."
        
  if (this_is_test == False):
    cp_command_string = "cp ../4_minimize/*.gro ."
    # there will be minimized_c_term_renamed_by_resnum_oc.gro
  
  libtbx.easy_run.fully_buffered(cp_command_string) #copy step_4 output
  
  command_string = "python runme_make_contact_potential.py *.gro " + str(command_path) + " " + str(cryo_fit_path)
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)

  check_whether_the_step_was_successfully_ran("Step 5", "disre2.itp")
  end = time.time()
  print "Step 5", (show_time(start, end))
  
  #color_print ((show_time("Step 5", start, end)), 'green')
  # [keep] looks as "[32mStep 5 finished in 10.66 seconds (wallclock).[0m" in GUI
  
  os.chdir( starting_dir )
  return this_is_test
# end of step_5 (make constraints) function


def step_6(command_path, starting_dir, model_file_without_pathways):
  show_header("Step 6: Make all charges of atoms be 0")
  remake_and_move_to_this_folder(starting_dir, "steps/6_make_0_charge")

  cp_command_string = "cp " + command_path + "steps/6_make_0_charge/changetop.awk ."
  libtbx.easy_run.fully_buffered(cp_command_string)
      
  cp_command_string = "cp " + command_path + "steps/6_make_0_charge/runme_make_0_charge.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  cp_command_string = ''

  if (model_file_without_pathways != "regression_pdb5khe.pdb"):
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        cp_command_string = "cp ../../data/input_for_step_6/* ."
        
  if (this_is_test == False):
    cp_command_string = "cp ../5_make_constraints/*including_disre2_itp.top ."
  # In normal case, there will be minimized_c_term_renamed_by_resnum_oc_including_disre2_itp.top
  
  libtbx.easy_run.fully_buffered(cp_command_string) #copy step_5 output
  
  command_string = "python runme_make_0_charge.py *.top"
  print "\tcommand: ", command_string
  start = time.time()
  libtbx.easy_run.fully_buffered(command_string)
  end = time.time()

  for check_this_file in glob.glob("*_0_charge.top"): # there will be only one file like this
    check_whether_the_step_was_successfully_ran("Step 6", check_this_file)
    
  print "Step 6", (show_time(start, end))
  os.chdir( starting_dir )
  return this_is_test
# end of step_6 (neutralize) function
    

def step_7(command_path, starting_dir, number_of_steps_for_cryo_fit, emweight_multiply_by, \
           emsteps, emwritefrequency, lincs_order, time_step_for_cryo_fit, \
           model_file_without_pathways, cryo_fit_path):
  show_header("Step 7 : Make a tpr file for cryo_fit")
  remake_and_move_to_this_folder(starting_dir, "steps/7_make_tpr_with_disre2")

  cp_command_string = "cp " + command_path + "steps/7_make_tpr_with_disre2/template_for_cryo_fit.mdp ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  cp1_command_string = ''
  cp2_command_string = ''
  cp3_command_string = ''
  
  if (model_file_without_pathways != "regression_pdb5khe.pdb"):
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        cp1_command_string = "cp ../../data/input_for_step_7/* ."
        
  if (this_is_test == False):
    cp1_command_string = "cp ../5_make_constraints/*.gro ." # there will be minimized_c_term_renamed_by_resnum_oc.gro only
    cp2_command_string = "cp ../5_make_constraints/disre2.itp ."
    libtbx.easy_run.fully_buffered(cp2_command_string)
    cp3_command_string = "cp ../6_make_0_charge/*0_charge.top ." # there is only one *0_charge.top file
    libtbx.easy_run.fully_buffered(cp3_command_string)
    
  libtbx.easy_run.fully_buffered(cp1_command_string) #copy step_5 output
  
  print "\tBe number_of_steps_for_cryo_fit as ", number_of_steps_for_cryo_fit
  with open("template_for_cryo_fit.mdp", "rt") as fin:
    with open("for_cryo_fit.mdp", "wt") as fout:
      for line in fin:
        splited = line.split()
        if splited[0] == "dt":
          new_line = "dt = " + str(time_step_for_cryo_fit) + "\n"
          fout.write(new_line)
        elif splited[0] == "emsteps":
          if (emsteps == None):
              new_line = "emsteps = " + str(int(number_of_steps_for_cryo_fit/80)) + "\n" # to make cryo_fit step 8 faster
              # when emsteps is too sparse, cc went to become worse
              fout.write(new_line)
          else:
            new_line = "emsteps = " + str(emsteps) + "\n"
            fout.write(new_line)
        elif splited[0] == "emweight":
          number_of_atoms_in_gro = return_number_of_atoms_in_gro()
          print "\temweight_multiply_by:", emweight_multiply_by
          new_line = "emweight = " + str(int(number_of_atoms_in_gro)*int(emweight_multiply_by)) + "\n"
          fout.write(new_line)
        elif splited[0] == "emwritefrequency":
          if (emwritefrequency == None): # default is 1,000,000, because I don't see any usefulness of writing intermediate .sit file
            fout.write(line)
          else:
            new_line = "emwritefrequency = " + str(emwritefrequency) + "\n"
            fout.write(new_line)
        elif splited[0] == "lincs-order":
          if (lincs_order == None):
            fout.write(line)
          else:
            new_line = "lincs-order  = " + str(lincs_order) + "\n"
            fout.write(new_line)
        elif splited[0] == "nsteps":
          new_line = "nsteps  = " + str(number_of_steps_for_cryo_fit) + " ; Maximum number of steps to perform cryo_fit\n"
          fout.write(new_line)
        else:
          fout.write(line)
    fout.close()
  fin.close()
  
  cp_command_string = "cp " + command_path + "steps/7_make_tpr_with_disre2/runme_make_tpr_with_disre2.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)

  start_make_tpr = time.time()
  command_string = "python runme_make_tpr_with_disre2.py " + str(command_path) + " " + str(cryo_fit_path)
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string) 
  end_make_tpr = time.time()
  
  check_whether_the_step_was_successfully_ran("Step 7", "for_cryo_fit.tpr")
  print "Step 7", (show_time(start_make_tpr, end_make_tpr))
  os.chdir( starting_dir )
  return this_is_test
# end of step_7 (make tpr for cryo_fit) function
           
           
def step_8(logfile, command_path, starting_dir, number_of_available_cores, number_of_cores_to_use, \
       map_file_with_pathways, no_rerun, devel, restart_w_longer_steps, re_run_with_higher_map_weight, \
       model_file_without_pathways, cryo_fit_path):
  show_header("Step 8: Run cryo_fit")
  print "\tmap_file_with_pathways:",map_file_with_pathways
  
  remake_and_move_to_this_folder(starting_dir, "steps/8_cryo_fit")
  
  command_string = "cp " + command_path + "steps/8_cryo_fit/* ."
  libtbx.easy_run.fully_buffered(command_string)
  
  if (str(restart_w_longer_steps) == "True"):
    command_string = "cp ../../state.cpt . "
    libtbx.easy_run.fully_buffered(command_string)
  
  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  if (model_file_without_pathways != "regression_pdb5khe.pdb"):
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        print "\tthis_is_test:", this_is_test
  
  command_string = "python runme_cryo_fit.py " + str(command_path) + " " + str(number_of_available_cores) \
              + " " + number_of_cores_to_use + " " + map_file_with_pathways + " " + str(starting_dir) \
              + " " + str(this_is_test) + " " + str(restart_w_longer_steps) + " " + str(cryo_fit_path)
  print "\n\tcommand: ", command_string
  print "\n\tA user can check progress at ", starting_dir + "/steps/8_cryo_fit\n"
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
  
  returned = check_whether_the_step_was_successfully_ran("Step 8", "cc_record")
  
  if (returned != "success"):
    print "Step 8 (Run cryo_fit) didn't run successfully"
    logfile.write("Step 8 (Run cryo_fit) didn't run successfully\n")
    searched = search_charge_in_md_log()
    if searched == 0: # no "charge group... " message in md.log
      return "failed"
    else:
      return "re_run_w_smaller_MD_time_step"
    
  if ((str(restart_w_longer_steps) == "False") and (str(re_run_with_higher_map_weight) == "False")):
    user_s_cc = get_users_cc("cc_record")
    user_s_cc = round(float(user_s_cc),3)
    write_this = "\nUser's provided atomic model had " + str(user_s_cc) + " cc\n\n"
    logfile.write(write_this)
    
  cc_has_been_increased = check_whether_cc_has_been_increased(logfile, "cc_record")
  print "\t\tVerdict of cc_has_been_increased function in the last 30 cc evaluations:", cc_has_been_increased
  
  if (devel == True):
      no_rerun = True
  if (no_rerun == False):
    if cc_has_been_increased == True:
      return "re_run_with_longer_steps"
    elif cc_has_been_increased == "re_run_with_higher_map_weight":
      return "re_run_with_higher_map_weight"
    else:
      print "\tcc has been saturated, so go ahead to the next step"
  
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
  
  this_is_test = False
  splited_starting_dir = starting_dir.split("/")
  for i in range(len(splited_starting_dir)):
    if splited_starting_dir[i] == "phenix_regression":
      print "\tthis_is_test = True"
      this_is_test = True
      
  os.chdir( starting_dir )
  if (this_is_test == True):
    return this_is_test
  else:
    return results
# end of step_8 (cryo_fit itself) function

def step_final(logfile, command_path, starting_dir, model_file_without_pathways, cryo_fit_path):
  os.chdir( starting_dir )
  time_start = time.time()
  show_header("Step 9 (final): Arrange output")
  remake_and_move_to_this_folder(starting_dir, "output")
  
  this_is_test = False # just initial value assignment
  splited_starting_dir = starting_dir.split("/")
  
  cp_command_string = "cp " + command_path + "steps/10_after_cryo_fit/*.py ."
  libtbx.easy_run.fully_buffered(cp_command_string)
  
  logfile.close() # to write user's cc for now
  
  if (model_file_without_pathways != "regression_pdb5khe.pdb"):
    for i in range(len(splited_starting_dir)):
      if splited_starting_dir[i] == "phenix_regression":
        this_is_test = True
        cp_command_string = "cp ../data/input_for_step_final/* ."
        libtbx.easy_run.fully_buffered(cp_command_string)
        make_trajectory_gro()
      
  if (this_is_test == False):
    cp_command_string = "mv ../cc_record_full_renumbered ."
    libtbx.easy_run.fully_buffered(cp_command_string)
  
    cp_command_string = "cp ../steps/8_cryo_fit/cc_record ." # needed for extract_3_highest_cc_gro
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    cp_command_string = "cp ../steps/8_cryo_fit/for_cryo_fit.tpr ."
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    cp_command_string = "cp ../steps/8_cryo_fit/traj.xtc ."
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    print "\n\tExtract .gro files from the 3 highest cc values."
    command_string = "python extract_3_highest_cc_gro.py " + str(this_is_test) + " " + str(cryo_fit_path)
    print "\t\tcommand: ", command_string
    libtbx.easy_run.call(command_string)
  
    os.remove("cc_record")

    print "\n\tConvert .gro -> .pdb"
    print "\t\t(.gro file is for chimera/gromacs/vmd)"
    print "\t\t(.pdb file is for chimera/pymol/vmd)"
    
    for extracted_gro in glob.glob("*.gro"):
      command_string = cryo_fit_path + "editconf -f " + extracted_gro + " -o " + extracted_gro[:-4] + ".pdb"
      print "\t\tcommand: ", command_string
      libtbx.easy_run.fully_buffered(command_string)
    
    make_trajectory_gro(cryo_fit_path)
  
  os.mkdir("trajectory")
  run_this = "mv for_cryo_fit.tpr trajectory.gro traj.xtc trajectory"
  print "\n\tMove trajectory files into trajectory folder "
  print "\t\tcommand: ", run_this
  libtbx.easy_run.call(run_this)
      
  pdb_file_with_original_chains = ''
  for pdb_with_original_chains in glob.glob("../steps/1_make_gro/*.pdb"):
    pdb_file_with_original_chains = pdb_with_original_chains
  
  log_file_name = "../cryo_fit.overall_log"
  logfile = open(log_file_name, "a+") # append
  
  if (this_is_test == False): # recover chain information
    print "\n\tRecover chain information (since gromacs erased it). "
    for pdb in glob.glob("*.pdb"):
      command_string = "python recover_chain.py " + pdb_file_with_original_chains + " " + pdb # worked perfectly with tRNA and Dieter's molecule
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
  
  print "\n\tChange OC1 and OC2 so that molprobity can run"
  for pdb in glob.glob("*.pdb"):
    run_this = "python clean_pdb_for_molprobity.py " + pdb
    print "\t\tcommand: ", run_this
    libtbx.easy_run.call(run_this)
  
  print "\n\tAdd element to each line"
  for pdb in glob.glob("*.pdb"):
    run_this = "python add_element_to_pdb.py " + pdb
    print "\t\tcommand: ", run_this
    libtbx.easy_run.call(run_this)
  
  returned = ''
  if (this_is_test == True):
    returned = check_whether_the_step_was_successfully_ran("Step final", "cryo_fitted_chain_recovered_cleaned_for_molprobity_new_element_added.pdb")
  else:    
    returned = check_whether_the_step_was_successfully_ran("Step final", "cc_record_full_renumbered")
  
  print "\n\tOutputs are in \"output\" folder"
  print "  \t\tIf cryo_fit fitted better than a user provided atomic model, a model with the highest cc value is cryo_fitted.pdb"
  print "  \t\tThis best fitted bio-molecule may not necessarily be the \"best\" atomic model depending on user's specific purposes."
  print "  \t\tTherefore, a user may use other slightly less fitted extracted_x_steps_x_ps.gro/pdb as well."
  print "\n\t\tTo draw a figure for cc,"
  print "  \t\t\tpython <phenix_path>/modules/cryo_fit/steps/9_draw_cc_commandline/draw_cc.py cc_record"
  print "  \t\t\t(Phenix GUI shows this figure automatically)."
  print "  \t\t\t(If a user is using ssh linked linux, set DISPLAY to avoid \"Unable to access the X Display, is $DISPLAY set properly?\")"
  print "\n\t\tTo see a trajectory movie (a user can open a map at the same time of course),"
  print "  \t\t\t<UCSF Chimera 1.13 or later> Tools -> MD/Ensemble analysis -> MD Movie -> Trajectory format = GROMACS,"
  print "  \t\t\t\t.tpr = (trajectory/for_cryo_fit.tpr), .xtc = (trajectory/traj.xtc) -> OK -> (click play button)"
  print "  \t\t\t<VMD 1.9.3 or later> File -> New Molecule -> Browse -> (trajectory/trajectory.gro) -> Load "

  if (returned != "success"):
    write_this = "Step final (arrange output) didn't run successfully"
    print write_this
    logfile.write("Step final (arrange output) didn't run successfully\n")
    logfile.close()
    exit(1)
  
  for py in glob.glob("*.py"): # most users will not need *.py
    run_this = "rm " + py
    libtbx.easy_run.call(run_this)
    
  logfile.write("Step final (arrange output) is successfully ran\n")
  time_end = time.time()
  print "\nStep final", (show_time(time_start, time_end))
  
  os.chdir( starting_dir )
  return this_is_test
# end of step_final (arrange output) function

''' not used now, but keep
def step_9(command_path, starting_dir):
  show_header("Step 9: Show Correlation Coefficient")
  remake_and_move_to_this_folder(starting_dir, "steps/9_draw_cc_commandline")
  
  command_string = "cp " + command_path + "steps/9_draw_cc_commandline/draw_cc.py ."
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
    cryo_fit_path = params.cryo_fit.Input.cryo_fit_path
    cryo_fit_path = cryo_fit_path + "/" # for later steps
  else:
    cryo_fit_path = mdrun_path
  print "cryo_fit_path:",cryo_fit_path
  
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
  bool_step_9 = params.cryo_fit.Steps.step_9
  
  returned = assign_model_name(params, starting_dir, inputs, params.cryo_fit.Input.model_file_name)
  model_file_with_pathways = returned[0]
  model_file_without_pathways = returned[1]
  
  # Options  
  constraint_algorithm_minimization = params.cryo_fit.Options.constraint_algorithm_minimization
  emsteps = params.cryo_fit.Options.emsteps
  emweight_multiply_by = params.cryo_fit.Options.emweight_multiply_by
  emwritefrequency = params.cryo_fit.Options.emwritefrequency
  time_step_for_cryo_fit = params.cryo_fit.Options.time_step_for_cryo_fit
  time_step_for_minimization = params.cryo_fit.Options.time_step_for_minimization
  user_entered_number_of_steps_for_cryo_fit = params.cryo_fit.Options.number_of_steps_for_cryo_fit
  user_entered_number_of_steps_for_minimization = params.cryo_fit.Options.number_of_steps_for_minimization
  
  # Development options
  devel = params.cryo_fit.devel
  no_rerun = params.cryo_fit.no_rerun
  force_field = params.cryo_fit.force_field
  ignh = params.cryo_fit.ignh
  kill_mdrun_mpirun_in_linux = params.cryo_fit.kill_mdrun_mpirun_in_linux
  lincs_order = params.cryo_fit.lincs_order
  missing = params.cryo_fit.missing
  ns_type = params.cryo_fit.ns_type
  number_of_cores_to_use = params.cryo_fit.number_of_cores_to_use
  perturb_xyz_by = params.cryo_fit.perturb_xyz_by
  remove_metals = params.cryo_fit.remove_metals
  
  number_of_steps_for_minimization = determine_number_of_steps_for_minimization(model_file_without_pathways,\
                                                                            model_file_with_pathways, \
                                                                            user_entered_number_of_steps_for_minimization, devel)
  params.cryo_fit.Options.number_of_steps_for_minimization = number_of_steps_for_minimization
  print "\tparams.cryo_fit.Options.number_of_steps_for_minimization (a real value that will be used eventually): ", \
    params.cryo_fit.Options.number_of_steps_for_minimization
  
  number_of_steps_for_cryo_fit = determine_number_of_steps_for_cryo_fit(model_file_without_pathways,\
                                                                            model_file_with_pathways, \
                                                                            user_entered_number_of_steps_for_cryo_fit, devel)
  params.cryo_fit.Options.number_of_steps_for_cryo_fit = number_of_steps_for_cryo_fit
  print "\tparams.cryo_fit.Options.number_of_steps_for_cryo_fit (a real value that will be used eventually): ", \
    params.cryo_fit.Options.number_of_steps_for_cryo_fit
  
  steps_list = [bool_step_1, bool_step_2, bool_step_3, bool_step_4, bool_step_5, bool_step_6, bool_step_7\
                , bool_step_8, bool_step_9]
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
  this_is_test = False # by default
  
  if ((platform.system() == "Linux") and (kill_mdrun_mpirun_in_linux == True)):
    kill_mdrun_mpirun_in_linux()    
  if (steps_list[0] == True):
    this_is_test = step_1(logfile, command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, force_field, ignh, missing, remove_metals, cryo_fit_path)
    logfile.write("Step 1 (Make gro and topology file by regular gromacs) is successfully ran\n")
    
  if (steps_list[1] == True):
    this_is_test = step_2(command_path, starting_dir, model_file_with_pathways, model_file_without_pathways, force_field, \
           perturb_xyz_by, remove_metals)
    logfile.write("Step 2 (Clean gro file to be compatible for amber03 forcefield) is successfully ran\n")
  
  if str(constraint_algorithm_minimization) != "none_default": # this is default for "regression"
    if (steps_list[2] == True):
      this_is_test = step_3(logfile, command_path, starting_dir, ns_type, constraint_algorithm_minimization, number_of_steps_for_minimization, \
             time_step_for_minimization, model_file_without_pathways, devel, cryo_fit_path)
      logfile.write("Step 3 (Make a tpr file for minimization) is successfully ran\n")
    if (steps_list[3] == True):
      this_is_test = step_4(command_path, starting_dir, ns_type, number_of_available_cores, \
                            number_of_cores_to_use, model_file_without_pathways, cryo_fit_path)
      logfile.write("Step 4 (Minimize a gro file (to prevent \"blowup\" during MD Simulation)) is successfully ran\n")
  else: #str(constraint_algorithm_minimization) = "none_default"
    if (steps_list[2] == True):
      step_3(logfile, command_path, starting_dir, ns_type, "none", number_of_steps_for_minimization, \
             time_step_for_minimization, model_file_without_pathways, devel, cryo_fit_path)
      logfile.write("Step 3 (Make a tpr file for minimization) is successfully ran\n")
      
    if (steps_list[3] == True):
      step_4(command_path, starting_dir, ns_type, number_of_available_cores, number_of_cores_to_use, \
             model_file_without_pathways, cryo_fit_path)
      logfile.write("Step 4 (Minimize a gro file (to prevent \"blowup\" during MD Simulation)) is successfully ran\n")
      
    cp_command_string = "cp steps/4_minimize/minimized_c_term_renamed_by_resnum_oc.gro . "
    libtbx.easy_run.fully_buffered(cp_command_string)
    
    shutil.rmtree("steps/3_make_tpr_to_minimize")
    shutil.rmtree("steps/4_minimize")
    
    if (steps_list[2] == True):
      step_3(logfile, command_path, starting_dir, ns_type, "none_default", number_of_steps_for_minimization, \
             time_step_for_minimization, model_file_without_pathways, devel,cryo_fit_path)
      logfile.write("Step 3 (Make a tpr file for minimization) is successfully ran\n")
      
    if (steps_list[3] == True):
      step_4(command_path, starting_dir, ns_type, number_of_available_cores, number_of_cores_to_use, \
             model_file_without_pathways, cryo_fit_path)
      logfile.write("Step 4 (Minimize a gro file (to prevent \"blowup\" during MD Simulation)) is successfully ran\n")
  
  if (steps_list[4] == True):
    this_is_test = step_5(command_path, starting_dir, model_file_without_pathways, cryo_fit_path)
    logfile.write("Step 5 (Make contact potential (constraints) and topology file with it) is successfully ran\n")
  
  if (steps_list[5] == True):
    this_is_test = step_6(command_path, starting_dir, model_file_without_pathways)
    logfile.write("Step 6 (Make all charges of atoms be 0) is successfully ran\n")
  
  if (this_is_test == True):
    return 0 # return early for regression of step 1~6
  
  cc_has_been_increased = True # just an initial value
  charge_group_moved = True # just an initial value
  re_run_with_higher_map_weight = True # just an initial value
  
  returned = assign_map_name(params, starting_dir, inputs, params.cryo_fit.Input.map_file_name)
  map_file_with_pathways = returned[0]
  map_file_without_pathways = returned[1]
  
  if (model_file_without_pathways == "tRNA_tutorial.pdb"):
    no_rerun = True
    number_of_steps_for_cryo_fit = 5000
    #number_of_steps_for_cryo_fit = 100 # temp for devel
  
  if (devel == True):
    no_rerun = True
    number_of_steps_for_cryo_fit = 100
    
  restart_w_longer_steps = False # this is a proper initial assignment
  re_run_with_higher_map_weight = False # this is a proper initial assignment
  
  while ((cc_has_been_increased == True) or (charge_group_moved == True) or (re_run_with_higher_map_weight == True)):
    if ((this_is_test == True) or (steps_list[0] == False and steps_list[1] == False and steps_list[2] == False \
                                and steps_list[3] == False and steps_list[4] == False and steps_list[5] == False \
                                and steps_list[6] == False and steps_list[7] == False)):
      break
    
    if (steps_list[6] == True):
      this_is_test = step_7(command_path, starting_dir, number_of_steps_for_cryo_fit, emweight_multiply_by, emsteps, \
             emwritefrequency, lincs_order, time_step_for_cryo_fit, model_file_without_pathways, cryo_fit_path)
      logfile.write("Step 7 (Make a tpr file for cryo_fit) is successfully ran\n")
      if (this_is_test == True):
        return 0 # return early for regression of step 7
    
    if (steps_list[7] == True):
      results = step_8(logfile, command_path, starting_dir, number_of_available_cores, number_of_cores_to_use, 
             map_file_with_pathways, no_rerun, devel, restart_w_longer_steps, \
             re_run_with_higher_map_weight, model_file_without_pathways, cryo_fit_path)
      
      if (results == True): # this is a test
        print "This is a test, so break early"
        break  
      
      if (model_file_without_pathways == "regression.pdb"): # for regression purpose
        logfile.write("Step 8 (cryo_fit itself) is successfully ran\n")
        return results
      
      if results == "failed":
        return "failed" # flat failed
      
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
          write_this = "re_run_with_longer_steps is recommended, but no_rerun = True, Step 8 (cryo_fit itself) is successfully ran\n"
          print write_this
          logfile.write(write_this)
          this_is_test = step_final(logfile, command_path, starting_dir, model_file_without_pathways) # just to arrange final output
          exit(1)
          #return results
        restart_w_longer_steps = True
        re_run_with_higher_map_weight = False
        
        # copy for a next restart step
        if (os.path.isfile("state.cpt") == False):
          write_this = "state.cpt not found, step_8 may be full of stepxb_nx.pdb. This means that map weight may be too high, exit now\n"
          print write_this
          logfile.write(write_this)
          exit(1)
        
        cp_command_string = "cp state.cpt ../.."
        libtbx.easy_run.fully_buffered(command=cp_command_string).raise_if_errors()
  
        charge_group_moved = False # just initial value
        number_of_steps_for_cryo_fit = number_of_steps_for_cryo_fit * 2
        
        write_this = "\nStep 8 (cryo_fit itself) is ran well, but correlation coefficient values tend to be increased over the last 30 steps\n"
        print write_this
        logfile.write(write_this)
        
        write_this = "Therefore, step 7 & 8 will re-run with longer steps (" + str(number_of_steps_for_cryo_fit) + ")\n\n"
        print write_this
        logfile.write(write_this)
        
        restart_record = open("../../restart_record.txt", "a+")
        write_this = str(number_of_steps_for_cryo_fit)+"\n"
        restart_record.write(write_this)
        restart_record.close()
        
        if (number_of_steps_for_cryo_fit > 1000000000000000 ): # to avoid infinite loop
          write_this = "number_of_steps_for_cryo_fit > 1000000000000000, exit now"
          print write_this
          logfile.write(write_this)
          break # no need to use exit(1) since this break will break this while loop
        os.chdir( starting_dir ) # needed for re-running
      
      elif results == "re_run_with_higher_map_weight":
        emweight_multiply_by = emweight_multiply_by * 2
        write_this = "\nStep 8 (cryo_fit itself) is ran well, but correlation coefficient values tend to decrease over the last 30 steps\n"
        print write_this
        logfile.write(write_this)
        write_this = "Therefore, step 7 & 8 will re-run with higher emweight_multiply_by (e.g. " + str(emweight_multiply_by) + ")\n\n"
        print write_this
        logfile.write(write_this)
        re_run_with_higher_map_weight = True
        if (emweight_multiply_by > 1000000000000000 ): # to avoid infinite loop
          write_this = "emweight_multiply_by > 1000000000000000, exit now"
          print write_this
          logfile.write(write_this)
          break # no need to use exit(1) since this break will break this while loop
        os.chdir( starting_dir ) # needed for re-running
        
      else: # normal ending of cryo_fit
        charge_group_moved = False
        cc_has_been_increased = False
        re_run_with_higher_map_weight = False
  logfile.write("Step 8 (cryo_fit itself) is successfully ran\n")
  this_is_test = step_final(logfile, command_path, starting_dir, model_file_without_pathways, \
                            cryo_fit_path) # just to arrange final output
  if (this_is_test == False):
    return results
  
  # keep for now for this cc draw
  #if (steps_list[8] == True):
  #  step_9(command_path, starting_dir, model_file_without_pathways, map_file_without_pathways)
# end of run_cryo_fit function

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
  
  return results
# end of cmd_run function

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
