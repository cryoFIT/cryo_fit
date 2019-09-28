import glob, os, subprocess, sys, time
from os.path import expanduser # to find home_dir
from subprocess import check_output

''' needed only for development
# don't use this remove_former_files in common funtion, step 2-1 is unique
def remove_former_files():
  current_directory = os.getcwd()
  print "\tRemove former files in ", current_directory
  for each_file in glob.glob("*"):
    if (each_file[:1] == "#") or (each_file[-1:] == "~") or (each_file[-4:] == ".gro") \
                 or (each_file == "inp_temp") or (each_file[-4:] == ".itp") or (each_file == "md.out") or (each_file[-4:] == ".top"):
      subprocess.call(["rm", each_file])  
'''

def make_gro_top(input_pdb_file_name, force_field, *args):
  print "\n\tMake gro and topology files from a given pdb file." 
  
  splited = input_pdb_file_name.split("/")
  input_pdb_file_name_without_folder = splited[len(splited)-1]

  output_gro_file_name_without_folder = input_pdb_file_name_without_folder[:-4] + "_by_pdb2gmx.gro"
  output_top_file_name_without_folder = input_pdb_file_name_without_folder[:-4] + "_by_pdb2gmx.top"
  
  if ((force_field == "amber03") or (force_field == "a")):
      run_this = "echo 1 > input_parameters" # amber03 forcefield (actually amber03_gdp_MIA.ff)
      libtbx.easy_run.call(run_this)
      
      run_this = "echo 6 >> input_parameters" # no water model
      libtbx.easy_run.call(run_this)
  else:
      os.system("echo 13 > input_parameters") # gromos96 53a6 forcefield
      os.system("echo 3 >> input_parameters") # no water model
  
  print "\tinput pdb file name: ", input_pdb_file_name
  
  f_out = open('log.step_1_2_runme_make_gro', 'wt')
  run_this = '' # initial
  
  common_command_script = cryo_fit_path + "pdb2gmx -f " + input_pdb_file_name_without_folder + \
                          " -o " + output_gro_file_name_without_folder + " -p " + \
                          output_top_file_name_without_folder  + " -merge all  "
  
   # gmx pdb2gmx reads a .pdb (or .gro) file, reads some database files, adds hydrogens to the molecules and \
   # generates coordinates in GROMACS (GROMOS), or optionally .pdb, format and a topology in GROMACS format. \
   # These files can subsequently be processed to generate a run input file.
  if (bool_missing == "True"): # default choice
    common_command_script = common_command_script + "-missing"
  else:
    common_command_script = common_command_script + "-nomissing"
  if (bool_ignh == "True"):
    run_this = common_command_script + " -ignh < input_parameters >> md.out"
  else:
    run_this = common_command_script + " -noignh < input_parameters >> md.out"
  
  print "\tcommand: ", run_this
  
  time_start = time.time()
  libtbx.easy_run.call(run_this) 
  time_end = time.time()
  
  write_this_input_command = run_this + "\n"
  f_out.write(write_this_input_command)
  
  write_this_time = show_time(time_start, time_end)
  write_this_time = "\n\nstep_1_2_runme_make_gro" + write_this_time + "\n"
  f_out.write(write_this_time)
  
  f_out.close()
########### end of make_gro_top function


if (__name__ == "__main__") :
  #remove_former_files() # only needed for development
  args=sys.argv[1:]
  if len(args) < 5:
    count = 0
    for pdb_file in glob.glob("*.pdb"):
      input_pdb_file_name = pdb_file # if there is only 1 pdb file in this folder, use it
      count +=1
      if count == 2:
        print "Specify one input PDB file"
        print "Example usage: clean_pdb_for_gromacs.py input.pdb"
        sys.exit("clean_pdb_for_gromacs exits now (expecting a pdb file at next run)")
    force_field = args[1] # pdb input file
    bool_ignh = args[2]
    bool_missing = args[3]
    make_gro_top(input_pdb_file_name, force_field, bool_ignh, bool_missing)
  else: # will be used by python
    for pdb_file in glob.glob("*.pdb"):
      input_pdb_file_name = pdb_file # if there is only 1 pdb file in this folder, use it
      command_path = args[0]
      force_field = args[1] # pdb input file
      bool_ignh = args[2]
      bool_missing = args[3]
      cryo_fit_path = args[4]
      
      common_functions_path = command_path + "/common_functions/"
      sys.path.insert(0, common_functions_path)
      from common_functions import *
      make_gro_top(input_pdb_file_name, force_field, bool_ignh, bool_missing, cryo_fit_path)
  bool_gro = 0
  for gro_file in glob.glob("*.gro"):
      bool_gro = 1
      print "\tSuccess. A gro file is made\n"
  if (bool_gro == 0):
      print "\tFailure. A gro file is not made\n"
########## end of if (__name__ == "__main__")