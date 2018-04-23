from elbow.utilities.phil_utils import master_phil
import glob, iotbx.pdb.hierarchy, os, platform, subprocess, sys, time
from iotbx import file_reader
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out
import mmtbx.utils
from subprocess import check_output, Popen, PIPE

args=sys.argv[1:]
command_path = args[0]
common_functions_path = command_path + "/common_functions/"
sys.path.insert(0, common_functions_path)
from common_functions import *

ns_type = args[1]
number_of_available_cores = int(args[2])
number_of_cores_to_use = args[3] # for mpi -> cores, for threads -> threads
target_map_with_pathways = args[4]
starting_dir = args[5]
output_file_name_prefix = args[6]
this_is_test = args[7]

if (__name__ == "__main__") :
  #remove_former_files() #needed for development only
  
  cp_command_string = ''
  if (str(this_is_test) == "False"):
    cp_command_string = "cp ../7_make_tpr_with_disre2/for_cryo_fit.tpr . "
  else:
    cp_command_string = "cp ../../data/input_for_step_8/* ."

  libtbx.easy_run.fully_buffered(command=cp_command_string).raise_if_errors()
  
  bool_minimization = False
  home_bin_cryo_fit_bin_dir = know_home_cryo_fit_bin_dir_by_ls_find()
  
  bool_just_get_input_command = False
  write_this_input_command, output_file_name = first_prepare_for_minimization_cryo_fit(bool_minimization, bool_just_get_input_command, \
                                                         home_bin_cryo_fit_bin_dir, \
                                                         ns_type, number_of_available_cores, \
                                                         number_of_cores_to_use,\
                                                         target_map_with_pathways, output_file_name_prefix)
  
  f_out = open('log.step_8_cryo_fit_used_command', 'wt')
  f_out.write(write_this_input_command)
  f_out.close()
#end of if (__name__ == "__main__")
