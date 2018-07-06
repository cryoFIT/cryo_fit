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

number_of_available_cores = int(args[1])
number_of_cores_to_use = args[2] # for mpi -> cores, for threads -> threads
target_map_with_pathways = args[3]
starting_dir = args[4]
this_is_test = args[5]
restart = args[6]
cryo_fit_path = args[7]

if (__name__ == "__main__") :
  
  cp_command_string = ''
  if (str(this_is_test) == "False"):
    cp_command_string = "cp ../7_make_tpr_with_disre2/for_cryo_fit.tpr . "
  else:
    cp_command_string = "cp ../../data/input_for_step_8/* ."
  libtbx.easy_run.fully_buffered(command=cp_command_string).raise_if_errors()
  
  write_this_input_command = first_prepare_cryo_fit(number_of_available_cores, \
                                                         number_of_cores_to_use, \
                                                         target_map_with_pathways, restart, cryo_fit_path)
  
  f_out = open('log.step_8_cryo_fit_used_command', 'wt')
  f_out.write(write_this_input_command)
  f_out.close()
#end of if (__name__ == "__main__")
