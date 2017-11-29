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
output_file_format = args[5]
starting_dir = args[6]
output_file_name_prefix = args[7]

if (__name__ == "__main__") :
  #remove_former_files() #needed for development only
  command_string = "cp ../6_make_tpr_with_disre2/for_cryo_fit.tpr . "
  color_print ("\tcommand: ", 'green')
  print "\t", command_string, "\n"
  libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors()
  
  f_out = open('log.step_7', 'wt')
  bool_minimization = False
  bool_just_get_input_command = True
  print "\ttarget_map_with_pathways:", target_map_with_pathways
  home_bin_cryo_fit_bin_dir = know_home_cryo_fit_bin_dir_by_ls_find()
  command_that_will_be_used, output_file_name = first_prepare_for_minimization_cryo_fit(bool_minimization, \
                                                              bool_just_get_input_command, \
                                                              home_bin_cryo_fit_bin_dir, \
                                                              ns_type, number_of_available_cores, \
                                                              number_of_cores_to_use,\
                                                              target_map_with_pathways, output_file_format,
                                                              output_file_name_prefix)
  write_this_input_command = str(command_that_will_be_used) + "\n"
  f_out.write("will_be_used_input_command\n")
  f_out.write(write_this_input_command)
  f_out.write("\n")
  
  bool_just_get_input_command = False
  command_used, output_file_name = first_prepare_for_minimization_cryo_fit(bool_minimization, bool_just_get_input_command, \
                                                         home_bin_cryo_fit_bin_dir, \
                                                         ns_type, number_of_available_cores, \
                                                         number_of_cores_to_use,\
                                                         target_map_with_pathways, output_file_format,
                                                         output_file_name_prefix)
  
  write_this_input_command = str(command_used) + "\n"
  f_out.write("really_used_input_command\n")
  f_out.write(write_this_input_command)
  f_out.write("\n")
  f_out.close()
  
  f_out = open('log.step_7_cryoFIT_real_command', 'wt')
  f_out.write(write_this_input_command)
  f_out.close()
#end of if (__name__ == "__main__")
