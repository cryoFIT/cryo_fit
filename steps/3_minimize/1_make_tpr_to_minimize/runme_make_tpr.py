import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from iotbx import file_reader
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out

''' # may not work at Karissa's macbook?
args=sys.argv[1:]
command_path = args[0]
common_functions_path = command_path + "/command_line/"
sys.path.insert(0, common_functions_path)
from common_functions import *
'''

cryo_fit_repository_dir = libtbx.env.dist_path("cryoFIT")
common_functions_path = cryo_fit_repository_dir + "/common_functions/"
print "\tcommon_functions_path:", common_functions_path

sys.path.insert(0, common_functions_path)
from common_functions import *

home_cryo_fit_bin_dir = know_home_cryo_fit_bin_dir_by_ls_find()
run_this = home_cryo_fit_bin_dir + "/grompp -f minimization.mdp -c *.gro -p *.top -o to_minimize.tpr \
           -maxwarn 10"
  # -f, -c, -p are for input files of grompp
  # -o is for output file
print "\tcommand: ", run_this

f_out = open('log.step_3_1_make_tpr', 'wt')
write_this_input_command = home_cryo_fit_bin_dir + "\n"
write_this_input_command = write_this_input_command + run_this + "\n"
f_out.write(write_this_input_command)

time_start = time.time()
libtbx.easy_run.call(command=run_this)
time_end = time.time()

#write_this_time = show_time("step_3_1_make_tpr", time_start, time_end)
write_this_time = "step_3_1_make_tpr"
write_this_time = write_this_time + show_time (time_start, time_end)
write_this_time = "\n\n" + write_this_time + "\n"
f_out.write(write_this_time)
f_out.close()
