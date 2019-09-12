import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from iotbx import file_reader
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out
from subprocess import check_output

# some header(s) among these are needed for libtbx.env.dist_path
from cctbx import maptbx
import iotbx.pdb
import iotbx.pdb.mmcif
import mmtbx.model
import mmtbx.utils

cryo_fit_repository_dir = libtbx.env.dist_path("cryo_fit")
common_functions_path = os.path.join(cryo_fit_repository_dir, 'common_functions')
print common_functions_path
sys.path.insert(0, common_functions_path)
from common_functions import *

args=sys.argv[1:]
cryo_fit_path = args[0]

run_this = cryo_fit_path + "grompp -f minimization.mdp -c *.gro -p *.top -o to_minimize.tpr -maxwarn 10"
  # -f, -c, -p are for input files of grompp
  # -o is for output file
print "\tcommand: ", run_this

f_out = open('log.step_3_make_tpr', 'wt')
write_this_input_command = run_this + "\n"
f_out.write(write_this_input_command)
f_out.close()

time_start = time.time()
libtbx.easy_run.call(command=run_this)
time_end = time.time()

f_out = open('log.step_3_make_tpr', 'a+')
write_this_time = "step_3_make_tpr"
write_this_time = write_this_time + show_time (time_start, time_end)
write_this_time = "\n" + write_this_time + "\n"
f_out.write(write_this_time)
f_out.close()
