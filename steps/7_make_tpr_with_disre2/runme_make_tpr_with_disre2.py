import glob, os, subprocess, sys, time
from subprocess import check_output

args=sys.argv[1:]
command_path = args[0]
common_functions_path = command_path + "/common_functions/"
sys.path.insert(0, common_functions_path)
from common_functions import *

#remove_former_files() # only needed for development
home_cryo_fit_bin_dir = know_home_cryo_fit_bin_dir_by_ls_find()
command_string = home_cryo_fit_bin_dir + "/grompp -f for_cryo_fit.mdp -c *.gro -p *0_charge.top \
                 -o for_cryo_fit.tpr -maxwarn 10"       
                  # -f, -c, -p are for input files of grompp
                  # -o is for an output file
                  
f_out = open('log.step_7', 'wt')
write_this_input_command = command_string + "\n"
f_out.write(write_this_input_command)
f_out.close() # close early so that it writes input command for sure, before writing running time

f_out = open('log.step_7', 'at+') # reopen here
time_start = time.time()
os.system(command_string)
time_end = time.time()

write_this = "\n" + show_time("step_7", time_start, time_end) + "\n"
color_print ((write_this), 'green')
f_out.write(write_this)
f_out.close()
