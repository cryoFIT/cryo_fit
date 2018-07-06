import glob, os, subprocess, sys, time
from subprocess import check_output

args=sys.argv[1:]
command_path = args[0]
cryo_fit_path = args[1]
common_functions_path = command_path + "/common_functions/"
sys.path.insert(0, common_functions_path)
from common_functions import *

command_string = cryo_fit_path + "grompp -f for_cryo_fit.mdp -c *.gro -p *0_charge.top \
                 -o for_cryo_fit.tpr -maxwarn 10"       
                  # -f, -c, -p are for input files of grompp
                  # -o is for an output file
                  
f_out = open('log.step_7_used_command', 'wt')
write_this_input_command = command_string + "\n"
f_out.write(write_this_input_command)

f_out = open('log.step_7_used_command', 'at+') # reopen here
time_start = time.time()
os.system(command_string)
time_end = time.time()

write_this = "\n" + show_time("step_7", time_start, time_end) + "\n"
color_print ((write_this), 'green')

f_out.close()
