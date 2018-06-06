# This code runs @steps/8_cryo_fit
import os, subprocess, sys
from os.path import expanduser # to find home_dir
args=sys.argv[1:]
this_is_test = args[0]

def adjust_step_number():
    print "\n\t\t\tAdjust step number due to restart."
    f = open('../../restart_record.txt', 'r')
    # count # of lines
    i = 0
    for line in f:
        i = i + 1
    #print "\t\t\t\t../../restart_record.txt has ", i, " number of lines"
    f.close()
    
    f = open('../../restart_record.txt', 'r')
    j = 0
    last_step_to_be_added = '' # initial
    for line in f:
        j = j + 1
        if (j == (i - 1)):
            last_step_to_be_added = int(line)
            break
    f.close()
    print "\t\t\t\tAdd this step number to each current step number:", last_step_to_be_added
    
    f_in = open('cc_record', 'r')
    f_out = open('cc_record_adjusted_step', 'w')
    for line in f_in:
       splited = line.split()
       #print "\t\t\t\told line:",line
       new_line = splited[0] + " " + str((int(last_step_to_be_added)+int(splited[1]))) + " " + splited[2] + " " + splited[3] + " " + splited[4] + "\n"
       #print "\t\t\t\tnew_line:",new_line
       f_out.write(new_line)
    f_in.close()
    f_out.close()
    command_string = "mv cc_record_adjusted_step cc_record"
    os.system(command_string)
# end of def adjust_step_number ()
    
def extract_gro(target_step, i):
    for_cryo_fit_mdp_location = ''
    
    if this_is_test == "False":
        for_cryo_fit_mdp_location = "../7_make_tpr_with_disre2/for_cryo_fit.mdp"
    else:
        print "\tthis_is_test in extract_gro:", this_is_test
        for_cryo_fit_mdp_location = "for_cryo_fit.mdp"
        
    grep_dt_string = "grep dt " + for_cryo_fit_mdp_location + " | grep -v when"
    print "\t\t\t\tcommand:", grep_dt_string
    result = os.popen(grep_dt_string).read()
    splited = result.split()
    dt = splited[2]
    
    grep_nsteps_string = "grep nsteps " + for_cryo_fit_mdp_location + " | grep -v when"
    result = os.popen(grep_nsteps_string).read()
    splited = result.split()
    nsteps = splited[2]
    
    total_ps = float(dt)*float(nsteps)
    print "\t\t\t\ttotal_ps = float(dt)*float(nsteps) = ", total_ps
    
    print "\t\t\t\tTherefore, total mdrun running time was: ", total_ps, "pico (10^-12) second"
    print "\t\t\t\tCryo_fit needs to extract a gro file from ", target_step, "steps"
    target_ps = (float(target_step)/float(nsteps))*float(total_ps)
    print "\t\t\t\ttarget_ps = (float(target_step)/float(nsteps))*float(total_ps)"
    print "\t\t\t\tTherefore, the cryo_fit will extract a gro file from ", target_ps, "ps"
    
    output_gro_name = "extracted_" + str(target_step) + "_steps_" + str(target_ps) + "_ps.gro"
    os.system("echo 0 > input_parameters") # to select system
    
    home_dir = expanduser("~")
    home_bin_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit/bin"
    cmd = home_bin_cryo_fit_bin_dir + "/trjconv -f traj.xtc -dump " + str(target_ps) + " -o " + str(output_gro_name) + \
          " -s for_cryo_fit.tpr < input_parameters"
    os.system(cmd)
    
    if (i == 0):
        print "\t\t\t\t", target_step, " step has the highest cc"
        if (target_step == "0"): # works as expected
            print "\t\t\t\tHowever, it was the initial model that a user provided, so don't rename it to cryo_fitted.gro"
            cmd = "mv " + output_gro_name + " user_provided.gro"
            print "\t\t\t\t\tcommand:", cmd, "\n"
            os.system(cmd)
        else:
            print "\t\t\t\tso rename it to cryo_fitted.gro"
            cmd = "mv " + output_gro_name + " cryo_fitted.gro"
            print "\t\t\t\t\tcommand:", cmd, "\n"
            os.system(cmd)
# end of extract_gro function

if (__name__ == "__main__") :
    print "\n\t\textract_3_highest_cc_gro_from_cryofit_md_log"
    cmd = "grep correlation md.log > cc_record"
    print "\t\t\tcommand:", cmd
    os.system(cmd)
    
    # adjust step number if needed
    #print "os.path.isfile(\"../../restart_record.txt\"):",os.path.isfile("../../restart_record.txt")
    if (os.path.isfile("../../restart_record.txt") == True):
        adjust_step_number ()
        os.remove("../../restart_record.txt") # only for development, keep this file
    
    result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    for i in range(len(result)):
        splited = result[i].split()
        target_step = splited[1]
        print "\t\t\ttarget_step for extracting a gro file: ", target_step
        extract_gro(target_step, i)
