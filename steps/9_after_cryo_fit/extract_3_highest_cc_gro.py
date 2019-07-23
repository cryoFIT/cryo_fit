import os, subprocess, sys

# It is ESSENTIAL to adjust step number if restarted
def adjust_step_number():
    print "\n\t\t\tAdjust step number due to restart."
    f = open('../restart_record.txt', 'r')
    # count # of lines
    number_of_lines = 0
    for line in f:
        number_of_lines = number_of_lines + 1
    f.close()
    
    f = open('../restart_record.txt', 'r')
    j = 0
    last_step_to_be_added = '' # initial
    for line in f:
        last_step_to_be_added = int(line) # just in case when there is only one value like nucleosome case
        j = j + 1
        if (j == (number_of_lines - 1)):
            last_step_to_be_added = int(line)
            break
    f.close()
    print "\t\t\t\tAdd this step number to each current step number:", last_step_to_be_added
    
    f_in = open('cc_record', 'r')
    f_out = open('cc_record_adjusted_step_use_for_extraction', 'w')
    for line in f_in:
       splited = line.split()
       new_line = splited[0] + " " + str((int(last_step_to_be_added)+int(splited[1]))) + " " + splited[2] + " " + splited[3] + " " + splited[4] + "\n"
       f_out.write(new_line)
    f_in.close()
    f_out.close()
    
    os.remove("cc_record") # no longer neeeded
################# end of def adjust_step_number ()


def extract_gro(target_step, i, cc, cryo_fit_path):
    for_cryo_fit_mdp_location = ''
    if (this_is_test == "False"):
        for_cryo_fit_mdp_location = "../steps/7_make_tpr_with_disre2/for_cryo_fit.mdp"
    else:
        print "\t This is a test in extract_gro"
        for_cryo_fit_mdp_location = "for_cryo_fit.mdp"
    
    grep_dt_string = "grep dt " + for_cryo_fit_mdp_location + " | grep -v when"
    
    print "\t\t\t\tcommand:", grep_dt_string
    result = os.popen(grep_dt_string).read()
    splited = result.split()
    dt = splited[2]

    print_this = "\t\t\t\tdt:" + dt + "\n"
    print print_this
    
    grep_nsteps_string = "grep nsteps " + for_cryo_fit_mdp_location + " | grep -v when"
    result = os.popen(grep_nsteps_string).read()
    splited = result.split()
    nsteps = splited[2]
    print_this = "\t\t\t\tnsteps:" + str(nsteps) + "\n"
    print print_this
    
    total_ps = float(dt)*float(nsteps)

    print_this = "\t\t\t\ttotal_ps = float(dt)*float(nsteps) = " + str(total_ps) + "\n"
    print print_this
    
    ''' reading_frame is relevant only when restarted
    print_this = "\t\t\t\tEstimated reading_frame = total_ps/((dt)*1,000) = " + str(float(total_ps)/((float(dt))*1000.0)) + " ps \n"
    print print_this
    f_out.write(print_this)
    '''
    
    print_this = "\t\t\t\tTherefore, total mdrun running time was: " + str(total_ps) + " pico (10^-12) second" + "\n"
    print print_this
    
    print_this = "\t\t\t\tCryo_fit needs to extract a gro file from " + str(target_step) + " step(s)" + "\n"
    print print_this
        
    print_this = "\t\t\t\ttarget_ps = (float(target_step)/float(nsteps))*float(total_ps)" + "\n"
    print print_this

    target_ps = (float(target_step)/float(nsteps))*float(total_ps)
    print_this = "\t\t\t\tTherefore, the cryo_fit will extract a gro file from " + str(target_ps) + " ps" + "\n"
    print print_this
    
    output_gro_name = "extracted_" + str(target_step) + "_steps_" + str(target_ps) + "_ps.gro"
    
    os.system("echo 0 > input_parameters") # to select system
    
    cmd = cryo_fit_path + "trjconv -f traj.xtc -dump " + str(target_ps) + " -o " + str(output_gro_name) + \
          " -s for_cryo_fit.tpr < input_parameters"
    print "\t\t\t\tcommand: ",cmd
    os.system(cmd)
    
    if (i == 0):
        print "\t\t\t\t", target_step, " step has the highest cc"
        if (target_step == "0"): # works as expected
           print "\t\t\t\tHowever, it was the initial model that a user provided, so don't rename it to cryo_fitted.gro"
           cmd = "mv " + output_gro_name + " user_provided.gro"
           print "\t\t\t\t\tcommand:", cmd, "\n"
           os.system(cmd)
        else:
            users_cc = get_users_cc_from_overall_log("../cryo_fit.overall_log")
            print "cc:",cc
            print "users_cc:",users_cc
            
            if (float(cc) > float(users_cc)):
                print "\t\t\t\tso rename it to cryo_fitted.gro"
                cmd = "mv " + output_gro_name + " cryo_fitted.gro"
                print "\t\t\t\t\tcommand:", cmd, "\n"
                os.system(cmd)
    os.remove("input_parameters")
################# end of extract_gro function


def get_users_cc_from_overall_log(log):
  f_in = open(log)
  for line in f_in:
    
    '''
    splited_by_apostrophe = line.split("'")
    if (splited_by_apostrophe[0] == "User"):
        splited = line.split(" ")
        cc = splited[5]
        f_in.close()
        print "\tUser provided atomic model's cc: ", cc
        return cc
    '''
    splited = line.split(" ")
    if (splited[0] == "A"):
        if (splited[1] == "user's"):
            cc = splited[7]
            f_in.close()
            print "\tUser provided atomic model's cc: ", cc
            return cc
    
################# end of get_users_cc(cc_record)


if (__name__ == "__main__") :
    args=sys.argv[1:]
    this_is_test = args[0]
    cryo_fit_path = args[1]
    no_rerun = args[2]

    print "\n\t\tExtract 3 highest cc gro (among the whole run, not just the last run)"
    # Although I assign number_of_steps_for_cryo_fit*2 as a new number_of_steps_for_cryo_fit,
    #due to state.cpt, mdrun runs only until a new number_of_steps_for_cryo_fit INCLUDING FORMERLY RAN STEPS  
    
    # The last run probably has the highest CC anyway, so let me extract only in the last run.
    # traj.xtc is overwritten every time when em_weight or number_of_steps_for_cryo_fit is reassigned.
    # Actually previous traj.xtc is erased (not keeping previous record) every time when em_weight or number_of_steps_for_cryo_fit is reassigned.
    # --> so cc_record_full_renumbered should NOT be used for extrqcting gro, should be used only for overall cc change

    result = '' # initial temporary assignment
    if (this_is_test == "False"): # default running
        # adjust step number if I restarted
        if (os.path.isfile("../restart_record.txt") == True):
            adjust_step_number ()
            os.remove("../restart_record.txt") # only for development, keep this file
        
        if (no_rerun == "False"): # default running
            # this cc_record is step_adjusted if restarted
            if (os.path.isfile("cc_record_adjusted_step_use_for_extraction") == False):
                print "cc_record_adjusted_step_use_for_extraction is not found, please email doonam@lanl.gov"
                exit(1)
            result = os.popen("cat cc_record_adjusted_step_use_for_extraction | sort -nk5 -r | head -3").readlines()
        else:
            result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    else: # test
        result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    print "3 highest cc steps that need to be extracted:", result
    
    
    for i in range(len(result)):
        splited = result[i].split()
        target_step = splited[1]
        cc = splited[4]
        print "\t\t\ttarget_step for extracting a gro file: ", target_step
        extract_gro(target_step, i, cc, cryo_fit_path)
