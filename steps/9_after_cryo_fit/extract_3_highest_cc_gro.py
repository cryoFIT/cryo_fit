import os, subprocess, sys

def file_size(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size
######## end of file_size(fname)

''' do not import these, to avoid "extract_3_highest_cc_gro.py:121: UserWarning: os.popen() is not safe: please use the subprocess module or libtbx.easy_run instead."
# this is needed to import all common functions
path = subprocess.check_output(["which", "phenix.cryo_fit"])
splited = path.split("/")
command_path = ''
for i in range(len(splited)-3):
  command_path = command_path + splited[i] + "/"
command_path = command_path + "modules/cryo_fit/"
common_functions_path = command_path + "common_functions/"
sys.path.insert(0, common_functions_path)
from common_functions import *
'''

# This is ESSENTIAL to adjust step number if restarted for longer steps to extract gro file
# Otherwise, there will be "WARNING no output, last frame read at t=xxx".
# Even with cc_record.txt, 
def adjust_step_number():
    print "\n\t\t\tAdjust step number due to restart for longer steps."
    f = open('../restart_record_for_longer_steps.txt', 'r')
    # count # of lines
    number_of_lines = 0
    for line in f:
        number_of_lines = number_of_lines + 1
    f.close()
    
    #last_step_to_be_added = '' # initial
    last_step_to_be_added = 0 # initial
    f = open('../restart_record_for_longer_steps.txt', 'r')
    j = 0
    for line in f:
        last_step_to_be_added = int(line) # just in case when there is only one value like nucleosome case
        j = j + 1
        if (j == (number_of_lines - 1)):
            last_step_to_be_added = int(line)
            break
    f.close()
    print "\t\t\t\tAdd this step number to each current step number in cc_record file:", last_step_to_be_added
    
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


def extract_gro(gro_extraction_note_file, cryo_fit_path, nsteps, total_ps, target_step, i, cc):

    print_this = "\n\ttarget_ps = (float(target_step)/float(nsteps))*float(total_ps)" + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)

    target_ps = (float(target_step)/float(nsteps))*float(total_ps)
    print_this = "\tTherefore, the cryo_fit will extract a gro file from " + str(target_ps) + " ps" + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    output_gro_name = "extracted_" + str(target_step) + "_steps_" + str(target_ps) + "_ps.gro"
    
    os.system("echo 0 > input_parameters") # to select system
    
    cmd = cryo_fit_path + "trjconv -f traj.xtc -dump " + str(target_ps) + " -o " + str(output_gro_name) + \
          " -s for_cryo_fit.tpr < input_parameters"
    write_this = "\t" + cmd + "\n"
    print write_this
    gro_extraction_note_file.write(write_this)
    os.system(cmd)
    
    returned_file_size = file_size(output_gro_name)
    if (returned_file_size == 0):
        write_this = "extracted gro file is empty, check step numbers, cryo_fit will exit soon."
        print write_this
        gro_extraction_note_file.write(write_this)
        gro_extraction_note_file.close()
        return "empty" 
        
    if (i == 0):
        print "\t", target_step, " step has the highest cc"
        if (target_step == "0"): # works as expected
           print "\tHowever, it was the initial model that a user provided, so don't rename it to cryo_fitted.gro"
           cmd = "mv " + output_gro_name + " user_provided.gro"
           write_this = "\t" + cmd + "\n"
           print write_this
           gro_extraction_note_file.write(write_this)
           
           os.system(cmd)
           
        else:
            users_cc = get_users_cc_from_overall_log("../cryo_fit.overall_log")
            # print "cc:",cc
            # print "users_cc:",users_cc
            
            if (float(cc) > float(users_cc)):
                print "\ttherefore rename it to cryo_fitted.gro"
                cmd = "mv " + output_gro_name + " cryo_fitted.gro"
                print "\t", cmd, "\n"
                gro_extraction_note_file.write(cmd)
                os.system(cmd)
                
    os.remove("input_parameters")
    return 1
################# end of extract_gro function



def get_nsteps_total_ps(gro_extraction_note_file, cryo_fit_path):
    for_cryo_fit_mdp_location = ''
    if (this_is_test == "False"):
        for_cryo_fit_mdp_location = "../steps/7_make_tpr_with_disre2/for_cryo_fit.mdp"
    else:
        print "\t This is a test in extract_gro"
        for_cryo_fit_mdp_location = "for_cryo_fit.mdp"
    
    grep_dt_string = "grep dt " + for_cryo_fit_mdp_location + " | grep -v when"
    
    print "\tcommand:", grep_dt_string
    gro_extraction_note_file.write(grep_dt_string)
    result = os.popen(grep_dt_string).read()
    splited = result.split()
    dt = splited[2]

    print_this = "\n\tdt:" + dt + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    grep_nsteps_string = "grep nsteps " + for_cryo_fit_mdp_location + " | grep -v when"
    
    result = os.popen(grep_nsteps_string).read()
    
    splited = result.split()
    nsteps = splited[2]
    print_this = "\tnsteps: " + str(nsteps) + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    total_ps = float(dt)*float(nsteps)

    print_this = "\ttotal_ps = float(dt)*float(nsteps) = " + str(total_ps) + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    ''' reading_frame is relevant only when restarted
    print_this = "\t\t\t\tEstimated reading_frame = total_ps/((dt)*1,000) = " + str(float(total_ps)/((float(dt))*1000.0)) + " ps \n"
    print print_this
    f_out.write(print_this)
    '''
    
    print_this = "\tTherefore, total mdrun running time was: " + str(total_ps) + " pico (10^-12) second" + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    return nsteps, total_ps
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
    
    gro_extraction_note_file = open("gro_extraction.txt","w+")
    
    write_this = "\nExtract 3 highest cc gro (among the whole run, not just the last run)\n"
    gro_extraction_note_file.write(write_this)
    print write_this
    
    # Although I assign number_of_steps_for_cryo_fit*2 as a new number_of_steps_for_cryo_fit,
    #due to state.cpt, mdrun runs only until a new number_of_steps_for_cryo_fit INCLUDING FORMERLY RAN STEPS  
    
    # Previous traj.xtc is erased (not keeping previous record) every time when em_weight or number_of_steps_for_cryo_fit is reassigned.
    # Therefore, cc_record_full_renumbered should NOT be used for extrqcting gro. It should be used only for overall cc change.

    
    #'''
    result = '' # initial temporary assignment
    if (this_is_test == "True"): # test
        result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    else: # default running
        
        # adjust step number if cryo_fit restarted for longer steps
        if (os.path.isfile("../restart_record_for_longer_steps.txt") == True): # this exists only when cryo_fit restarted with longer steps, not with higher map
            adjust_step_number ()
            #os.remove("../restart_record_for_longer_steps.txt") # only for development, keep this file
        
        if (no_rerun == "False"): # default running
            # this cc_record is step_adjusted if restarted
            if (os.path.isfile("cc_record_adjusted_step_use_for_extraction") == False):
                print "cc_record_adjusted_step_use_for_extraction is not found, please email doonam@lanl.gov"
                exit(1)
            result = os.popen("cat cc_record_adjusted_step_use_for_extraction | sort -nk5 -r | head -3").readlines()
        else:
            result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
            
    #'''
        
    #result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    
    write_this = "3 highest cc steps that need to be extracted:" + str(result) + "\n\n"
    gro_extraction_note_file.write(write_this)
    print write_this
    
    if (len(result) == 0):
        print "no steps to be extracted, please email doonam@lanl.gov"
        exit(1)
    
    nsteps, total_ps = get_nsteps_total_ps(gro_extraction_note_file, cryo_fit_path)
    
    for i in range(len(result)):
        splited = result[i].split()
        target_step = splited[1]
        cc = splited[4]
        
        write_this = "\n\nCryo_fit will extract a gro file from this target_step: " + str(target_step)
        gro_extraction_note_file.write(write_this)
        print write_this
        
        #extract_gro(gro_extraction_note_file, cryo_fit_path, nsteps, total_ps, target_step, i, cc)
        
        returned = extract_gro(gro_extraction_note_file, cryo_fit_path, nsteps, total_ps, target_step, i, cc)
        if (returned == "empty"):
            exit(1)
    
    gro_extraction_note_file.close()
