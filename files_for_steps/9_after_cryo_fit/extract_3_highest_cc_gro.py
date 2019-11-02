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
    f = open('../aim_this_step_when_restart.txt', 'r')
    # count # of lines
    number_of_lines = 0
    for line in f:
        number_of_lines = number_of_lines + 1
    f.close()
    
    #last_step_to_be_added = '' # initial
    last_step_to_be_added = 0 # initial
    f = open('../aim_this_step_when_restart.txt', 'r')
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
    
    #os.remove("cc_record") # no longer needed, but keep for development
################# end of def adjust_step_number ()


def extract_gro(gro_extraction_note_file, cryo_fit_path, nsteps, nsteps_from_state_cpt, dt, total_ps, target_step, i, cc):

    # old style when state.cpt is not considered
    #target_ps = (float(target_step)/float(nsteps))*float(total_ps)
    
    # new style when state.cpt is considered
    target_ps = ''
    print_this = ''
    if (nsteps_from_state_cpt != ''):
        print_this = "nsteps_from_state_cpt:",nsteps_from_state_cpt
        print print_this
        gro_extraction_note_file.write(print_this)
        
        target_ps = (float(target_step)/float(nsteps))*float(total_ps) + float(dt)*float(nsteps_from_state_cpt)
        print_this = "\n\ttarget_ps = (float(target_step)/float(nsteps))*float(total_ps) + float(dt)*float(nsteps_from_state_cpt)" + "\n"
        print print_this
        gro_extraction_note_file.write(print_this)
        
        print_this = "target_step:" + str(target_step) + " nsteps:" + str(nsteps) + " total_ps:" + str(total_ps) + "dt: " + str(dt)
        
         
    else:
        target_ps = (float(target_step)/float(nsteps))*float(total_ps)    
        print_this = "\n\ttarget_ps = (float(target_step)/float(nsteps))*float(total_ps)" + "\n"

    print print_this
    gro_extraction_note_file.write(print_this)

    
    print_this = "\tTherefore, the cryo_fit will extract a gro file from " + str(target_ps) + " ps" + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    output_gro_name = "extracted_" + str(target_step) + "_target_step_" + str(target_ps) + "_target_ps.gro"
    #output_gro_name = "extracted_" + str(target_ps) + "_ps.gro"
    
    os.system("echo 0 > input_parameters") # to select system
    
    cmd = cryo_fit_path + "trjconv -f traj.xtc -dump " + str(target_ps) + " -o " + str(output_gro_name) + \
          " -s for_cryo_fit.tpr < input_parameters"
    write_this = "\t" + cmd + "\n"
    print write_this
    gro_extraction_note_file.write(write_this) # "ValueError: I/O operation on closed file"
    os.system(cmd)
    
    returned_file_size = file_size(output_gro_name)
    if (returned_file_size == 0):
        write_this = "Extracted gro file is empty, check step numbers, cryo_fit will exit soon.\n"
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
            if ((users_cc == '') or (users_cc == None)):
                write_this = "User's cc can't be retrieved. Please email doonam.kim@gmail.com"
                print write_this
                gro_extraction_note_file.write(write_this)
                return "empty"
            
            print_this = "\tusers_cc:", users_cc
            print print_this
            gro_extraction_note_file.write(print_this)
            
            if (float(cc) > float(users_cc)):
                print "\tTherefore, rename this gro file to cryo_fitted.gro"
                cmd = "mv " + output_gro_name + " cryo_fitted.gro"
                print "\t", cmd, "\n"
                gro_extraction_note_file.write(cmd)
                os.system(cmd)
                
    os.remove("input_parameters")
    return 1
################# end of extract_gro function


def get_nsteps_total_ps(gro_extraction_note_file, cryo_fit_path):
    md_log_location = "../steps/8_cryo_fit/md.log"
    grep_step_colon = "grep step: " + md_log_location
    print "\tcommand:", grep_step_colon
    gro_extraction_note_file.write(grep_step_colon)
    result = os.popen(grep_step_colon).read()
    state_cpt_used = True
    used_step = ''
    nsteps = ''
    nsteps_from_state_cpt = ''
    try:
        splited = result.split()
        nsteps = splited[1] # actually used step_number
        print_this = "\n\tnsteps when state.cpt was used: " + str(nsteps) + "\n" # this \n at the end is needed for gro_extraction.txt
        gro_extraction_note_file.write(print_this)
        print print_this
        nsteps_from_state_cpt = nsteps
    except:
        state_cpt_used = False
    
    ################ <begin> extract dt
    for_cryo_fit_mdp_location = ''
    if (this_is_test == "False"):
        for_cryo_fit_mdp_location = "../steps/7_make_tpr_with_disre2/for_cryo_fit.mdp"
    else:
        print "\t This is a test in extract_gro"
        for_cryo_fit_mdp_location = "for_cryo_fit.mdp"
    
    grep_dt_string = "grep dt " + for_cryo_fit_mdp_location + " | grep -v when"
    
    #print "\tcommand:", grep_dt_string
    #gro_extraction_note_file.write(grep_dt_string)
    result = os.popen(grep_dt_string).read()
    splited = result.split()
    dt = splited[2]

    print_this = "\ndt:" + dt + "\n" # this \n at the end is needed for gro_extraction.txt
    print print_this
    gro_extraction_note_file.write(print_this)
    ############# <end> extract dt
    
    if (state_cpt_used == False):
        grep_nsteps_string = "grep nsteps " + for_cryo_fit_mdp_location + " | grep -v when"
        
        result = os.popen(grep_nsteps_string).read()
        splited = result.split()
        nsteps = splited[2]
        
        print_this = "\t\nnsteps when state.cpt is not used: " + str(nsteps) + "\n" # this \n at the end is needed for gro_extraction.txt
        gro_extraction_note_file.write(print_this)
        print print_this
        
        gro_extraction_note_file.write(print_this)
    
    total_ps = float(dt)*float(nsteps)
    print_this = "\ttotal_ps = float(dt)*float(nsteps) = " + str(total_ps) + "\n" # this \n at the end is needed for gro_extraction.txt
    print print_this
    gro_extraction_note_file.write(print_this)
    
    print_this = "\tTherefore, total mdrun running time was: " + str(total_ps) + " pico (10^-12) second" + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    print_this = "\tHowever, when extracting gro, cryo_fit may need to consider a fact that whether it was restarted" + "\n"
    print print_this
    gro_extraction_note_file.write(print_this)
    
    return nsteps, nsteps_from_state_cpt, dt, total_ps
################# end of get_nsteps_total_ps ()


def get_users_cc_from_overall_log(log):
  f_in = open(log)
  for line in f_in:
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

    ''' # old style which results in an error
    result = '' # initial temporary assignment
    if (this_is_test == "True"): # test
        result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    else: # default running
        # adjust step number if cryo_fit restarted for longer steps
        if (os.path.isfile("../aim_this_step_when_restart.txt") == True):
            # this exists for only a case when cryo_fit restarted with longer steps, not with higher map
            adjust_step_number()
            #os.remove("../aim_this_step_when_restart.txt") # only for development, keep this file
        
        if (no_rerun == "False"): # default running
            
            if (os.path.isfile("cc_record_adjusted_step_use_for_extraction") == False):
                ##### need to use cc_record when cryo_fit bumped up map_weight only
                result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
                #print "cc_record_adjusted_step_use_for_extraction is not found, please email doonam.kim@gmail.com"
                #exit(1)
            else:
                # number of steps of cc_record is adjusted if it was restarted for longer steps
                result = os.popen("cat cc_record_adjusted_step_use_for_extraction | sort -nk5 -r | head -3").readlines()
        else:
            result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    '''
    
    # new style    
    result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
        
    write_this = "3 highest cc steps that need to be extracted:" + str(result) + "\n\n"
    gro_extraction_note_file.write(write_this)
    print write_this
    
    if (len(result) == 0):
        print "no steps to be extracted, please email doonam.kim@gmail.com"
        exit(1)
    
    nsteps, nsteps_from_state_cpt, dt, total_ps = get_nsteps_total_ps(gro_extraction_note_file, cryo_fit_path)
    
    for i in range(len(result)):
        splited = result[i].split()
        target_step = splited[1]
        cc = splited[4]
        
        write_this = "\n\nCryo_fit will extract a gro file from this target_step: " + str(target_step)
        gro_extraction_note_file.write(write_this)
        print write_this
        
        returned = extract_gro(gro_extraction_note_file, cryo_fit_path, nsteps, nsteps_from_state_cpt, dt, total_ps, target_step, i, cc)
        if (os.path.isfile("extract_gro_failed.txt") == True):
            os.remove("extract_gro_failed.txt")
            
        if (returned == "empty"):
            f= open("extract_gro_failed.txt","w+")
            f.close()
    
    gro_extraction_note_file.close()
