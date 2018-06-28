import glob, os, subprocess, sys, time

def prepare_minimize(number_of_available_cores, number_of_cores_to_use, ns_type):
    home_bin_cryo_fit_bin_dir = know_home_cryo_fit_bin_dir_by_ls_find()
    
    command_used = first_prepare_minimization(home_bin_cryo_fit_bin_dir, ns_type, number_of_available_cores, \
                                                  number_of_cores_to_use)
    
    f_out = open('log.step_4_1_minimization_used_command', 'wt')
    f_out.write(command_used)
    f_out.write("\n")
    f_out.close()
# end of minimize function

if (__name__ == "__main__"):
    args=sys.argv[1:]
    input_tpr_name = args[0] # input_tpr_name not used in this .py, but specify for former calling
    command_path = args[1]
    
    common_functions_path = command_path + "/common_functions/"
    sys.path.insert(0, common_functions_path)
    from common_functions import *
    
    ns_type = args[2]
    number_of_available_cores = int(args[3])
    number_of_cores_to_use = args[4]
    
    prepare_minimize(number_of_available_cores, number_of_cores_to_use, ns_type)
# end of if (__name__ == "__main__")