import glob, os, subprocess, sys, time

def minimize(input_tpr_name, number_of_available_cores, number_of_cores_to_use, ns_type):
    bool_enable_mpi = know_output_bool_enable_mpi_by_ls()
    print "\tbool_enable_mpi:", bool_enable_mpi
    home_bin_cryo_fit_bin_dir = know_home_cryo_fit_bin_dir_by_ls_find()
    bool_minimization = True
    target_map = '' # no map for minimization
    output_file_format = '' # no map for minimization
    
    f_out = open('log.step_3_2', 'wt')
    bool_just_get_input_command = True
    
    command_that_will_be_used = first_prepare_for_minimization_cryo_fit(bool_minimization, \
                                                                   bool_just_get_input_command, \
                                                                home_bin_cryo_fit_bin_dir, \
                                                                ns_type, number_of_available_cores, \
                                                                number_of_cores_to_use,\
                                                                target_map, output_file_format, "bogus_output_file_name_prefix")
    f_out.write("will_be_used_input_command\n")
    write_this_input_command = str(command_that_will_be_used) + "\n\n"
    f_out.write(write_this_input_command)
    f_out.close()
    
    bool_just_get_input_command = False
    
    command_used = first_prepare_for_minimization_cryo_fit(bool_minimization,\
                                                           bool_just_get_input_command, \
                                                           home_bin_cryo_fit_bin_dir,\
                                                           ns_type, \
                                                          number_of_available_cores, number_of_cores_to_use,\
                                                          target_map, output_file_format, "bogus_output_file_name_prefix")
    
    f_out = open('log.step_3_2_minimization_real_command', 'wt')
    f_out.write(str(command_used))
    f_out.close()
# end of minimize function

if (__name__ == "__main__"):
    args=sys.argv[1:]
    input_tpr_name = args[0]
    command_path = args[1]
    
    common_functions_path = command_path + "/command_line/"
    sys.path.insert(0, common_functions_path)
    from common_functions import *
    
    ns_type = args[2]
    number_of_available_cores = int(args[3])
    number_of_cores_to_use = args[4]
    
    minimize(input_tpr_name, number_of_available_cores, number_of_cores_to_use, ns_type)
# end of if (__name__ == "__main__")