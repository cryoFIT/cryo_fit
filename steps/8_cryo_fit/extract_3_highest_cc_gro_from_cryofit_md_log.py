# run @steps/8_cryo_fit
import os, subprocess, sys
from os.path import expanduser # to find home_dir
args=sys.argv[1:]
this_is_test = args[0]

def extract_gro(target_step, i):
    for_cryo_fit_mdp_location = ''
    
    if this_is_test == "False":
        for_cryo_fit_mdp_location = "../7_make_tpr_with_disre2/for_cryo_fit.mdp"
    else:
        print "\tthis_is_test in extract_gro:", this_is_test
        for_cryo_fit_mdp_location = "for_cryo_fit.mdp"
        
    grep_dt_string = "grep dt " + for_cryo_fit_mdp_location + " | grep -v when"
    print "\tgrep_dt_string:", grep_dt_string
    result = os.popen(grep_dt_string).read()
    splited = result.split()
    dt = splited[2]
    
    grep_nsteps_string = "grep nsteps " + for_cryo_fit_mdp_location + " | grep -v when"
    result = os.popen(grep_nsteps_string).read()
    splited = result.split()
    nsteps = splited[2]
    
    total_ps = float(dt)*float(nsteps) 
    print "\tTherefore, total mdrun running time was: ", total_ps, "ps (10^-12) second"
    print "\tUser wants to extract a gro file from ", target_step, "steps"
    target_ps = (float(target_step)/float(nsteps))*float(total_ps)
    print "\tTherefore, I will extract a gro file from ", target_ps, "ps"
    
    output_gro_name = "extracted_" + str(target_step) + "_steps_" + str(target_ps) + "_ps.gro"
    os.system("echo 0 > input_parameters") # to select system
    
    home_dir = expanduser("~")
    home_bin_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit/bin"
    cmd = home_bin_cryo_fit_bin_dir + "/trjconv -f traj.xtc -dump " + str(target_ps) + " -o " + str(output_gro_name) + \
          " -s for_cryo_fit.tpr < input_parameters"
    os.system(cmd)
    
    if (i == 0):
        print "\tThis has the highest cc"
        if (target_step == "0"): # works as expected
            print "\tHowever, it was the initial model that a user provided"
            print "\tso don't rename it to cryo_fitted.gro"
            #cmd = "rm cryo_fitted.gro"
            cmd = "mv " + output_gro_name + " user_provided.gro"
            print "\tcmd:", cmd
            os.system(cmd)
        else:
            print "\tso rename it to cryo_fitted.gro"
            cmd = "mv " + output_gro_name + " cryo_fitted.gro"
            print "\tcmd:", cmd
            os.system(cmd)
# end of extract_gro function

if (__name__ == "__main__") :
    cmd = "grep corre md.log > cc_record"
    os.system(cmd)
    
    result = os.popen("cat cc_record | sort -nk5 -r | head -3").readlines()
    #print result
    for i in range(len(result)):
        splited = result[i].split()
        target_step = splited[1]
        print "\ttarget_step: ", target_step
        extract_gro(target_step, i)
