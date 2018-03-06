# run @steps/7_cryo_fit
import os, subprocess, sys
from os.path import expanduser # to find home_dir

def extract_gro(target_step):
    result = os.popen("grep dt ../6_make_tpr_with_disre2/for_cryo_fit.mdp | grep -v when").read()
    splited = result.split()
    dt = splited[2]
    print "\tdt: ", dt
    
    result = os.popen("grep nsteps ../6_make_tpr_with_disre2/for_cryo_fit.mdp | grep -v when").read()
    splited = result.split()
    nsteps = splited[2]
    print "\tnsteps: ", nsteps
    
    total_ps = float(dt)*float(nsteps) 
    print "\tTherefore, total mdrun running time was: ", total_ps, "ps (10^-12) second"
    print "\tUser wants to extract a gro file from ", target_step, "steps"
    target_ps = (float(target_step)/float(nsteps))*float(total_ps)
    print "\tTherefore, we will extract a gro file from ", target_ps, "ps"
    
    output_gro_name = "extracted_" + str(target_step) + "_steps_" + str(target_ps) + "_ps.gro"
    os.system("echo 0 > input_parameters") # to select system
    
    home_dir = expanduser("~")
    home_bin_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit/bin"
    cmd = home_bin_cryo_fit_bin_dir + "/trjconv -f traj.xtc -dump " + str(target_ps) + " -o " + str(output_gro_name) + \
          " -s for_cryo_fit.tpr < input_parameters"
    os.system(cmd)
    
    # output_pdb_name = "extracted_" + str(target_step) + "_steps_" + str(target_ps) + "_ps.pdb"
    # cmd = "editconf -f " + str(output_gro_name) + " -o " + str(output_pdb_name)
    # os.system(cmd)
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
        extract_gro(target_step)
