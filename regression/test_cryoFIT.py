# LIBTBX_SET_DISPATCHER_NAME phenix.cryoFIT
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from libtbx import phil
import libtbx.phil.command_line
from subprocess import check_output
import libtbx.load_env

cryo_fit_repository_dir = libtbx.env.dist_path("cryoFIT")

def run ():
    # running this test_cryoFIT is not recommended to be ran at /Users/doonam/bin/phenix-dev-2747/modules/cryo_fit to avoid \
    #git related changes
    
    # copy input files to a current folder (although it may take longer time by copying these files, it is more organized \
    # with respect to development in the long term)

    if (os.path.isfile("minimized_c_term_renamed_by_resnum_oc.gro")):
        print "remove former result file"
        os.remove("minimized_c_term_renamed_by_resnum_oc.gro")
        
    command_string = "cp " + cryo_fit_repository_dir + "/steps/2_clean_gro/*rename*.py ."
    print "\tcommand: ", command_string
    libtbx.easy_run.fully_buffered(command_string)
  
    command_string = "cp " + cryo_fit_repository_dir + "/regression/input_file/minimized.gro ."
    print "\tcommand: ", command_string
    libtbx.easy_run.fully_buffered(command_string)
    
    command_string = "python 1_rename_term_res_to_Cres.py"
    print "\tcommand_string: ", command_string
    libtbx.easy_run.fully_buffered(command_string)
    if (os.path.isfile("minimized_c_term_renamed_by_resnum_oc.gro")):
        print "OK\n"
    else:
        print "regression for cryo_fit didn't run successfully, exit now.\n"
        exit(1)

if __name__=="__main__":
    run()
