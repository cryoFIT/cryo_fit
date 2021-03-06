# LIBTBX_SET_DISPATCHER_NAME cryo_fit.run_tests
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from subprocess import check_output
import libtbx.load_env
import shutil
import time


cryo_fit_repository_dir = libtbx.env.dist_path("cryo_fit") # # Locate phenix.cryo_fit.run_tests executable

if (__name__ == "__main__") :

    assert len(os.listdir(os.getcwd()))==0, 'run in an empty directory' # added by Nigel so that this test runs in a clear path

    print "This phenix.cryo_fit.run_tests executable comes from ", cryo_fit_repository_dir

    ############# Test 1, tutorial adenylate_kinase, each steps ###############    
    regression_path = os.path.join(cryo_fit_repository_dir,
                                     'regression',
                                     'Adenylate_Kinase')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    for i in range (1,9):
        time_start = time.time()
        command_string = "python tst_step_" + str(i) + ".py"
        print "command_string:", command_string
        rc1 = libtbx.easy_run.call(command=command_string)
        time_end = time.time()
        print "Minutes took for tst_step_" + str(i) + ".py:" + str(( round(((time_end-time_start)/60),2)   ))
        assert rc1==0

    
    time_start = time.time()
    
    command_string = "python tst_step_final.py"
    print "command_string:", command_string
    rc1 = libtbx.easy_run.call(command=command_string)
    
    time_end = time.time()
    print "Minutes took for tst_step_final.py:", ( round(((time_end-time_start)/60),2)   )
    assert rc1==0
    
    
    ############# Test 2, tutorial with GTPase_activation_center, each step ###############
    time_start = time.time()
    regression_path = os.path.join(cryo_fit_repository_dir,
                                     'regression',
                                     'GTPase_activation_center')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    for i in range (1,9):
        time_start = time.time()    
    
        command_string = "python tst_step_" + str(i) + ".py"
        print "command_string:", command_string
        rc2 = libtbx.easy_run.call(command=command_string)
    
        time_end = time.time()
        print "Minutes took for tst_step_" + str(i) + ".py:" + str(( round(((time_end-time_start)/60),2)   ))
        assert rc2==0

    time_start = time.time()
    
    command_string = "python tst_step_final.py"
    print "command_string:", command_string
    rc2 = libtbx.easy_run.call(command=command_string)
    
    time_end = time.time()
    print "Minutes took for tst_step_final.py:", ( round(((time_end-time_start)/60),2)   )
    assert rc2==0
    

    
    ############# Test 3, simple biomolecule all steps without restart ###############
    time_start = time.time()
    regression_path = os.path.join(cryo_fit_repository_dir,
                                     'regression',
                                     'emd_8249')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    
    command_string = "python tst_emd_8249.py" % locals()
    rc3 = libtbx.easy_run.call(command=command_string) # if failed (such as gromacs is not installed), rc = 1
    
    time_end = time.time()
    print "Minutes took for tst_emd_8249.py:", ( round(((time_end-time_start)/60),2)   )
    
    assert rc3==0    
    
    
    ''' # sometimes, exceeds 1 min
    ############# test 4, simple protein all steps allowing restart ###############
    time_start = time.time()
    regression_path = os.path.join(cryo_fit_repository_dir,
                                     'regression',
                                     'emd_8249_restart')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    
    command_string = "python tst_emd_8249_restart.py" % locals()
    rc = libtbx.easy_run.call(command=command_string)
    assert rc==0
    
    time_end = time.time()
    print "Minutes took for tst_emd_8249_restart.py:", ( round(((time_end-time_start)/60),2)   )
    '''
    
    
    
    ############# Test 4, simple protein all steps allowing restart ###############
    time_start = time.time()
    regression_path = os.path.join(cryo_fit_repository_dir,
                                     'regression',
                                     'emd_8249_10_AA_restart')
    print "regression_path:", regression_path
    os.chdir(regression_path)

    time_start = time.time()
    
    command_string = "python tst_emd_8249_10_AA_restart.py" % locals()
    
    rc4 = libtbx.easy_run.call(command=command_string)
    assert rc4==0
    
    time_end = time.time()
    print "Minutes took for tst_emd_8249_10_AA_restart.py:", ( round(((time_end-time_start)/60),2)   )
    
    
    
    
    ############# Test 5, tRNA+EFTU all steps allowing restart ###############
    time_start = time.time()
    regression_path = os.path.join(cryo_fit_repository_dir,
                                     'regression',
                                     'tRNA_EFTU_partial')
    print "regression_path:", regression_path
    os.chdir(regression_path)
    for i in range (1,8):
        time_start = time.time()    
    
        command_string = "phenix.python tst_step_" + str(i) + ".py"
        print "command_string:", command_string
        rc5 = libtbx.easy_run.call(command=command_string)
        assert rc5==0
        
        time_end = time.time()
        print "Minutes took for tst_step_" + str(i) + ".py:" + str(( round(((time_end-time_start)/60),2)   ))
        
    
    

    ##########  don't run this test, tutorial_GTPase_activation_center for all steps, since it takes 2 minutes #####
    ########## tutorial_GTPase_activation_center for each steps run individually anyway ##########
    '''
    pdb_file_name = 'GTPase_activation_center_tutorial.pdb'
    map_file_name = 'GTPase_activation_center_tutorial.map'

    # added by Nigel? temporarily disabled since on Doonam's macbook pro can't run this

    #shutil.copyfile(os.path.join(cryo_fit_repository_dir,
    #                             'tutorial_input_files',
    #                             pdb_file_name), pdb_file_name)
    #shutil.copyfile(os.path.join(cryo_fit_repository_dir,
    #                             'tutorial_input_files',
    #                             sit_file_name), sit_file_name)

    pdb_file_name_w_path = os.path.join(cryo_fit_repository_dir,
                                 'tutorial_input_files',
                                 pdb_file_name)

    map_file_name_w_path = os.path.join(cryo_fit_repository_dir,
                                 'tutorial_input_files',
                                 map_file_name)

    #command_string = "phenix.cryo_fit %(pdb_file_name)s %(map_file_name)s" % locals()
    command_string = "phenix.cryo_fit %(pdb_file_name_w_path)s %(map_file_name_w_path)s devel=True" % locals()

    print "command that will be executed: ", command_string
    print "(for your information) this run_tests took 2 minutes on Doonam's macbook pro"
    print '\n ~> %s\n' % command_string


    # as of 06/28/2018, below printout on the screen didn't print on Doonam's mac screen
    # so temporarily disabled

    # rc = libtbx.easy_run.go(command=command_string)
    # print '*'*80
    # for line in rc.stdout_lines:
    #   print line
    #
    # print '*'*80
    # print rc.stderr_lines

    # temporarily use this instead
    libtbx.easy_run.call(command=command_string)
    '''

