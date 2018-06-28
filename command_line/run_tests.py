# LIBTBX_SET_DISPATCHER_NAME phenix.cryo_fit.run_tests
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from subprocess import check_output
import libtbx.load_env
import shutil

cryo_fit_repository_dir = libtbx.env.dist_path("cryoFIT")

if (__name__ == "__main__") :
  assert len(os.listdir(os.getcwd()))==0, 'run in an empty directory'

# Locate phenix executable
  print "This cryo_fit.run executable comes from ", cryo_fit_repository_dir

# copy input files to a current folder (although it may take longer time by copying these files, it is more organized \
# with respect to development in the long term)

  pdb_file_name = 'transmin1_gro.pdb'
  sit_file_name = 'H40-H44_0.5A.sit'
  shutil.copyfile(os.path.join(cryo_fit_repository_dir,
                               'tutorial_input_files',
                               pdb_file_name), pdb_file_name)
  shutil.copyfile(os.path.join(cryo_fit_repository_dir,
                               'tutorial_input_files',
                               sit_file_name), sit_file_name)
  
  command_string = "phenix.cryoFIT %(pdb_file_name)s %(sit_file_name)s" % locals()
  
# Start the simplest testing
  print "command that will be executed: ", command_string
  print "(for your information) this run_tests took 1.5 minutes on Doonam's laptop"
  print '\n ~> %s\n' % command_string
  rc = libtbx.easy_run.go(command=command_string)
  print '*'*80
  for line in rc.stdout_lines:
    print line

  print '*'*80
  print rc.stderr_lines
