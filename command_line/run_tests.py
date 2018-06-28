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

cryo_fit_repository_dir = libtbx.env.dist_path("cryo_fit")

if (__name__ == "__main__") :
  
  # added by Nigel? temporarily disabled
  assert len(os.listdir(os.getcwd()))==0, 'run in an empty directory'

# Locate phenix executable
  print "This cryo_fit.run executable comes from ", cryo_fit_repository_dir

# copy input files to a current folder (although it may take longer time by copying these files, it is more organized \
# with respect to development in the long term)

  pdb_file_name = 'tRNA_tutorial.pdb'
  map_file_name = 'tRNA_tutorial.map'
  
  # added by Nigel? temporarily disabled since on Doonam's macbook pro can't run this
  '''
  shutil.copyfile(os.path.join(cryo_fit_repository_dir,
                               'tutorial_input_files',
                               pdb_file_name), pdb_file_name)
  shutil.copyfile(os.path.join(cryo_fit_repository_dir,
                               'tutorial_input_files',
                               sit_file_name), sit_file_name)
  '''
  
  pdb_file_name_w_path = os.path.join(cryo_fit_repository_dir,
                               'tutorial_input_files',
                               pdb_file_name)
  
  map_file_name_w_path = os.path.join(cryo_fit_repository_dir,
                               'tutorial_input_files',
                               map_file_name)
  
  #command_string = "phenix.cryo_fit %(pdb_file_name)s %(map_file_name)s" % locals()
  command_string = "phenix.cryo_fit %(pdb_file_name)s %(map_file_name)s devel=True" % locals()
  
# Start the simplest testing
  print "command that will be executed: ", command_string
  print "(for your information) this run_tests took 2 minutes on Doonam's macbook pro"
  print '\n ~> %s\n' % command_string

  
  # as of 06/28/2018, below printout on the screen didn't print on Doonam's mac screen
  # so temporarily disabled
  '''
  rc = libtbx.easy_run.go(command=command_string)
  print '*'*80
  for line in rc.stdout_lines:
    print line

  print '*'*80
  print rc.stderr_lines
  '''
  
  # temporarily use this instead
  libtbx.easy_run.call(command=command_string)