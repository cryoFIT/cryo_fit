# LIBTBX_SET_DISPATCHER_NAME cryo_fit.run_tests
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

import glob, iotbx.pdb.hierarchy, os, subprocess, sys, time
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from subprocess import check_output
import libtbx.load_env

cryo_fit_repository_dir = libtbx.env.dist_path("cryo_fit")

if (__name__ == "__main__") :
  # running this run_tests is not recommended to be ran at /Users/doonam/bin/phenix-dev-2747/modules/cryo_fit to avoid \
  #git related changes

# Locate phenix executable
  print "This cryo_fit.run executable comes from ", cryo_fit_repository_dir

# copy input files to a current folder (although it may take longer time by copying these files, it is more organized \
# with respect to development in the long term)

  command_string = "cp " + cryo_fit_repository_dir + "/tutorial_input_files/* ."
  print "\tcommand: ", command_string
  libtbx.easy_run.fully_buffered(command_string)
  
  command_string = "cryo_fit.run transmin1_gro_translated.pdb H40-H44_0.5A.sit"
  
# Start the simplest testing
  # command_string = "cryo_fit.run model=%s map=%s number_of_steps=100" % (
  #   os.path.join(command_path, 'tutorial_input_files', 'transmin1_gro_translated.pdb'),
  #   os.path.join(command_path, 'tutorial_input_files', 'H40-H44_0.5A.sit'),
  #   )
  print "command that will be executed: ", command_string
  print "(for your information) this run_tests took 1.5 minutes on Doonam's laptop"
  print '\n ~> %s\n' % command_string
  rc = libtbx.easy_run.go(command=command_string)
  print '*'*80
  for line in rc.stdout_lines:
    print line

  print '*'*80
  print rc.stderr_lines
