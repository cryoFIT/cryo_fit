import glob, iotbx.pdb, os, shutil, time
from libtbx import easy_run
from subprocess import check_output


def check_whether_the_step_was_successfully_ran(step_name, check_this_file):
  if (os.path.isfile(check_this_file)):
    returned_file_size = file_size(check_this_file)
    if (returned_file_size > 0):
      if (step_name != "Step 8"): # for step_8 (drawing a graph), determining a success now is early
        print step_name, " successfully ran"
      return 1
  print step_name, " didn't successfully ran"
  exit(1)
############## end of check_whether_the_step_was_successfully_ran function


def file_size(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size
############## end of file_size()


def locate_Phenix_executable():
    path = check_output(["which", "phenix.cryo_fit"])
    splited = path.split("/")
    command_path = ''
    for i in range(len(splited)-3):
      command_path = command_path + splited[i] + "/"
    command_path = command_path + "modules/cryo_fit/"
    print "\tUser's phenix.cryo_fit executable comes from ", command_path
    return command_path
########### end of locate_Phenix_executable function


def run(prefix="tst_step_1"):
  """
  Exercise phenix.cryo_fit step_1 with all defaults
  """
  
  assert (os.path.isfile("data/input_for_all/tRNA_EFTU_within_10.ccp4") == True)
  assert (os.path.isfile("data/input_for_all/regression_tRNA_EFTU_within_10.pdb") == True)
  
  cmd = " ".join([
    "phenix.cryo_fit",
    "data/input_for_all/regression_tRNA_EFTU_within_10.pdb",
    "data/input_for_all/tRNA_EFTU_within_10.ccp4",
    "step_2=False",
    "step_3=False",
    "step_4=False",
    "step_5=False",
    "step_6=False",
    "step_7=False",
    "step_8=False"
    ])
  print cmd
  easy_run.call(cmd)
    
  starting_dir = os.getcwd()
  new_path = starting_dir + "/steps/1_make_gro"
  os.chdir( new_path )
  
  for check_this_file in glob.glob("*_by_pdb2gmx.gro"): # there will be only one *_by_pdb2gmx.gro file
    this_step_was_successfully_ran = check_whether_the_step_was_successfully_ran("Step 1", check_this_file)
    if (this_step_was_successfully_ran != 1):
        print "failed, sleep for 10,000 seconds"
        time.sleep(10000) # so that it is recognized instantly
    assert (this_step_was_successfully_ran != 0)
  os.chdir(starting_dir)
  shutil.rmtree("steps")
############# end of run function


if (__name__ == "__main__"):
  t0=time.time()
  run()
  print "Seconds for this tst_step_1: %6.4f"%(time.time()-t0)
  print "OK"
