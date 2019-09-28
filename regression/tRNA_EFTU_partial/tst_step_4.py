import iotbx.pdb
from libtbx import easy_run
import glob, os, shutil, time

def file_size(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size
# end of file_size()

def check_whether_the_step_was_successfully_ran(step_name, check_this_file):
  if (os.path.isfile(check_this_file)):
    returned_file_size = file_size(check_this_file)
    if (returned_file_size > 0):
      if (step_name != "Step 9"): # for step_9 (drawing a graph), determining a success now is early
        print step_name, " successfully ran"
      return 1
  print step_name, " didn't successfully ran"
  exit(1)
# end of check_whether_the_step_was_successfully_ran function

def run(prefix="tst_step_4"):
  """
  Exercise phenix.cryo_fit step_4 with all defaults
  """
  
  assert (os.path.isfile("data/input_for_all/tRNA_EFTU_within_10.ccp4") == True)
  assert (os.path.isfile("data/input_for_all/regression_tRNA_EFTU_within_10.pdb") == True)
  
  cmd = " ".join([
    "phenix.cryo_fit",
    "data/input_for_all/regression_tRNA_EFTU_within_10.pdb",
    "data/input_for_all/tRNA_EFTU_within_10.ccp4",
    "step_1=False",
    "step_2=False",
    "step_3=False",
    "step_5=False",
    "step_6=False",
    "step_7=False",
    "step_8=False"
    ])
  print cmd
  easy_run.call(cmd)
  
  starting_dir = os.getcwd()
  new_path = starting_dir + "/steps/4_minimize"
  os.chdir( new_path )
  
  the_step_was_successfully_ran = 0
  for check_this_file in glob.glob("*.gro"): # there will be only "will_be_minimized_cleaned.gro"
    the_step_was_successfully_ran = check_whether_the_step_was_successfully_ran("Step 4", check_this_file)  
  
  if (the_step_was_successfully_ran != 1):
        print "failed, sleep for 10,000 seconds"
        time.sleep(10000) # so that it is recognized instantly
  assert (the_step_was_successfully_ran != 0)
  
  os.chdir(starting_dir)
  shutil.rmtree("steps")
############# end of run function


if (__name__ == "__main__"):
  t0=time.time()
  run()
  print "Time: %6.4f"%(time.time()-t0) #on Doonam's lanl 13inch laptop, it took 5 seconds
  print "OK"
