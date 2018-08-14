import glob, os, subprocess, sys, time
from os.path import exists
from subprocess import check_output, Popen, PIPE # for FFTW_INSTALL
from os.path import expanduser # to find home_dir
import platform

# some header(s) among these are needed for libtbx.env.dist_path
#from cctbx import maptbx # commented for Jun Dong at Centos 7
import iotbx.pdb
import iotbx.pdb.mmcif
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out
import mmtbx.model
import mmtbx.utils

cryo_fit_repository_dir = libtbx.env.dist_path("cryo_fit")
print "cryo_fit_repository_dir:", cryo_fit_repository_dir

common_functions_path = cryo_fit_repository_dir + "/common_functions/"
sys.path.insert(0, common_functions_path)
print "common_functions_path:",common_functions_path
from common_functions import  * # (sometimes) ImportError: No module named libtbx at doonam's newest personal macbookpro

def add_path(GMX_MD_INSTALL, shell):
  
  if (shell == "bash"):
    home_dir = expanduser("~")
    add_this = "\n\nexport PATH=\"" + str(GMX_MD_INSTALL) + "/bin\":$PATH # added by cryo_fit installation\n\n"
    
    path_file = os.path.join(home_dir, '.bash_profile')
    if (os.path.isfile(path_file) == True):
      print "~/.bash_profile exists"
      print "\ncryo_fit_installation will add ", add_this, " to ~/.bash_profile"
    else:
      path_file = os.path.join(home_dir, '.bashrc')
      if (os.path.isfile(path_file) == True):
        print "~/.bashrc exists"
        print "\ncryo_fit_installation will add ", add_this, " to ~/.bashrc"
      else:
        print "both ~/.bashrc and ~/.bash_profile do not exist"
        path_file = None
    if (path_file != None):
      f = open(path_file, 'a') # append
      f.write(add_this)
      f.close()
    return path_file
  else:
    print "User may use cshell or zshell"
    print "Cshell needs to edit ~/.csh"
    return GMX_MD_INSTALL
  # adding PATH for pdb2gmx at phenix/build/bin and cryo_fit/bin at the same causes forcefield error or segfault (7/6/2018)
# end of add_path (GMX_MD_INSTALL, shell)

def clean ():
  color_print ("Hit enter key to clean", 'green')
  raw_input()
  
  command_string = "make distclean" # use distclean, not clean #http://www.gromacs.org/Documentation/Installation_Instructions_4.5
  color_print ("command: ", 'green')
  print command_string
  libtbx.easy_run.call(command=command_string)  
  
  color_print ("\nIf you see a message of ", 'green')
  print "    \"make: *** No rule to make target `clean'.  Stop.\""
  color_print ("it is OK.", 'green')
  color_print ("It just means that you never configured/compiled gromacs_cryo_fit at here before.", 'green')
  
  color_print ("\nOr, if you see the \"make disclean\" ended with", 'green')
  print "    rm -rf .libs _libs"
  print "    rm -f *.lo "
  color_print ("it is OK as well.", 'green')
  color_print ("It just means that you just removed former residual partial/full compilation.\n", 'green')
  color_print ("\nMake sure that there was no error", 'green')
# end of clean function

def configure_cryo_fit (GMX_MD_INSTALL, GMX_MD_SRC, enable_mpi, enable_fftw, enter_all):
  home_dir = expanduser("~")
  print "\ncommand: cd ", GMX_MD_SRC
  os.chdir(GMX_MD_SRC)

  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
  
  start_time_configure = time.time()
  if (enable_mpi == "Y" and enable_fftw == "Y"):
    FFTW_INSTALL = get_FFTW_INSTALL_path (home_dir)
    command_string = "./configure --prefix=" + GMX_MD_INSTALL + " --enable-mpi --enable-float LDFLAGS=-L" + FFTW_INSTALL + "/lib" + " CPPFLAGS=-I" + FFTW_INSTALL + "/include"
    color_print ("command: ", 'green')
    print command_string
    color_print ("\nHit enter key to configure.", 'green')
    raw_input()
    libtbx.easy_run.call(command=command_string)
  elif (enable_mpi == "Y" and enable_fftw == "N"):
    command_string = "./configure --prefix=" + GMX_MD_INSTALL + " --enable-mpi --enable-float --with-fft=fftpack"
    color_print ("command: ", 'green')
    print command_string
    color_print ("\nHit enter key to configure.", 'green')
    raw_input()
    libtbx.easy_run.call(command=command_string)
  elif (enable_mpi == "N" and enable_fftw == "N"):
    command_string = "./configure --prefix=" + GMX_MD_INSTALL + " --enable-float --with-fft=fftpack"
    color_print ("command: ", 'green')
    print command_string
    
    if (str(enter_all) != "True"):
      color_print ("\nHit enter key to configure.", 'green')
      raw_input()
    libtbx.easy_run.call(command=command_string)
  
  print '#'*105
  color_print ("\n\nCheck whether it was configured without any error.", 'green')
  color_print ("Was your configuration ended with something like this?", 'green')
  color_print ("\t...config.status: creating src/config.h", 'green')
  color_print ("\t   config.status: executing depfiles commands", 'green')
  if (enable_mpi == "Y"):
    color_print ("   ...", 'green')
    color_print ("   WARNING:", 'green')
    color_print ("      There are known problems with some MPI implementations:", 'green')
    color_print ("      OpenMPI version < 1.4.1", 'green')
    color_print ("      MVAPICH2 version <= 1.4.1", 'green')
  
  if (str(enter_all) != "True"):
    color_print ("Press Y or N and hit enter", 'green')
    configure_result = raw_input()
  else:
    configure_result = "Y"
    
  if (configure_result != "Y" and configure_result != "y"):
    color_print ("I'm sorry to hear that your configuration didn't go well", 'red')
    color_print ("\nWhen Doonam saw this error", 'red')
    color_print ("\t...checking whether the C compiler works... no", 'red')
    color_print ("\t   configure: error: in `/Users/doonam/src/gromacs-4.5.5_cryo_fit_added':", 'red')
    color_print ("\t   configure: error: C compiler cannot create executables", 'red')
    color_print ("\t   See `config.log' for more details", 'red')
    color_print ("\nInstalling commandline tool helped", 'red')
    color_print ("\thttp://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/ worked well for his macOS 10.12.5\n", 'red')
    color_print ("exit now", 'red')
    exit(1)
  else:
    color_print ("OK, I'm glad to hear that your configuration went well.", 'green')
    
  end_time_configure = time.time()
  print "configuration"
  color_print ((show_time(start_time_configure, end_time_configure)), 'green')
  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
# end of configure_cryo_fit function


''' # keep for now
def get_FFTW_INSTALL_path (home_dir):
  color_print ("Hit enter key to locate FFTW install path", 'green')
  raw_input()

  command_script = "ls -d " + home_dir + "/bin/fftw*"
  color_print ("command: ", 'green')
  print command_script
  FFTW_INSTALL = subprocess.check_output(command_script,
                                   stderr=subprocess.STDOUT,
                                   shell=True)
  print "A folder that installed FFTW: ", FFTW_INSTALL
  if ((FFTW_INSTALL) == False):
      print "Your computer seem to not have installed fftw yet."
      print "Please FULLY install fftw first before installing this cryo_fit, because the cryo_fit installation needs to specify fftw installation folder location."
      print "You are welcome to download the fftw at http://www.fftw.org/download.html"
      print "Then, you may install the fftw, for example, python /Users/doonam/bin/phenix-dev-2747/modules/cryo_fit/steps/0_install_cryo_fit/2_runme_to_install_fftw.py fftw-3.3.6-pl2.tar.gz"
      print "Hit enter key to exit."
      raw_input()
      exit(1)

  FFTW_INSTALL = FFTW_INSTALL[:-1]
  print "User's computer seems to have installed fftw already at: ", FFTW_INSTALL
  
  command_string = "export CPPFLAGS=-I" + FFTW_INSTALL + "/include/"
  color_print ("command: ", 'green')
  print command_string
  libtbx.easy_run.call(command=command_string)

  command_string = "export LDFLAGS=-L" + FFTW_INSTALL + "/lib/"
  color_print ("command: ", 'green')
  print command_string
  libtbx.easy_run.call(command=command_string)
  
  return FFTW_INSTALL
# end of get_FFTW_INSTALL_path function
'''

    
def install_gromacs_cryo_fit(zipped_file, *args):
  color_print ("If you need a troubleshooting, either try to run each sentence in this script or contact Doo Nam Kim (doonam@lanl.gov)\n", 'green')
  
  starting_dir = os.getcwd()
  color_print ("\nCurrent start working directory: ", 'green')
  print starting_dir
    
  splited = zipped_file.split("/")
  wo_ext = str(splited[len(splited)-1])
  splited = wo_ext.split(".zip")
  gromacs_cryo_fit_file_name = splited[0]
  
  GMX_MD_INSTALL = os.path.abspath(install_path) # abs path is needed for configure  
  
  '''
  # to avoid network related warning message during GUI based cryo_fit in MacOS,
  # and to avoid some version related mpi failure (like Karissa's case),
  # enable_mpi = "N" by default, enable_mpi = "N" will support most cores
  # (like 32 cores) anyway, and runs with same speed as with enable_mpi=Y
  '''
  enable_mpi = "N" # "mpi will not be enabled during configuration. Threads will be used instead."
  enable_fftw = "N"
  
  #make_this_folder_if_not_exists(GMX_MD_INSTALL)
  #GMX_MD_INSTALL  = os.path.join(GMX_MD_INSTALL, "bin_keep_this_for_executables")
  #make_this_folder_if_not_exists(GMX_MD_INSTALL)
  
  GMX_MD_SRC  = os.path.join(install_path, "source_keep_this_otherwise_segfault") # so that any case can be handled
  make_this_folder_if_not_exists(GMX_MD_SRC)
  
  command_string = "cp " + zipped_file + " " + GMX_MD_SRC
  os.chdir(GMX_MD_SRC)


  ################### unzip ###################
  command_string = "unzip " + zipped_file
  color_print ("\ncommand: ", 'green')
  print command_string

  color_print ("\nIf you see", 'green')
  print "   replace __MACOSX/gromacs_cryo_fit/._.compile2.bat.swp? [y]es, [n]o, [A]ll, [N]one, [r]ename"
  color_print ("Doonam recommends to press A\n", 'green')
  
  print "enter_all:", enter_all
  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
  
  start_time_unzip = time.time()
  libtbx.easy_run.call(command=command_string)
  end_time_unzip = time.time()
  message = "unzipping " + zipped_file
  print message
  color_print ((show_time(start_time_unzip, end_time_unzip)), 'green')

  
  ################### configure ###################
  GMX_MD_SRC = os.path.join(GMX_MD_SRC, gromacs_cryo_fit_file_name) # redefine GMX_MD_SRC for configure
  configure_cryo_fit (GMX_MD_INSTALL, GMX_MD_SRC, enable_mpi, enable_fftw, enter_all)
        

  ################### Make ###################
  core_numbers_to_use = ''
  if (str(enter_all) != "True"):
    core_numbers_to_use = decide_number_of_cores_to_use(1)
  else:
    core_numbers_to_use = 3
  
  make_command_string = ''
  try:
    debug
  except:
    debug = 0
  if (debug == True):
    make_command_string = "make --debug -j " + str(core_numbers_to_use)
  else:
    make_command_string = "make -j " + str(core_numbers_to_use)
  color_print ("make command: ", 'green')
  print make_command_string
  
  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to do \"Make\"", 'green')
    raw_input()
  
  start_time_make = time.time()
  libtbx.easy_run.call(command=make_command_string)
  end_time_make = time.time()
  
  print '#'*105
  color_print ("\nWas your Make done without any error?.", 'green')
  color_print ("The \"Make\" should have been ended with this kind of message.", 'green')
  color_print ("\t...make[2]: Nothing to be done for `all-am'.", 'green')
  color_print ("\t   make[1]: Nothing to be done for `all-am'.\n", 'green')
  
  make_result = ''
  if (str(enter_all) != "True"):
    color_print ("Press Y or N and enter key.", 'green')
    make_result = raw_input()
  else:
    make_result = "Y"
  
  if (make_result != "Y" and make_result != "y"):
    color_print ("I'm sorry to hear that your Make didn't go well", 'red')
    if (enable_mpi == "Y"):
      color_print ("If you see a message of", 'red')
      color_print ("\t\"configure: error: Cannot compile and link MPI code with cc\",", 'red')
      color_print ("consider to install mpicc by openmpi first", 'red')
      color_print ("\nYou can install openmpi by python /Users/ftuser/bin/phenix-dev-2880/modules/cryo_fit/steps/0_install_cryo_fit/1_runme_to_install_openmpi.py ~/Downloads/openmpi-2.1.1.tar.gz", 'red')
      color_print ("\nIf", 'red')
      color_print ("    \"an error ./Xstuff.h:51:10: fatal error: 'X11/Xresource.h' file not found\" ", 'red')
      color_print ("occurred, it needs to be fixed for proper installation.", 'red')
      color_print ("\nWhen doonam's macbook used /Users/doonam/EMAN2/bin/mpicc, this error occurred.", 'red')
      color_print ("Using /usr/local/bin/mpicc solved this error.\n", 'red')
    color_print ("exit now", 'red')
    exit(1)
  else:
    color_print ("OK, I'm glad to hear that your Make went well.", 'green')
    
  print "Make"
  color_print ((show_time (start_time_make, end_time_make)), 'green')
  
  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()

  ################### Installation of all gromacs executables ###################
  start_time_install = time.time()
  make_install_command_string = ''
  
  if (debug == True):
    make_install_command_string = "make install -g" # this command inflates only, doesn't make executables
  else:
    make_install_command_string = "make install"
  color_print ("make_install_command_string: ", 'green')
  print make_install_command_string
  
  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to install cryo_fit.", 'green')
    raw_input()
  libtbx.easy_run.call(command=make_install_command_string)
  
  print '#'*105
  color_print ("\n\nCheck whether the installation was done without any error.\n", 'green')
  color_print ("Was your installation ended with this kind of message?", 'green')
  color_print ("\t...\"GROMACS is installed under ...", 'green')
  color_print ("\t   \"Make sure to update your PATH and MANPATH to find the\"", 'green')
  color_print ("\t...\"If you want links to the executables in /usr/local/bin,", 'green')  
  color_print ("\t   you can issue \"make links\" now.", 'green')   
  color_print ("\t   make[2]: Nothing to be done for `install-data-am'.\"\n\n", 'green')
  
  install_result = '' # just initial value
  
  if (str(enter_all) != "True"):
    color_print ("Press Y or N and enter key.", 'green')
    install_result = raw_input()
  else:
    install_result = "Y"
  
  if (install_result != "Y" and install_result != "y"):
    color_print ("I'm sorry to hear that your Installation didn't go well", 'red')
    color_print ("The installation SHOULD NOT have been ended with this kind of message", 'red')
    color_print ("\t...\"3 errors generated.", 'red')
    color_print ("\t...\"make[3]: *** [copyrite.lo] Error 1", 'red')
    color_print ("\t...\"make[2]: *** [install-recursive] Error 1", 'red')
    color_print ("\t...\"make[1]: *** [install-recursive] Error 1", 'red')
    color_print ("\t...\"make: *** [install-recursive] Error 1\n", 'red')
    color_print ("exit now", 'red')
    exit(1)
  else:
    check_this_file_w_path = GMX_MD_INSTALL + "/bin/mdrun"
    check_whether_install_is_done(check_this_file_w_path)
  
  end_time_install = time.time()
  
  print "\nThe final installation of cryo_fit"
  color_print ((show_time (start_time_install, end_time_install)), 'green')
  
  path_file = add_path(GMX_MD_INSTALL, shell)
  
  if (shell == "bash"):
    print_this = "\nPlease source " + str(path_file) + " or open a new terminal so that " + str(path_file) + " can recognize cryo_fit path (which is " + str(GMX_MD_INSTALL) + ")"
    print print_this
  else:
    print "Please add " , GMX_MD_INSTALL, " into your PATH to run cryo_fit"
    
  if (str(enter_all) != "True"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
    
  '''
  if (enable_mpi == "Y"):
    command_string = "\nexport TMPDIR=/tmp"
    color_print ("\ncommand: ", 'green')
    print command_string
    color_print ("\nThis command was needed to avoid a unexpected error", 'green')
    color_print ("\t\"PMIx has detected a temporary directory name that results in a path that is too long for the Unix domain socket:\"", 'green')
    color_print ("when running cryo_fit in mpi mode", 'green')
    
    color_print ("\nHit enter key to export like this (Edition of ~/.bashrc or ~/.bash_profile is recommended for better convenience).", 'green')
    ######## Doonam needs to code to edit .bashrc automatically
    raw_input()
  '''
##################### end of install_gromacs_cryo_fit ()


def make_this_folder_if_not_exists(GMX_MD_INSTALL):
  if os.path.isdir(GMX_MD_INSTALL):
    print GMX_MD_INSTALL, "already exists"
  else:
    command_string = "mkdir -p " + GMX_MD_INSTALL
    color_print ("\ncommand: ", 'green')
    print command_string, "\n"
    libtbx.easy_run.call(command=command_string)
#################### end of make_this_folder_if_not_exists ()


def id_shell():
  from os import environ
  print "User is using ", environ['SHELL'] , " shell"
  splited = environ['SHELL'].split("/")
  shell = splited[2]
  return shell
#################### end of id_shell ()


if (__name__ == "__main__") :
  total_start_time = time.time()
  shell = id_shell()
  
  args=sys.argv[1:]
  
  if len(args) < 2:
      print "\nPlease specify your downloaded gromacs_cryo_fit zip file and path that you want to install gromacs_cryo_fit"
      print "Usage: python install_cryo_fit.py <gromacs_cryo_fit.zip> <install_path>"
      print "Example usage: python ~/bin/phenix-1.13rc1-2961/modules/cryo_fit/steps/0_install_cryo_fit/install_cryo_fit.py ~/Downloads/gromacs_cryo_fit.zip ~/cryo_fit"
      print "\nWith 2013 macbook pro, this installation took 4 ~ 11 minutes"
      sys.exit("install_cryo_fit.py exits now.")
  
  else:
    zipped_file = args[0] # input cryo_fit zip file
    if zipped_file[len(zipped_file)-4:len(zipped_file)] != ".zip":
      print "please provide .zip file as gromacs_cryo_fit installation file"
      print "exit now"
      exit(1)
    install_path = args[1]
    enter_all = True
    if (len(args) >= 3):
      enter_all = args[2] # if True,e enter Y to all Y/N questions
    debug = False # if True, -g will be added, with a hope that gdb can be ran
    install_gromacs_cryo_fit(zipped_file, install_path, enter_all, shell)
  
  total_end_time = time.time()
  print "\nTotal cryo_fit installation"
  color_print ((show_time(total_start_time, total_end_time)), 'green')
