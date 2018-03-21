# This script was built to be ran outside of phenix directory.
# Therefore, it uses os.system instead of libtbx.easy_run.fully_buffered
import glob, os, subprocess, sys, time
from os.path import exists
from subprocess import check_output, Popen, PIPE # for FFTW_INSTALL
from os.path import expanduser # to find home_dir
import platform

# this is needed to import all common functions
path = check_output(["which", "phenix.cryo_fit"])
splited = path.split("/")
command_path = ''
for i in range(len(splited)-3):
  command_path = command_path + splited[i] + "/"

command_path = command_path + "modules/cryo_fit/"
# capital letter does matter, so I'm uniting into cryo_fit
#command_path = command_path + "modules/cryofit/"


common_functions_path = command_path + "common_functions/"
sys.path.insert(0, common_functions_path)
#from common_functions import * # ImportError: No module named libtbx at doonam's newest personal macbookpro
from common_functions_without_libtbx import  *

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

def clean ():
  color_print ("Hit enter key to clean", 'green')
  raw_input()
  
  command_string = "make distclean"
  # use distclean, not clean #http://www.gromacs.org/Documentation/Installation_Instructions_4.5
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

def configure_cryo_fit (home_dir, GMX_MD_INSTALL, GMX_MD_SRC, enable_mpi, enable_fftw, enter_all):
  print "\ncommand: cd ", GMX_MD_SRC
  os.chdir(GMX_MD_SRC)

  if (enter_all != "1"):
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
    
    if (enter_all != "1"):
      color_print ("\nHit enter key to configure.", 'green')
      raw_input()
    #libtbx.easy_run.call(command=command_string)
    os.system(command_string)
  
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
  
  if (enter_all != "1"):
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
  color_print ((show_time("configuration", start_time_configure, end_time_configure)), 'green')
  if (enter_all != "1"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
# end of configure_cryo_fit function

def remake_GMX_MD_INSTALL(GMX_MD_INSTALL):
  if os.path.isdir(GMX_MD_INSTALL):
    command_script = "rm -rf " + GMX_MD_INSTALL
    color_print ("command: ", 'green')
    print command_script, "\n"
    os.system(command_script)
  command_script = "mkdir -p " + GMX_MD_INSTALL
  color_print ("\ncommand: ", 'green')
  print command_script, "\n"
  os.system(command_script)
    
def install_gromacs_cryo_fit(zipped_file, *args):
  color_print ("If you need a troubleshooting, either try to run each sentence in this script or contact Doo Nam Kim (doonam@lanl.gov)\n", 'green')
  starting_dir = os.getcwd()
  color_print ("\nCurrent working directory: ", 'green')
  print starting_dir

#  bypass_unzipping = 0
#  if zipped_file[:-4] != ".zip":
#    bypass_unzipping = 1
    
  splited = zipped_file.split("/")
  wo_ext = str(splited[len(splited)-1])
  splited = wo_ext.split(".zip")
  gromacs_cryo_fit_file_name = splited[0]
  color_print ("\ngromacs_cryo_fit_file_name: ", 'green')
  print gromacs_cryo_fit_file_name

  home_dir = expanduser("~")
  GMX_MD_INSTALL = home_dir + "/bin/" + gromacs_cryo_fit_file_name
    
  # color_print ("\nDo you want to install with mpi enabled? (Hi Nigel, please type N)", 'green')
  # color_print ("Unless users are confident with relatively newer version of \
  #              openmpi and plan to use more than 16 cores, \"N\" is recommended", 'green')
  # color_print ("Type either Y or N and hit enter.", 'green')
  # enable_mpi = raw_input()
  
  # to avoid network related warning message during GUI based cryo_fit in MacOS,
  # and to avoid some version related mpi failure (like Karissa's case),
  # enable_mpi = "N" by default, enable_mpi = "N" will support most cores
  # (like 32 cores) anyway, and runs with same speed as with enable_mpi=Y
  enable_mpi = "N"
  
  ''' #deprecated
  if (enable_mpi == "Y"):
    color_print ("mpi will be enabled during configuration. Threads will not be used.", 'green')
  else:
    color_print ("mpi will not be enabled during configuration. Threads will be used instead.", 'green')
  
  color_print ("\nHit enter key to continue.", 'green')
  raw_input()
  '''
  
  ''' # for development purpose only
  color_print ("\nDo you want to install with FFTW enabled? (Hi Nigel, please type N)", 'green')
  color_print ("Type either Y or N and hit enter.", 'green')
  enable_fftw = raw_input()
  if (enable_fftw == "N"):
    color_print ("FFTW will not be enabled during configuration, fftpack will be used instead.", 'green')
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
  else:
    color_print ("FFTW will be enabled during configuration.", 'green')
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
  '''
  enable_fftw = "N"
  
  if (enable_mpi == "Y"):
    GMX_MD_INSTALL = GMX_MD_INSTALL + "_mpi" 
    if (enable_fftw == "Y"):
      GMX_MD_INSTALL = GMX_MD_INSTALL + "_fftw" 
      remake_GMX_MD_INSTALL(GMX_MD_INSTALL)
    else:
      remake_GMX_MD_INSTALL(GMX_MD_INSTALL)
  else:
    if (enable_fftw == "Y"):
      GMX_MD_INSTALL = GMX_MD_INSTALL + "_fftw" 
      remake_GMX_MD_INSTALL(GMX_MD_INSTALL)
    else:
      remake_GMX_MD_INSTALL(GMX_MD_INSTALL)
  
  src_dir = home_dir + "/src"
  GMX_MD_SRC = home_dir + "/src/" + gromacs_cryo_fit_file_name
  color_print ("gromacs source code path: ", 'green')
  print GMX_MD_SRC, "\n"
  if os.path.isdir(GMX_MD_SRC):
    print GMX_MD_SRC
    color_print ("exists, so remove it\n", 'green')
    command_string = "rm -rf " + GMX_MD_SRC
    color_print ("command: ", 'green')
    print command_string, "\n"
#    libtbx.easy_run.call(command=command_string)
    os.system(command_string)

  print "\n", GMX_MD_INSTALL
  color_print ("was made", 'green')
    
  check_this_zip_file = home_dir + "/src/" + gromacs_cryo_fit_file_name + ".zip"  
  if (exists(check_this_zip_file) == False):
    command_string = "cp " + zipped_file + " ~/src"
    color_print ("command: ", 'green')
    print command_string
    #libtbx.easy_run.call(command=command_string)
    os.system(command_string)

  start_time_unzip = time.time()
  color_print ("\ncommand:  cd ~/src", 'green')
  os.chdir(src_dir)

  command_string = "unzip " + zipped_file
  color_print ("\ncommand: ", 'green')
  print command_string

  color_print ("\nIf you see", 'green')
  print "   replace __MACOSX/gromacs_cryo_fit/._.compile2.bat.swp? [y]es, [n]o, [A]ll, [N]one, [r]ename"
  color_print ("Doonam recommends to press A\n", 'green')
  
  print "enter_all:", enter_all
  if (enter_all != "1"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()
  #libtbx.easy_run.call(command=command_string)
  os.system(command_string)
  
  end_time_unzip = time.time()
  message = "unzipping " + zipped_file
  color_print ((show_time(message, start_time_unzip, end_time_unzip)), 'green')
  
  
  
  configure_cryo_fit (home_dir, GMX_MD_INSTALL, GMX_MD_SRC, enable_mpi, enable_fftw, enter_all)
        
  # Make
  core_numbers_to_use = ''
  if (enter_all != "1"):
    core_numbers_to_use = decide_number_of_cores_to_use(1)
  else:
    core_numbers_to_use = 4
    
  command_script = "make -j " + str(core_numbers_to_use)
  
  color_print ("\ncommand: ", 'green')
  print command_script
  
  if (enter_all != "1"):
    color_print ("\nHit enter key to do \"Make\"", 'green')
    raw_input()
  
  start_time_make = time.time()
  os.system(command_script)
  end_time_make = time.time()
  
  print '#'*105
  color_print ("\nWas your Make done without any error?.", 'green')
  color_print ("The \"Make\" should have been ended with this kind of message.", 'green')
  color_print ("\t...make[2]: Nothing to be done for `all-am'.", 'green')
  color_print ("\t   make[1]: Nothing to be done for `all-am'.\n", 'green')
  
  make_result = ''
  if (enter_all != "1"):
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
    
  color_print ((show_time ("Make", start_time_make, end_time_make)), 'green')
  
  if (enter_all != "1"):
    color_print ("\nHit enter key to continue.", 'green')
    raw_input()

  # Installation of cryo_fit (all gromacs executables)
  start_time_install = time.time()
  command_string = "make install"
  color_print ("command: ", 'green')
  print command_string
  
  if (enter_all != "1"):
    color_print ("\nHit enter key to install cryo_fit.", 'green')
    raw_input()
  #libtbx.easy_run.call(command=command_string)
  os.system(command_string)
  print '#'*105
  color_print ("\n\nCheck whether the installation was done without any error.\n", 'green')
  color_print ("Was your installation ended with this kind of message?", 'green')
  color_print ("\t...\"GROMACS is installed under ...", 'green')
  color_print ("\t   \"Make sure to update your PATH and MANPATH to find the\"", 'green')
  color_print ("\t...\"If you want links to the executables in /usr/local/bin,", 'green')  
  color_print ("\t   you can issue \"make links\" now.", 'green')   
  color_print ("\t   make[2]: Nothing to be done for `install-data-am'.\"\n\n", 'green')
  
  if (enter_all != "1"):
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
    color_print ("OK, I'm glad to hear that your installation went well.", 'green')
    
  end_time_install = time.time()

  color_print ((show_time ("The final installation of cryo_fit", start_time_install, end_time_install)), 'green')
  
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
    
  if (enter_all != "1"):
    color_print ("\nHit enter key to finish", 'green')
    raw_input()
# end of install_gromacs_cryo_fit function

if (__name__ == "__main__") :
  total_start_time = time.time()
  
  args=sys.argv[1:]
  if len(args) < 1:
      print "Please specify your downloaded gromacs_cryo_fit zip file"
      print "Usage: python install_cryo_fit.py <gromacs_cryo_fit.zip> <enter_all>"
      print "Example usage: python install_cryo_fit.py ~/gromacs_cryo_fit.zip 0"
      print "If \"enter_all\" equals 1, then all manual checkpoints will be bypassed to facilitate installation"
      print "With 2013 macbook pro, the installation took 9.6 minutes"
      sys.exit("install_cryo_fit.py exits now.")
  elif len(args) == 1:
      zipped_file = args[0] # input cryo_fit zip file
      enter_all = "1" # enter to all Y/N questions
      color_print ("input gromacs_cryo_fit.zip file: ", 'green')
      print zipped_file
      if zipped_file.find("openmpi") != -1:
        color_print ("\nPlease provide cryo_fit installation file, not openmpi installation file.", 'green')
        exit(1)
      install_gromacs_cryo_fit(zipped_file, enter_all)
  else: # len(args) >= 2:
      zipped_file = args[0] # input cryo_fit zip file
      enter_all = args[1] # enter to all Y/N questions
      if (enter_all != "1"):
        color_print ("Hit enter key to continue.", 'green')
        raw_input()
      color_print ("input gromacs_cryo_fit.zip file: ", 'green')
      print zipped_file
      if zipped_file.find("openmpi") != -1:
        color_print ("\nPlease provide cryo_fit installation file, not openmpi installation file.", 'green')
        exit(1)
      install_gromacs_cryo_fit(zipped_file, enter_all)
  total_end_time = time.time()
  color_print ((show_time("Total cryo_fit installation", total_start_time, total_end_time)), 'green')
