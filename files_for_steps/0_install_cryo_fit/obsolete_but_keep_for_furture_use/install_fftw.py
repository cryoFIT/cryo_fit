# This script was built to be ran outside of phenix directory, therefore it uses os.system instead of libtbx.easy_run.fully_buffered...
import glob, iotbx.pdb.hierarchy, os, platform, subprocess, sys, time
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from subprocess import check_output, Popen, PIPE
from os.path import expanduser # to find home_dir

termcolor_installed = '' # initial value
try:
    from termcolor import colored
    termcolor_installed = True
except Exception:
    termcolor_installed = False
    print "Your computer has no termcolor"
    print "If you want to see cryo_fit install helper's comment in color, download from https://pypi.python.org/pypi/termcolor, then do python setup.py install"
    print "Hit enter key to continue"
    raw_input()

# this is needed to import all common functions
path = check_output(["which", "cryo_fit.run"])
splited = path.split("/")
command_path = ''
for i in range(len(splited)-3):
  command_path = command_path + splited[i] + "/"
command_path = command_path + "modules/cryo_fit/"
common_functions_path = command_path + "command_line/"
sys.path.insert(0, common_functions_path)
from common_functions import *
 
def install_fftw(fftw_tar_gz):

# Check current working directory.
    starting_dir = os.getcwd()
    #print "Current working directory: %s" % starting_dir
    color_print ("Current working directory: ", 'green')
    print starting_dir

    splited = fftw_tar_gz.split("/")
    wo_ext = str(splited[len(splited)-1])
    splited = wo_ext.split(".tar")
    fftw_file_name = splited[0]
    color_print ("fftw installation file name: ", 'green')
    print fftw_file_name

    splited = starting_dir.split("/")
    home_dir = expanduser("~")
    home_bin_fftw_dir = home_dir + "/bin/" + fftw_file_name
    color_print ("os.path.isdir(home_bin_fftw_dir): ", 'green')
    print os.path.isdir(home_bin_fftw_dir)
    if os.path.isdir(home_bin_fftw_dir):
        command_script = "rm -rf " + home_bin_fftw_dir
        color_print ("command: ", 'green')
        print command_script
        os.system(command_script)
    command_script = "mkdir -p " + home_bin_fftw_dir
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)

    home_src_dir = home_dir + "/src"
    home_src_fftw_dir = home_dir + "/src/" + fftw_file_name
    if os.path.isdir(home_src_fftw_dir):
        command_script = "rm -rf " + home_src_fftw_dir
        color_print ("command: ", 'green')
        print command_script
        os.system(command_script)

    if (os.path.isdir(home_src_dir) == False):
        command_script = "mkdir " + home_src_dir
        color_print ("command: ", 'green')
        print command_script
        os.system(command_script)

    command_script = "cp " + fftw_tar_gz + " ~/src"
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)

    color_print ("command:  cd ~/src", 'green')
    os.chdir(home_src_dir)

    command_script = "tar -xvf " + fftw_tar_gz
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)

    command_script = "FFTW_INSTALL=" + home_bin_fftw_dir
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)

    command_script = "FFTW_SRC=" + home_src_fftw_dir
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)

    dir_that_mpicc_was_already_installed_at = '' # just initial value
    try:
        command_script = "which mpicc"
        dir_that_mpicc_was_already_installed_at = subprocess.check_output(command_script, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        dir_that_mpicc_was_already_installed_at = False
    if (dir_that_mpicc_was_already_installed_at == False):
        color_print ("Your computer can't find mpicc.", 'green')
        print "Either add the pathway for your mpicc to your .bashrc/.bash_profile's PATH, or install the mpicc."
        print "Hit enter key to exit this fftw installation."
        raw_input()
    else:
        print "Your installed mpicc lives at ", dir_that_mpicc_was_already_installed_at
        command_script = "export MPICC=" + dir_that_mpicc_was_already_installed_at
        print "command that will be executed: ", command_script
        os.system(command_script)

    color_print ("User's computer's operating system: ", 'green')
    print platform.system()
    if (platform.system() == "Darwin"):
        command_script = "export CFLAGS=-m64"
        color_print ("command: ", 'green')
        print command_script
        os.system(command_script)
    else: # for Ubuntu
        command_script = "export CFLAGS=-fPIC"
        color_print ("command: ", 'green')
        print command_script
        os.system(command_script)

    print "command: cd ", home_src_fftw_dir
    os.chdir(home_src_fftw_dir)

    cwd = os.getcwd()
    print "Current working directory: %s" % cwd

    command_script = "make clean"
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)    
    
    color_print ("\nA message of ", 'green')
    color_print ("       \"make: *** No rule to make target `clean'.  Stop.\"", 'green')
    color_print ("is OK. It just means that you never configured/compiled fftw at here before.", 'green')
    color_print ("\nHit Enter key to configure", 'green')
    raw_input()

    start_time_configure = time.time()
    command_script = "./configure --enable-threads --enable-sse2 --enable-shared --enable-float --enable-mpi  --prefix=" + home_bin_fftw_dir
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)

    color_print ("\nCheck whether it was configured without error", 'green')
    color_print ("Last sentences should have looked like something like these", 'green')
    color_print ("   ...", 'green')
    color_print ("   config.status: creating m4/Makefile", 'green')
    color_print ("   config.status: creating fftw.pc", 'green')
    color_print ("   ...", 'green')
    color_print ("   config.status: executing depfiles commands", 'green')
    color_print ("   config.status: executing libtool commands", 'green')
    end_time_configure = time.time()
    color_print ("\nIf there was no error during configuration,", 'green')
    color_print ((show_time("FFTW configuration", start_time_configure, end_time_configure)), 'green')
    color_print ("\nHit Enter key to continue", 'green')
    raw_input()

    dir_that_gcc_6_was_already_installed_at = '' # just initial value
    try:
        command_script = "which gcc-6"
        dir_that_gcc_6_was_already_installed_at = subprocess.check_output(command_script, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        dir_that_gcc_6_was_already_installed_at = False
    if (dir_that_gcc_6_was_already_installed_at == False):
        print "Your computer can't find gcc-6 (Serdal's macbook uses gcc-6)."
        print "Your computer may use /usr/bin/gcc instead (Doonam's macbook uses gcc and had no problem)"
        cores = decide_nproc(1)
        
        command_script = "make -j " + cores
    else: # reference: https://stackoverflow.com/questions/39548840/issue-compiling-fftw
        print "Your installed gcc-6 is at ", dir_that_gcc_6_was_already_installed_at
        cores = decide_nproc(1)
        splited = dir_that_gcc_6_was_already_installed_at.split("\n")
        command_script = "OMPI_CC=" + str(splited[0]) + " make -j " + cores

    color_print ("command: ", 'green')
    print command_script
    color_print ("\nHit enter key to do \"Make\".", 'green')
    raw_input()
    start_time_make = time.time()
    os.system(command_script)

    color_print ("\nCheck whether \"Make\" was done without any error.", 'green') 
    color_print ("Last lines should look like something like this.", 'green')
    print """
      ...
      libtool: link: gcc -D_THREAD_SAFE -O3 -fomit-frame-pointer -mtune=native -fstrict-aliasing -ffast-math -o .libs/fftwf-wisdom fftwf_wisdom-fftw-wisdom.o ../tests/bench-bench.o ../tests/bench-fftw-bench.o  ../threads/.libs/libfftw3f_threads.dylib /Users/doonam/src/fftw-3.3.6-pl2/.libs/libfftw3f.dylib ../.libs/libfftw3f.dylib ../libbench2/libbench2.a -lm
      Making all in m4
      make[2]: Nothing to be done for `all'.
          """

    color_print ("An error message like below is not ok, so consider to contact Doonam.", 'red')
    print """
     "clangclang: : error: error: unknown argument: '-malign-double'unknown argument: '-malign-double'
      clangclang: : warning: warning: -Wl,-no_compact_unwind: 'linker' input unused
      -Wl,-no_compact_unwind: 'linker' input unused
      clangclang: error: : unsupported argument '-q' to option 'Wa,'
      error: unsupported argument '-q' to option 'Wa,'
      ...
      make[1]: *** [all-recursive] Error 1
      make: *** [all] Error 2"
     """
    end_time_make = time.time()
    color_print ("If there was no error during \"Make\",", 'green')
    color_print ((show_time("Make", start_time_configure, end_time_configure)), 'green')
    color_print ("\nHit Enter key to install", 'green')
    raw_input()

    start_time_install = time.time()
    command_script = "make install"
    color_print ("command: ", 'green')
    print command_script
    os.system(command_script)
    end_time_install = time.time()

    color_print ("\nCheck whether compilation was done without error", 'green')
    color_print ("It should have been ended with this kind of message", 'green')
    color_print ("   ...", 'green')
    color_print ("   make[2]: Nothing to be done for `install-exec-am'.", 'green')
    color_print ("   make[2]: Nothing to be done for `install-data-am'.", 'green')
    color_print ("\nIf there was no error,", 'green')
    color_print ((show_time("installation", start_time_configure, end_time_configure)), 'green')

def check_whether_mpicc_is_installed():
    check_whether_mpicc_was_already_installed = '' # just initial value
    try:
        command_script = "which mpicc"
        check_whether_mpicc_was_already_installed = subprocess.check_output(command_script, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        check_whether_mpicc_was_already_installed = False
    if (check_whether_mpicc_was_already_installed != False):
        splited = check_whether_mpicc_was_already_installed.split("/")
        dir_that_mpicc_is_installed_at = '' # initial value
        for i in range(len(splited)-1):
            dir_that_mpicc_is_installed_at = dir_that_mpicc_is_installed_at + splited[i] + "/"
        color_print ("Your computer had already installed mpicc at ", 'green')
        print dir_that_mpicc_is_installed_at
        return "mpicc_installed"
    return "mpicc_not_installed"

def check_whether_openmpi_is_installed():
    check_whether_openmpi_is_installed = '' # just initial value
    try:
        command_script = "which ompi_info"
        check_whether_openmpi_is_installed = subprocess.check_output(command_script, stderr=subprocess.STDOUT,shell=True)
    except Exception:
        check_whether_openmpi_is_installed = False
    if (check_whether_openmpi_is_installed != False):
        splited = check_whether_openmpi_is_installed.split("/")
        dir_that_openmpi_is_installed_at = '' # initial value
        for i in range(len(splited)-1):
            dir_that_openmpi_is_installed_at = dir_that_openmpi_is_installed_at + splited[i] + "/"
        color_print ("\nUser's computer already had installed openmpi at ", 'green')
        print dir_that_openmpi_is_installed_at
        return "openmpi_installed"
    return "openmpi_not_installed"
# end of check_whether_openmpi_is_installed function

def check_whether_openmpi_mpicc_are_installed():
    returned = check_whether_mpicc_is_installed()
    if (returned != "mpicc_installed"):
        color_print ("mpicc is not installed, install it before you try to install fftw, exit now", 'red')
        exit(1)
    
    # returned = check_whether_openmpi_is_installed()
    # if (returned != "openmpi_installed"):
    #     color_print ("openmpi is not installed, consider to install it before you try to install fftw", 'red')
            
if (__name__ == "__main__") :
    check_whether_openmpi_mpicc_are_installed()
    start_time = time.time()
    color_print ("\nThis will install fftw.", 'green')
    color_print ("It will remove existing ~/bin/fftw... and ~/src/fftw... folders (if any) and install fftw there.", 'green')
    color_print ("In order to install fftw, your computer needs to have mpicc (from openmpi or other means) accessible (by PATHWAY)", 'green')
    color_print ("If you need troubleshooting, either try to run each sentence in this script or contact Doo Nam Kim (doonam.kim@pnnl.gov)", 'green')
    color_print ("Hit Enter key to continue.", 'green')
    raw_input()
    
    args=sys.argv[1:]
    if len(args)<1:
        print "Specify your downloaded fftw tar.gz file"
        print "Example usage: python runme_to_install_fftw.py ~/fftw-3.3.6-pl2.tar.gz"
        sys.exit("runme_to_install_fftw.py exits now (expecting a fftw tar file at next run)")
    else:
        fftw_tar_gz = args[0] # pdb input file
        print "input fftw.tar.gz file: ", fftw_tar_gz
        install_fftw(fftw_tar_gz)
  
    end_time = time.time()
    show_time("In total, installation of fftw", start_time, end_time)
    
