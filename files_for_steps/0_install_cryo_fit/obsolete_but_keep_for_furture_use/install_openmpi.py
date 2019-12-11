import glob, os, platform, subprocess, sys, time
from subprocess import check_output, Popen, PIPE
from os.path import expanduser # to find home_dir

# for some unknown reason "ImportError: No module named libtbx" at both
# sparky and lanl 13inch
from common_functions_without_libtbx import *
  
def check_whether_mpicc_is_installed(bool_after_openmpi_installation_trial):
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
        if (bool_after_openmpi_installation_trial == 0):
            color_print ("Your computer had already installed mpicc at ", 'green')
            print dir_that_mpicc_is_installed_at
            color_print ("Therefore, your computer may not need to install this openmpi", 'green')
        if (bool_after_openmpi_installation_trial == 1):
            color_print ("Your computer just installed mpicc at ", 'green')
            print dir_that_mpicc_is_installed_at
        return "mpicc_installed"
    color_print ("mpicc is not installed at your computer", 'red')
    return "mpicc_not_installed"

def check_whether_openmpi_is_installed(bool_after_openmpi_installation_trial, write_this_bin_path, write_this_lib_path):
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
        if (bool_after_openmpi_installation_trial == 0):
            color_print ("\nUser's computer already had installed openmpi at ", 'green')
            print dir_that_openmpi_is_installed_at
            color_print ("Therefore, there is no need to install this openmpi", 'green')
        ''' # doesn't work for macOS
        if (bool_after_openmpi_installation_trial == 1):
            color_print ("\nUser's computer just installed openmpi at ", 'green')
            print dir_that_openmpi_is_installed_at
        '''
        if (bool_after_openmpi_installation_trial == 1):
            color_print ("\nUser's computer just installed openmpi.", 'green')
            color_print ("\nInstall_openmpi.py just added PATH and LD_LIBRARY like", 'green')
            color_print (write_this_bin_path, 'green')
            color_print ("and", 'green')
            color_print (write_this_lib_path, 'green')
            color_print ("Therefore, please open a new tab in your terminal to execute ~/.bash_profile or ~/.bashrc", 'green')
        return "openmpi_installed"
    color_print ("openmpi is not installed at your computer", 'red')
    return "openmpi_not_installed"
# end of check_whether_openmpi_is_installed function

def clean_openmpi():
    command_string = "make clean"
    color_print ("command: ", 'green')
    print command_string
    os.system(command_string)

    color_print ("\nA message of ", 'green')
    color_print ("\t\"make: *** No rule to make target `clean'.  Stop.\"", 'green')
    color_print ("is OK. It just means that you never configured/compiled openmpi at here before", 'green')
    if (check_at_each_step == 1):
        color_print ("\nHit enter key to continue.", 'green')
        raw_input()
# end of clean_openmpi function

def configure_openmpi(home_bin_openmpi_dir):
    configure_start_time = time.time()
    command_string = "./configure --prefix=" + home_bin_openmpi_dir
    color_print ("command: ", 'green')
    print command_string
    if (check_at_each_step == 1):
        color_print ("(This configuration took 3 minutes in Doonam's linux machine and macbookpro)", 'green')
        color_print ("Hit enter key to configure", 'green')
        raw_input()

    color_print ("Not trying to copy screen message to a file", 'green')
    os.system(command_string)

    #configure_log = libtbx.easy_run.call(command_script).raise_if_errors().stdout_lines
#AttributeError: 'int' object has no attribute 'raise_if_errors'

 #   print "use call wo stdout wo raise"
#    configure_log = libtbx.easy_run.call(command_script)
    # configure_log = 0

 #   print "using call wo raise"
#    configure_log = libtbx.easy_run.call(command_script).stdout_lines
#AttributeError: 'int' object has no attribute 'stdout_lines'

#    configure_log = libtbx.easy_run.fully_buffered(command_script).raise_if_errors().stdout_lines

#    print "using subprocess"
    #configure_log = subprocess.check_output(command_script,
     #                                stderr=subprocess.STDOUT,
      #                               shell=True)

#    color_print ("Hit enter key to show log ", 'green')
 #   raw_input()
  #  print "configure_log: ", configure_log
    color_print ("Check whether it was configured without any error", 'green')
    color_print ("If there's no error, it may have been ended with this kind of message", 'green')
    color_print ("    ...", 'green')
    color_print ("    \"Slurm: yes", 'green')
    color_print ("    ssh/rsh: yes", 'green')
    color_print ("    Torque: no\"", 'green')
    configure_end_time = time.time()
    if (check_at_each_step == 1):
        color_print ("Hit enter key to continue.", 'green')
        raw_input()

    print show_time("openmpi installation configuration", configure_start_time, configure_end_time)
    #color_print ("openmpi installation is configured in ", 'green')
    #print round((configure_end_time-configure_start_time)/60, 2)
    #color_print ("minutes (wallclock).", 'green')
    if (check_at_each_step == 1):
        color_print ("Hit enter key to continue.", 'green')
        raw_input()
# end of configure_openmpi function

def prepare_to_install_openmpi(openmpi_tar_gz):
    if (check_at_each_step == 1):
        color_print ("\nPrepare_to_install_openmpi", 'green')
        color_print ("\nHit enter key to continue.", 'green')
        raw_input()
    starting_dir = os.getcwd()
    color_print ("Current working directory: ", 'green')
    print starting_dir, "\n"
    
    splited = openmpi_tar_gz.split("/")
    wo_ext = str(splited[len(splited)-1])
    splited = wo_ext.split(".tar")
    openmpi_file_name = splited[0]
    color_print ("openmpi_file_name: ", 'green')
    print openmpi_file_name, "\n"
    
    home_dir = expanduser("~")
    # ref: https://stackoverflow.com/questions/4028904/how-to-get-the-home-directory-in-python
    color_print ("home_dir:", 'green')
    print home_dir
    
    home_bin_dir = home_dir + "/bin/"
    print "home_bin_dir:", home_bin_dir, "\n"
    if os.path.isdir(home_bin_dir):
        print home_bin_dir
        color_print("already exists, so no need to make it again", 'green')
    else:
        print "command: mkdir ", home_bin_dir
        os.mkdir(home_bin_dir)

    if (check_at_each_step == 1):
        color_print ("Hit enter key to install", 'green')
        raw_input()
    
    home_bin_openmpi_dir = home_bin_dir + openmpi_file_name
    if os.path.isdir(home_bin_openmpi_dir):
        print home_bin_openmpi_dir
        color_print ("exists already, so remove it and remake it", 'green')
        command_string = "rm -rf " + home_bin_openmpi_dir
        color_print ("command: ", 'green')
        print command_string, "\n"
        os.system(command_string)
    print "command: mkdir ", home_bin_openmpi_dir
    os.mkdir(home_bin_openmpi_dir)

    home_src_dir = home_dir + "/src"
    if os.path.isdir(home_src_dir):
        print home_src_dir, " already exists"
    else:
        print "command: mkdir ", home_src_dir
        os.mkdir(home_src_dir)
    
    home_src_openmpi_dir = home_src_dir + "/" + openmpi_file_name
    color_print ("home_src_openmpi_dir:", 'green')
    print home_src_openmpi_dir
    if (check_at_each_step == 1):
        color_print ("Hit enter key to install", 'green')
        raw_input()
    
    if os.path.isdir(home_src_openmpi_dir):
        print home_src_openmpi_dir
        color_print ("exists already, so remove it and remake it", 'green')
        command_string = "rm -rf " + home_src_openmpi_dir
        color_print ("command: ", 'green')
        print command_string
        os.system(command_string)
    
    command_string = "cp " + openmpi_tar_gz + " ~/src"
    color_print ("command: ", 'green')
    print command_string
    os.system(command_string)

    color_print ("command:  cd ~/src", 'green')
    os.chdir(home_src_dir)
    
    command_string = "tar -xvf " + openmpi_tar_gz
    color_print ("command: ", 'green')
    print command_string
    os.system(command_string)

    print "home_src_openmpi_dir: ", home_src_openmpi_dir
    color_print ("command: ", 'green')
    print "cd ",home_src_openmpi_dir
    os.chdir(home_src_openmpi_dir)

    cwd = os.getcwd()
    color_print ("Current working directory: ", 'green')
    print cwd
    return home_bin_openmpi_dir
# end of prepare_to_install_openmpi function

def install_openmpi(openmpi_tar_gz, check_at_each_step):
    color_print("check_at_each_step:", 'green')
    print check_at_each_step
    home_bin_openmpi_dir = prepare_to_install_openmpi(openmpi_tar_gz)
    clean_openmpi()
    configure_openmpi(home_bin_openmpi_dir)

    start_time_install = time.time()
    cores = decide_nproc(check_at_each_step)
    command_string = "make -j " + cores + " all install"
    color_print ("command: ", 'green')
    print command_string
    if (check_at_each_step == 1):
        color_print ("Hit enter key to install", 'green')
        raw_input()

    os.system(command_string)

    color_print ("\nCheck whether compilation was done without error", 'green')
    
    color_print ("\nAn error message\n", 'red')
    color_print ("\tlibtool:   error: 'asm.lo' is not a valid libtool object\n", 'red')
    color_print ("\tmake[2]: *** [libasm.la] Error 1\n", 'red')
    color_print ("\tmake[1]: *** [all-recursive] Error 1\n", 'red')
    color_print ("\tmake[1]: *** [all-recursive] Error 1\n", 'red')
    color_print ("seems to be OK (sometimes openmpi running failed but seems not related to this error message).\n", 'red')
    
    color_print ("If there's no error, it may have been ended with this kind of messsage", 'green')
    color_print ("    ...", 'green')
    color_print ("    /Library/Developer/CommandLineTools/usr/bin/make  install-exec-hook", 'green')
    color_print ("    make[2]: Nothing to be done for `install-data-am'.", 'green')
    color_print ("    ...", 'green')
  
    color_print ("While, an error which is related to Fortran compiler (like below) may/may not be OK. Let's finish this installation trial first.", 'red')
    color_print ("    *** Fortran compiler", 'red')
    color_print ("    checking for gfortran... gfortran", 'red')
    color_print ("    ...", 'red')
    color_print ("    checking if Fortran compiler works... no", 'red')
    color_print ("    It appears that your Fortran compiler is unable to produce working executables.", 'red')
    color_print ("    ...", 'red')
    color_print ("    configure: error: Could not run a simple Fortran program.  Aborting.", 'red')
    end_time_install = time.time()

    color_print ("Trial of openmpi installation took ", 'green')
    if (round((end_time_install-start_time_install)/60/60, 1) < 1):
        print round((end_time_install-start_time_install)/60, 2)
        color_print ("minutes with ", 'green')
    else:
        print round((end_time_install-start_time_install)/60/60, 2)
        color_print ("hours with ", 'green')
    print cores
    color_print ("cores\n", 'green')

    if (check_at_each_step == 1):
        color_print ("Hit enter key to continue.", 'green')
        raw_input()
# end of install_openmpi function

def add_openmpi_pathway(openmpi_file_name):
    # Nigel may change some phenix code, so that cryo_fit and its related openmpi and fftw do no depend on PATH and LD_LIBRARY_PATH
    starting_dir = os.getcwd()
    color_print ("Current working directory: ", 'green')
    print " %s" % starting_dir
    splited = starting_dir.split("/")
    home_dir = expanduser("~")
    home_bin_openmpi_dir = home_dir + "/bin/" + openmpi_file_name
    home_bin_openmpi_bin_dir = home_bin_openmpi_dir + "/bin"
    home_bin_openmpi_lib_dir = home_bin_openmpi_dir + "/lib"
    write_this_bin_path = "export PATH=\"" + home_bin_openmpi_bin_dir + "\"" + ":$PATH\n"
    color_print ("write_this_bin_path: ", 'green')
    print write_this_bin_path
    write_this_lib_path = "export LD_LIBRARY_PATH=\"" + home_bin_openmpi_lib_dir + "\"" + ":$LD_LIBRARY_PATH\n"
    color_print ("write_this_lib_path: ", 'green')
    print write_this_lib_path
    command_script = '' # initial value
    bash_file_with_path = '' # initial value
    bashrc_with_path = home_dir + "/.bashrc" # to know whether user's computer has .bashrc or not
    if (os.path.isfile(bashrc_with_path) == False):
        print "There is no ~/.bashrc"
        print "openmpi bin pathway will be added to ~/.bash_profile's PATH"
        bash_file_with_path = home_dir + "/.bash_profile"
    else:
        print "openmpi bin pathway will be added to ~/.bashrc's PATH"
        bash_file_with_path = home_dir + "/.bashrc"
    current_bash = open( bash_file_with_path , "a+" )
    current_bash.write("\n#Added during openmpi installation.\n")
    current_bash.write(write_this_bin_path)
    current_bash.write("\n#Added during openmpi installation\n")
    current_bash.write(write_this_lib_path)
    current_bash.close()
    command_string = "source " + bash_file_with_path
    # "bash .bashrc" didn't add $PATH in both linux and mac, "source .bashrc" did
    color_print ("command: ", 'green')
    print command_string
    os.system(command_string)
    #libtbx.easy_run.call(command_string)
    
    color_print("If user's computer uses linux, he/she may have observed ", 'green')
    color_print("\"/home/doonam/.bashrc: line 8: return: can only `return' from a function or sourced script\"", 'green')
    color_print("But it is OK to ignore that", 'green')
    
    return write_this_bin_path, write_this_lib_path
# end of add_openmpi_pathway function

def check_whether_openmpi_mpicc_are_installed(check_at_each_step):
    bool_after_openmpi_installation_trial = 0
    returned = check_whether_openmpi_is_installed(bool_after_openmpi_installation_trial, "bogus", "bogus")
    print "check_at_each_step:", check_at_each_step
    if (check_at_each_step == 1):
        if (returned == "openmpi_installed"):
            color_print ("Do you still want to install openmpi again although you already have openmpi? (Y/N)", 'green')
            Y_N = raw_input()
            #print Y_N
            if (Y_N == "N" or Y_N == "n"):
                color_print ("OK, let's not install openmpi", 'green')
                exit(1)

    returned = check_whether_mpicc_is_installed(bool_after_openmpi_installation_trial)
    if (check_at_each_step == 1):
        if (returned == "mpicc_installed"):
            color_print ("Do you still want to install openmpi although you already have mpicc? (Y/N)", 'green')
            Y_N = raw_input()
            #print Y_N
            if (Y_N == "N" or Y_N == "n"):
                color_print ("OK, let's not install openmpi", 'green')
                exit(1)

if (__name__ == "__main__") :
    
    start_time = time.time()
    color_print ("This will install openmpi by python script", 'green')
    color_print ("It will remove existing ~/bin/openmpi... and ~/src/openmpi... folders (if any) and install openmpi there.", 'green')
    color_print ("If you need troubleshooting, either try to run each sentence in this script or contact Doo Nam Kim (doonam.kim@pnnl.gov)", 'green')
    
    args=sys.argv[1:]
    print "len(args):", len(args)
    openmpi_tar_gz = '' # initial value
    check_at_each_step = '' # temporary
    if len(args)<1:
        color_print ("Specify your downloaded openmpi tar file", 'green')
        color_print ("How to use: python runme_to_install_openmpi.py <openmpi.tar file> <check_at_each_step>", 'green')
        color_print ("Example usage 1: python runme_to_install_openmpi.py ~/openmpi-2.1.1.tar.gz", 'green')
        color_print ("Example usage 2: python runme_to_install_openmpi.py ~/openmpi-2.1.1.tar.gz 1", 'green')
        color_print ("Exits now", 'green')
        exit(1)
    elif len(args)==1: 
        openmpi_tar_gz = args[0] # pdb input file
        color_print ("input openmpi.tar.gz file: ", 'green')
        print openmpi_tar_gz
        check_at_each_step = 1
    else:
        openmpi_tar_gz = args[0] # pdb input file
        check_at_each_step = args[1]
        color_print ("input openmpi.tar.gz file: ", 'green')
        print openmpi_tar_gz
    
    if (platform.system() == "Darwin"):
      color_print ("For macOS 10.12, openmpi-2.1.1 works fine.", 'green')
      color_print ("However, for some macOS 10.11.6, openmpi-2.1.1 didn't work like below.", 'red')
      color_print ("\t\"dyld: lazy symbol binding failed: Symbol not found: _clock_gettime\"", 'red')
      color_print ("\t\"Referenced from: /Users/ftuser/bin/openmpi-2.1.1/lib/libopen-pal.20.dylib \
                   (which was built for Mac OS X 10.12)\"", 'red')
      color_print ("\t\"Expected in: /usr/lib/libSystem.B.dylib\"", 'red')
      color_print ("For macOS 10.11.6, openmpi-2.1.0 seems to not work as well.", 'red')
      
      color_print ("However, if you didn't see below messages at macOS 10.11.6, \
                   openmpi-2.1.1 worked well. \
                   (Doonam checked with lanl 13 inch macbook)", 'green')
      color_print ("\t\"dyld: lazy symbol binding failed: Symbol not found: _clock_gettime\"", 'red')
      color_print ("\t\"Referenced from: /Users/ftuser/bin/openmpi-2.1.1/lib/libopen-pal.20.dylib \
                   (which was built for Mac OS X 10.12)\"", 'red')
      color_print ("\t\"Expected in: /usr/lib/libSystem.B.dylib\"", 'red')
      
      color_print ("Hit enter key to continue.", 'green')
      raw_input()
        
    check_whether_openmpi_mpicc_are_installed(check_at_each_step)
    install_openmpi(openmpi_tar_gz, check_at_each_step)
    end_time = time.time()

    splited = openmpi_tar_gz.split("/")
    wo_ext = str(splited[len(splited)-1])
    splited = wo_ext.split(".tar")
    openmpi_file_name = splited[0]

    write_this_bin_path, write_this_lib_path = add_openmpi_pathway(openmpi_file_name)

    bool_after_openmpi_installation_trial = 1
    returned = check_whether_openmpi_is_installed(bool_after_openmpi_installation_trial, write_this_bin_path, write_this_lib_path)
    if (returned != "openmpi_installed"):
        exit(1)
    if (returned == "openmpi_installed"):
        print show_time("Installation of openmpi", start_time, end_time)
        
