import glob, os, platform, subprocess
from subprocess import check_output, Popen, PIPE
from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from libtbx.utils import multi_out
import mmtbx.utils
from os.path import expanduser # to find home_dir

termcolor_installed = '' # just initial value
try:
    from termcolor import colored
    termcolor_installed = True
    #print "User's computer has termcolor, installed"
except Exception:
    termcolor_installed = False
    ''' # disable this for now, so that phenix launch will not show this message
    print "\n\tUser's computer has no termcolor"
    print "\tIf you want to see cryo_fit installation helper's comments in color..."
    print "\t1. Download termcolor-1.1.0.tar.gz from https://pypi.python.org/pypi/termcolor"
    print "\t2. Extract termcolor-1.1.0.tar.gz (for example, tar -xvf termcolor-1.1.0.tar.gz)"
    print "\t3. Run \"python setup.py install\" at the extracted folder"
    print "Press any key to continue"
    '''
    #raw_input() # disable this for now, so that Phenix GUI will work


def cif_as_pdb(file_name):  
    try:
      assert os.path.exists(file_name)
      print "Converting %s to PDB format." %file_name
      cif_input = iotbx.pdb.mmcif.cif_input(file_name=file_name)
      hierarchy = cif_input.construct_hierarchy()
      basename = os.path.splitext(os.path.basename(file_name))[0]
      iotbx.pdb.write_whole_pdb_file(
          file_name=basename+".pdb",
          output_file=None,
          processed_pdb_file=None,
          pdb_hierarchy=hierarchy,
          crystal_symmetry=cif_input.crystal_symmetry(),
          ss_annotation=cif_input.extract_secondary_structure(),
          append_end=True,
          atoms_reset_serial_first_value=None,
          link_records=None)
    except Exception, e:
      print "Error converting %s to PDB format:" %file_name
      print " ", str(e)
# end of cif_as_pdb()
    
def color_print(text, color):
    if (termcolor_installed == True):
        print colored (text, color)
    else:
        print text

def decide_number_of_cores_to_use(check_at_each_step):
    number_of_total_cores = know_total_number_of_cores()
    color_print ("User's computer has ", 'green')
    print number_of_total_cores
    color_print ("number of cores in total", 'green')
    print "\n"
    cores = 0 # temporary value
    if check_at_each_step == 1:
        color_print ("Enter how many cores you want to use:", 'green')
        cores = raw_input()
    else:
        if number_of_total_cores > 38:
            cores = 35
        else:
            cores = 2
    return cores
# end of decide_number_of_cores_to_use function


def final_prepare_for_minimization_cryo_fit(bool_just_get_input_command, bool_minimization, \
                                         ns_type, number_of_available_cores, \
                                         number_of_cores_to_use, common_command_string):
    command_used = '' #just initial value
    if (number_of_cores_to_use == "max"):
        if (number_of_available_cores < 4):
            command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
                                                     bool_minimization, 2, \
                                                     ns_type, common_command_string)
        elif (number_of_available_cores < 8):
            command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
                                                     bool_minimization, 4, \
                                                     ns_type, common_command_string)
        elif (number_of_available_cores < 12):
            command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
                                                     bool_minimization, 8, \
                                                     ns_type, common_command_string)
        elif (number_of_available_cores < 16):
            command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
                                                     bool_minimization, 12, \
                                                     ns_type, common_command_string)
        else:
            command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
                                                     bool_minimization, 16, \
                                                     ns_type, common_command_string)
        # based on ribosome benchmark result, maximum number of cores to use is 16
        # elif (number_of_available_cores <= 32):
        #     command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
        #                                              bool_minimization, 24, \
        #                                              ns_type, common_command_string)
        # else:
        #     command_used = minimize_or_cryo_fit(bool_just_get_input_command, \
        #                                              bool_minimization, 32, \
        #                                              ns_type, common_command_string)
    else:
        command_used = minimize_or_cryo_fit(bool_just_get_input_command, bool_minimization, \
                                                 int(number_of_cores_to_use), ns_type, common_command_string)
    return command_used
# end of final_prepare_for_minimization_cryo_fit function


def first_prepare_for_minimization_cryo_fit(bool_minimization, bool_just_get_input_command, \
                                            home_bin_cryo_fit_bin_dir, ns_type, \
                                            number_of_available_cores, number_of_cores_to_use, target_map, \
                                            output_file_format, output_file_name_prefix):
    common_command_string = '' # initial value
    output_file_name = '' # initial value
    if (bool_minimization == True):
        common_command_string = home_bin_cryo_fit_bin_dir + "/mdrun -v -s to_minimize.tpr -c minimized.gro "
    else:
        #print "\toutput_file_name_prefix:", output_file_name_prefix
        if (output_file_name_prefix == "None"):
            output_file_name = "cryo_fitted." + output_file_format
        else:
            output_file_name = output_file_name_prefix + "_cryo_fitted." + output_file_format
        #print "\toutput_file_name: ", output_file_name
        common_command_string = home_bin_cryo_fit_bin_dir + "/mdrun -v -s for_cryo_fit.tpr -mmff -emf " + \
                                target_map + " -c " + output_file_name + " -nosum  -noddcheck "
        
        # -c       : confout.gro  Output       Structure file: gro g96 pdb etc
        # mmff     : Merck Molecular ForceField
        # noddcheck: When inter charge-group bonded interactions are beyond the bonded cut-off distance, \
        #            mdrun terminates with an error message. For pair interactions and tabulated bonds \
        #            that do not generate exclusions, this check can be turned off with the option -noddcheck.
        #-rdd      : real   0  The maximum distance for bonded interactions with DD (nm), \
        #           0 is determined from initial coordinates.
        #           Option -rdd can be used to set the required maximum distance for inter charge-group bonded interactions. \
        #           Communication for two-body bonded interactions below the non-bonded cut-off distance always comes for \
        #           free with the non-bonded communication. Atoms beyond the non-bonded cut-off are only communicated \
        #           when they have missing bonded interactions; this means that the extra cost is minor and nearly independent \
        #           of the value of -rdd. With dynamic load balancing option -rdd also sets the lower limit \
        #           for the domain decomposition cell sizes. By default -rdd is determined by mdrun based on the initial coordinates. \
        #           The chosen value will be a balance between interaction range and communication cost.
    command_used = final_prepare_for_minimization_cryo_fit(bool_just_get_input_command, \
                                                        bool_minimization, \
                                                     ns_type, number_of_available_cores, \
                                                     number_of_cores_to_use, common_command_string)
    return command_used, output_file_name
# end of first_prepare_for_minimization_cryo_fit function

def file_size(fname):
    statinfo = os.stat(fname)
    return statinfo.st_size

def get_fc(complete_set, xray_structure):
  f_calc = complete_set.structure_factors_from_scatterers(
    xray_structure=xray_structure).f_calc()
  return f_calc

def get_fft_map(map_coeffs=None):
    from cctbx import maptbx
    from cctbx.maptbx import crystal_gridding
    ccs=map_coeffs.crystal_symmetry()
    fft_map = map_coeffs.fft_map( resolution_factor = 0.25,
       symmetry_flags=maptbx.use_space_group_symmetry)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded().as_double()
# end of get_fft_map function

# not used for now, but will be used in future
def get_structure_factor_from_pdb_string () :
  prefix = "tmp_iotbx_map_tools"
  pdb_file = prefix + ".pdb"
  mtz_file = prefix + ".mtz"
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string="""\
ATOM      1  N   GLY P  -1     -22.866  -2.627  15.217  1.00  0.00           N
ATOM      2  CA  GLY P  -1     -22.714  -3.068  16.621  1.00  0.00           C
ATOM      3  C   GLY P  -1     -21.276  -3.457  16.936  1.00  0.00           C
ATOM      4  O   GLY P  -1     -20.538  -3.887  16.047  1.00  0.00           O
ATOM      5  H1  GLY P  -1     -22.583  -3.364  14.590  1.00  0.00           H
ATOM      6  H2  GLY P  -1     -22.293  -1.817  15.040  1.00  0.00           H
ATOM      7  H3  GLY P  -1     -23.828  -2.392  15.027  1.00  0.00           H
""")
  xrs = pdb_in.input.xray_structure_simple()
# x-ray structure

#  open(pdb_file, "w").write(pdb_in.hierarchy.as_pdb_string(xrs))
  fc = xrs.structure_factors(d_min=1.5).f_calc()
  #print dir(fc).statistical_mean
# end of get_structure_factor_from_pdb_string function


def kill_mdrun_mpirun_in_linux():
    color_print ("\tkill any existing mdrun jobs (gromacs)", 'green')
    command_string = "top -b -d 1 | head -200 > top_200"
    libtbx.easy_run.call(command=command_string) 
    
    f = open('top_200', 'r')
    for line in f:
      splited = line.split()
      if len(splited) == 12:
        if splited[11] == "mdrun" or splited[11] == "mpirun":
          command_string = "kill " + splited[0]
          print command_string
          libtbx.easy_run.call(command=command_string) 
    f.close()
# end of kill_mdrun_mpirun_in_linux function

def know_number_of_atoms_in_input_pdb(starting_pdb):
    command_string = "cat " + starting_pdb + " | grep ATOM | wc -l"
    print "\tcommand: ", command_string
    num_ATOMs = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    number_of_atoms_in_input_pdb = int(num_ATOMs[0])
    print "\n\tUser's input pdb file, ", starting_pdb, ", has ", number_of_atoms_in_input_pdb, " atoms"
    return number_of_atoms_in_input_pdb
  

def know_output_bool_enable_mpi_by_ls():
    # used exit early for users who didn't install cryofit yet as well
    output_bool_enable_mpi = ''
    home_dir = expanduser("~")
    home_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit"
    #print "\thome_cryo_fit_bin_dir:", home_cryo_fit_bin_dir
    if (os.path.exists(home_cryo_fit_bin_dir) == False):
        print "\nInstall cryo_fit first. Refer http://www.phenix-online.org/documentation/reference/cryo_fit.html"
        print "exit now"
        exit(1)
    output_bool_enable_mpi = False
    return output_bool_enable_mpi
    ''' # now we use non-mpi version only, below is not necessary, not working for CentOS machine as well.
    command_string = "ls ~/bin | grep gromacs-4.5.5_cryo_fit"
    #print "\n\tcommand: ", command_string
    folder_of_cryo_fit = ''
    try: # this try-except seems to be needed for CentOS machine
        folder_of_cryo_fit = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    except:
        print "\nInstall cryo_fit first. Refer http://www.phenix-online.org/documentation/reference/cryo_fit.html"
        print "exit now"
        exit(1)
    #print "folder_of_cryo_fit[0]:", folder_of_cryo_fit[0]
  
    if folder_of_cryo_fit[0] == "gromacs-4.5.5_cryo_fit_added":
        print "\tUser's cryo_fit was installed with enable_mpi=False, so the current cryo_fit will run as enable_mpi = False"
        output_bool_enable_mpi = False  
    elif folder_of_cryo_fit[0] == "gromacs-4.5.5_cryo_fit_added_mpi":
        print "folder_of_cryo_fit[0] = gromacs-4.5.5_cryo_fit_added_mpi"
        output_bool_enable_mpi = True
    else:
        print "\nInstall cryo_fit first. Refer http://www.phenix-online.org/documentation/reference/cryo_fit.html"
        print "exit now"
        exit(1)
    return output_bool_enable_mpi
    '''
# end of know_output_bool_enable_mpi_by_ls function

'''
def know_cryo_fit_bin_dir_by_options(home_dir, bool_enable_mpi, bool_enable_fftw):
    cryo_fit_bin_dir = '' # initial value
    if (bool_enable_mpi == "True" and bool_enable_fftw == "True"):
      cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit_added_mpi_fftw/bin"
    elif (bool_enable_mpi == "True" and bool_enable_fftw == "False"):
      cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit_added_mpi/bin"
    else:
      cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit_added/bin"
    return cryo_fit_bin_dir
# end of know_cryo_fit_bin_dir function
'''

def know_home_cryo_fit_bin_dir_by_ls():
    home_dir = expanduser("~")
    home_cryo_fit_bin_dir = ''
    command_string = "ls ~/bin | grep gromacs-4.5.5_cryo_fit"
    print "\n\tcommand: ", command_string
    folder_of_cryo_fit = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    #print "\tfolder_of_cryo_fit[0]:", folder_of_cryo_fit[0]
    
    # needed only for debug
    # f_out = open('log.folder_of_cryo_fit', 'wt')
    # f_out.write(folder_of_cryo_fit[0])
    # f_out.close()
    
    if folder_of_cryo_fit[0] == "gromacs-4.5.5_cryo_fit_added":
        #print "\tUser's cryo_fit was installed with enable_mpi=False, so the current cryo_fit will run as enable_mpi = False"
        home_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit_added/bin"
    elif folder_of_cryo_fit[0] == "gromacs-4.5.5_cryo_fit_added_mpi":
        home_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit_added_mpi/bin"
    else:
        print "Install cryo_fit first. Refer http://www.phenix-online.org/documentation/reference/cryo_fit.html"
    return home_cryo_fit_bin_dir
# end of know_output_bool_enable_mpi_by_ls function

def know_home_cryo_fit_bin_dir_by_ls_find():
    home_dir = expanduser("~")
    home_cryo_fit_bin_dir = ''
    command_string = "ls ~/bin | grep gromacs-4.5.5_cryo_fit"
    print "\n\tcommand: ", command_string
    folder_of_cryo_fit = libtbx.easy_run.fully_buffered(command=command_string).raise_if_errors().stdout_lines
    
    if folder_of_cryo_fit[0].find("mpi") == -1:
        print "\tUser's cryo_fit was installed with enable_mpi=False, so the cryo_fit will run as enable_mpi = False"
        home_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit/bin"
    else: # folder_of_cryo_fit[0] == "gromacs-4.5.5_cryo_fit_added_mpi":
        home_cryo_fit_bin_dir = home_dir + "/bin/gromacs-4.5.5_cryo_fit_mpi/bin"
    return home_cryo_fit_bin_dir
# end of know_output_bool_enable_mpi_by_ls_find function


def know_total_number_of_cores():
    if ((platform.system() != "Darwin") and (platform.system() != "Linux")):
        color_print ("User's computer's operating system could be windows")
        number_of_total_cores = 1
        return number_of_total_cores
        
    number_of_total_cores = '' # just initial value
    if (platform.system() == "Darwin"):
        command_string = "sysctl -n hw.ncpu "
        number_of_total_cores = subprocess.check_output(command_string, stderr=subprocess.STDOUT,shell=True)
    elif (platform.system() == "Linux"):
        command_string = "nproc"
        number_of_total_cores = subprocess.check_output(command_string, stderr=subprocess.STDOUT,shell=True)
    else: # maybe Windows
        number_of_total_cores = 2
    
    #color_print ("\tUser's computer's operating system: ", 'green')
    print "\tUser's computer's operating system: "
    print "\t", platform.system(), "\n"
    return number_of_total_cores
# end of know_total_number_of_cores function

def locate_Phenix_executable():
    path = check_output(["which", "phenix.cryo_fit"])
    splited = path.split("/")
    command_path = ''
    for i in range(len(splited)-3):
      command_path = command_path + splited[i] + "/"
    command_path = command_path + "modules/cryo_fit/"
    print "\tUser's phenix.cryo_fit executable comes from ", command_path
    return command_path
# end of locate_Phenix_executable function


def minimize_or_cryo_fit(bool_just_get_input_command, bool_minimization, cores_to_use, \
                              ns_type, common_command_string):
    command_string = '' # just initial
    if (bool_minimization == True and ns_type == "simple"):
        if (bool_enable_mpi == True):
            command_string = "mpirun -np 1 " + common_command_string + " -dd 1 1 1 "
        else:
            command_string = common_command_string + " -nt 1 -dd 1 1 1 "
    else:
        command_string = common_command_string + " -nt " + str(cores_to_use) + " -dd 0 "
        # [keep this comment] for -nt 12, -dd 3 2 2 is needed instead of 2 2 3
    
    if bool_just_get_input_command == False:
        color_print ("\tcommand: ", 'green')
        print "\t", command_string
    #    libtbx.easy_run.call(command=command_string)
    return command_string
# end of minimize_or_cryo_fit function


def remove_former_files():
    current_directory = os.getcwd()
    print "\tRemove former files in ", current_directory
    for each_file in glob.glob("*"):
      if (each_file[:1] == "#") or (each_file[-1:] == "~") or (each_file[-4:] == ".edr") \
        or (each_file == "cryo_fit_log") or (each_file[-4:] == ".log") or (each_file == "md.log") \
        or (each_file[-4:] == ".trr") or (each_file[-4:] == ".xtc") or (each_file == "md.out"):
          subprocess.call(["rm", each_file])
# end of remove_former_files function



def show_time(time_start, time_end):
    time_took = 0 # temporary of course
    if (round((time_end-time_start)/60, 1) < 1):
      time_took = " finished in " + str(round((time_end-time_start), 2)) + " seconds (wallclock)."
    elif (round((time_end-time_start)/60/60, 1) < 1):
      time_took = " finished in " + str(round((time_end-time_start)/60, 2)) + " minutes (wallclock)."
    else:
      time_took = " finished in " + str(round((time_end-time_start)/60/60, 1)) + " hours (wallclock)."
    return time_took
# end of show_time function