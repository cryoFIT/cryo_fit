# first make .top file from your .pdb file at http://davapc1.bioch.dundee.ac.uk/prodrg

import glob, os, platform, subprocess, sys, time
from subprocess import check_output, Popen, PIPE
from os.path import expanduser # to find home_dir

termcolor_installed = '' # just initial value
try:
    from termcolor import colored
    termcolor_installed = True
except Exception:
    termcolor_installed = False
    print "Your computer has no termcolor"
    print "If you want to see cryo_fit install helper's comments in color, download it from https://pypi.python.org/pypi/termcolor"
    print "Then run \"python setup.py install\" in the extracted folder"
    print "Press any key to continue"
    raw_input()

def color_print(text, color):
    if (termcolor_installed == True):
        print colored (text, color)
    else:
        print text

def show_time(process, time_start, time_end):
  time_took = 0 # temporary of course
  if (round((time_end-time_start)/60, 1) < 1):
    time_took = process + " finished in " + str(round((time_end-time_start), 2)) + " seconds (wallclock)."
  elif (round((time_end-time_start)/60/60, 1) < 1):
    time_took = process + " finished in " + str(round((time_end-time_start)/60, 2)) + " minutes (wallclock)."
  else:
    time_took = process + " finished in " + str(round((time_end-time_start)/60/60, 1)) + " hours (wallclock)."
  return time_took

def top2rtp(top_file_name):
    f_in = open(top_file_name, 'r')
    output_rtp_file_name = top_file_name[:-4] + "_top_by_prodrg.rtp"
    f_out = open(output_rtp_file_name, 'w')
    use_next_next_line_as_name = 0
    write_as_atom = 0
    write_as_bonds = 0
    write_as_angles = 0
    write_as_dihedrals = 0
    improper_name_written = 0
    dihedrals_name_written = 0
    for line in f_in:
        line = line[:-1] # remove \n
        splited = line.split()
        if (len(splited) < 2):
            continue
        if (splited[0] == ";"):
            continue
        if (splited[1] == "moleculetype"):
            use_next_next_line_as_name = 1
            continue
        if (use_next_next_line_as_name == 1):
            use_next_next_line_as_name = 0
            cmd = "[ " + splited[0] + " ] ; added by top2rtp"
            f_out.write(cmd)
            f_out.write("\n")
            continue
        if (splited[1] == "atoms"):
            cmd = " [ " + splited[1] + " ]"
            f_out.write(cmd)
            f_out.write("\n")
            write_as_atom = 1
            continue
        if (write_as_atom == 1):
            if (splited[1] == "bonds"):
                cmd = " [ " + splited[1] + " ]"
                f_out.write(cmd)
                f_out.write("\n")
                write_as_atom = 0
                write_as_bonds = 1
                continue
            cmd = " \t" + splited[4] + "\t" + splited[1] + "\t" + splited[6] + "\t" + splited[0]
            f_out.write(cmd)
            f_out.write("\n")
        if (write_as_bonds == 1):
            if (splited[1] == "pairs"):
                write_as_bonds = 0
                continue
            cmd = " \t" + splited[8] + "\t" + splited[9]
            f_out.write(cmd)
            f_out.write("\n")
        if (splited[1] == "angles"):
            cmd = " [ " + splited[1] + " ]"
            f_out.write(cmd)
            f_out.write("\n")
            write_as_angles = 1
            continue
        if (write_as_angles == 1):
            if (splited[1] == "dihedrals"):
                write_as_angles = 0
                write_as_dihedrals = 1
                continue
            cmd = " \t" + splited[9] + "\t" + splited[10] + "\t" + splited[11]
            f_out.write(cmd)
            f_out.write("\n")
        if (write_as_dihedrals == 1):
            if (splited[10] == "imp"):
                if (improper_name_written == 0):
                    cmd = " [ impropers ]"
                    f_out.write(cmd)
                    f_out.write("\n")
                    improper_name_written = 1
                cmd = " \t" + splited[11] + "\t" + splited[12] + "\t" + splited[13] + "\t" + splited[14]
                f_out.write(cmd)
                f_out.write("\n")
            else:
                if (dihedrals_name_written == 0):
                    cmd = " [ dihedrals ]"
                    f_out.write(cmd)
                    f_out.write("\n")
                    dihedrals_name_written = 1
                cmd = " \t" + splited[13] + "\t" + splited[14] + "\t" + splited[15] + "\t" + splited[16]
                f_out.write(cmd)
                f_out.write("\n")
    f_in.close()
    f_out.close()
            

if (__name__ == "__main__") :
    start_time = time.time()
    args=sys.argv[1:]
    if len(args)<1:
        color_print ("How to use   : python runme_top2rtp.py <.top file>", 'green')
        color_print ("Example usage: python runme_top2rtp.py SEP_prodrg.top", 'green')
        color_print ("runme_top2rtp.py exits now", 'green')
        exit(1)
    else:
        top_file_name = args[0] # pdb input file
        color_print ("input .top file: ", 'green')
        print top_file_name
        top2rtp(top_file_name)
    end_time = time.time()
    show_time("top2rtp", start_time, end_time)
