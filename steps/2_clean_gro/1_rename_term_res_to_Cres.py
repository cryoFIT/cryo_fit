# How to identify terminal residues
# 1st method: If a residue number decreases, I assume that a new chain starts (so a former residues should be CXXX). It is possible that I need to identify a terminal residue by counting number of atoms, but so far the current concept works well.
# 2nd method: If a residue in gro file has both OC1 and OC2, then it is a terminal one

import glob, os, subprocess, sys

def remove_old_files():
    for each_file in glob.glob("*"):
        if (each_file[:1] == "#") or (each_file[-1:] == "~"):
            subprocess.call(["rm", each_file])
    for to_be_removed_file in glob.glob("*_c_term_renamed*"):
        print "to_be_removed_file: ", to_be_removed_file
        cmd = "rm " + to_be_removed_file
        os.system(cmd)
    for to_be_removed_file in glob.glob("*_temp.gro"):
        print "to_be_removed_file: ", to_be_removed_file
        cmd = "rm " + to_be_removed_file
        os.system(cmd)
################# end of remove_old_files function

if (__name__ == "__main__") :
    #remove_old_files()

    args=sys.argv[1:]
    print "len(args): ", len(args)
    if len(args)<1: # both step 1-2 and 2-3 doesn't need argument, because each folder will have only 1 gro file
        count = 0
        for gro_file in glob.glob("*.gro"):
            input_gro_file_name = gro_file # if there is only 1 gro file in this folder, use it
            count +=1
            if count == 2:
                print "Please specify one input gro file"
                print "example usage: rename.py input.gro"
                sys.exit("rename_term exits now (expecting a gro file at next run)")

        cmd = "python 2_rename_term_res_to_Cres_by_resnum.py " + input_gro_file_name
        os.system(cmd)
        
        new_input_gro_file_name = input_gro_file_name[:-4] + "_c_term_renamed_by_resnum.gro"
        cmd = "python 3_rename_term_res_to_Cres_by_oc.py " + new_input_gro_file_name
        os.system(cmd)
    else:
        command_path = args[0]
        input_gro_file_name = args[1] # gro input file
        command_script = "python " + command_path + "steps/2_clean_gro/2_rename_term_res_to_Cres_by_resnum.py"
        print "command that will be executed: ", command_script
        os.system(command_script)

        new_input_gro_file_name = input_gro_file_name[:-4] + "_c_term_renamed_by_resnum.gro"
        command_script = "python " + command_path + "steps/2_clean_gro/3_rename_term_res_to_Cres_by_oc.py "\
                         + new_input_gro_file_name
        print "command that will be executed: ", command_script
        os.system(command_script)
        