# How to identify terminal residues
# 2nd method: If a residue in gro file has both OC1 and OC2, then it is a terminal one

import glob, os, sys

def decide_HID_HIE_HIP(temp_gro_array):
    his_atom_array = []
    for temp_gro_array_line in temp_gro_array:
        atom_name = str(temp_gro_array_line[12:15])
        his_atom_array.append(atom_name)        
    if his_atom_array[1] == " H1":
        return "CHID"
    elif his_atom_array[9] == "CE1":
        return "CHIE"
    else:
        return "CHIP"
# end of decide_HID_HIE_HIP function

def rename_term_by_OC(input_gro_file_name):
    output_gro_file_name = input_gro_file_name[:-29] + "_c_term_renamed_by_resnum_oc.gro"
    temp_gro_file_name = input_gro_file_name[:-4] + "_temp.gro"
    temp_gro_array = []
    very_first_residue = 0
    
    with open(output_gro_file_name, "a+") as f_out:
        with open(temp_gro_file_name, "w") as f_temp:
            with open(input_gro_file_name) as f_in:
                old_residue_num = -1
                first_OC_found = 0 #initial value
                ignore_the_first_line = 1 # important initial value
                for line in f_in:
                    if (ignore_the_first_line == 1):
                        ignore_the_first_line = 0
                        f_out.write(line) # we need to write the first line anyway
                        continue
                    splited = line.split( )
                    if len(splited) == 6 or len(splited) == 5:
                        residue_num = int(line[1:5])
#                        temp_gro_array.append(line)
                        OC_index = line.find("OC")
                        if (OC_index != -1 and first_OC_found == 0):
                            first_OC_found = 1
                            temp_gro_array.append(line)
                            continue
                        if (OC_index != -1 and first_OC_found == 1):
                            first_OC_found = 0 # reinitialization
                            temp_gro_array.append(line)
                            old_residue_num = residue_num
    #                        print "new chain & residue start"
                            old_residue_num = -1
                            for temp_gro_array_line in temp_gro_array:
                                former_residue_name = str(temp_gro_array_line[5:9]) # to better marshal, use 5:9 rather than 5:8
                                print "former_residue_name: ", former_residue_name
                                print "former_residue_name[0:1]: ", former_residue_name[0:1]
                                if former_residue_name[0:1] == "C": #already renamed with C prefix
                                    f_out.write(temp_gro_array_line)
                                elif former_residue_name == "RG  " or former_residue_name == "G   " or former_residue_name == "C   " or former_residue_name == "A   " or former_residue_name == "T   " or former_residue_name == "U   ": 
                                    # if we are dealing with a nucleic acid, don't add C prefix
                                    f_out.write(temp_gro_array_line)
                                elif former_residue_name == "HIS ": # decide among CHID, CHIE, CHIP
                                    new_former_residue_name = decide_HID_HIE_HIP(temp_gro_array)
                                    new_temp_gro_array_line = temp_gro_array_line.replace(former_residue_name, new_former_residue_name)
                                    f_out.write(new_temp_gro_array_line)
                                else:
                                    new_former_residue_name = "C" + former_residue_name[:3]
                                    new_temp_gro_array_line = temp_gro_array_line.replace(former_residue_name, new_former_residue_name)
                                    f_out.write(new_temp_gro_array_line)

                            #reinitialize temp_gro_array
                            temp_gro_array = []
                            #temp_gro_array.append(line)

                        elif (old_residue_num != residue_num):
                            old_residue_num = residue_num
                            if very_first_residue == 0:
                                temp_gro_array.append(line)
                                very_first_residue = 1
                            else:
       #                         print "new residue starts"
                            
                                # it is safe to write to output_gro as is
                                for array_line in temp_gro_array:
                                    f_out.write(array_line)

                                #reinitialize temp_gro_array
                                temp_gro_array = []
                                temp_gro_array.append(line)
    
                        else:
                            temp_gro_array.append(line)
                            old_residue_num = residue_num
                    else:
                        if len(splited) != 3:
                            f_out.write(line)
                        else: # this is for the last line in gro file (3 sums of each)
                            temp_gro_array.append(line)
                # a very last residue is obviously a c-terminal one, but we already renamed by rename_term_res_to_Cres_by_resnum.py
                for temp_gro_array_line in temp_gro_array:
                    f_out.write(temp_gro_array_line)
    f_in.close()
    f_temp.close()
    f_out.close()
    for to_be_removed_file in glob.glob("*_c_term_renamed_by_resnum.gro"):
        cmd = "rm " + to_be_removed_file
        os.system(cmd)
    for to_be_removed_file in glob.glob("*_temp.gro"):
        cmd = "rm " + to_be_removed_file
        os.system(cmd)
# end of rename_term_by_OC function

if (__name__ == "__main__") :
    args=sys.argv[1:]
    if len(args)<1:
        count = 0
        for gro_file in glob.glob("*.gro"):
            input_gro_file_name = gro_file # if there is only 1 gro file in this folder, use it
            count +=1
            if count == 2:
                print "Please specify one input gro file"
                print "example usage: rename.py input.gro"
                sys.exit("rename_term_res_to_Cres_by_oc.py exits now (expecting a gro file at a next run)")
        rename_term_by_OC(input_gro_file_name)
    else:
        input_gro_file_name = args[0] # gro input file
        rename_term_by_OC(input_gro_file_name)
