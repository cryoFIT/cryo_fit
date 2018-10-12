# this is needed to prevent initial cell size smaller error
import glob, os, sys

def remove_these_residues (input_gro_file_name, residues_to_be_removed):
    f_in = open(input_gro_file_name)
    output_gro_file_name = input_gro_file_name[:-4] + "_removed_residues_w_000_xyz.gro"
    f_out = open(output_gro_file_name, 'wt')
    for line in f_in:
        splited = line.split()
        #result = residues_to_be_removed.index(str(splited[0]))
        try:
            result = residues_to_be_removed.index(str(splited[0]))
        except:
            f_out.write(line)
    f_in.close()
    f_out.close()
    cmd = "rm " + input_gro_file_name
    os.system(cmd)
    return output_gro_file_name
# end of remove_these_residues function
    
def identify_000_xyz_residues(input_gro_file_name):
    f_in = open(input_gro_file_name)
    residues_to_be_removed = []
    old_residue = '' #initial value
    first_atom_identifier = '' # initial value
    for line in f_in:
        splited = line.split()
        #residue = splited[0]
        # if (residue != old_residue):
        #     first_atom_identifier = splited[1]
        x = line[20:28]
        y = line[28:36]
        z = line[36:44]
        if (x == "   0.000" and y == "   0.000" and z == "   0.000"):
            residues_to_be_removed.append(splited[0])
        #    first_atom_identifier.append(splited[1])
        #old_residue = residue
    f_in.close()
    return residues_to_be_removed
# end of identify_000_xyz_residues function

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
                sys.exit("rename_term_res_to_Cres_by_resnum.py exits now (expecting a gro file at next run)")
        residues_to_be_removed = identify_000_xyz_residues(input_gro_file_name)
        remove_these_residues (input_gro_file_name, residues_to_be_removed)
    else: # is used by cryo_fit python code
        input_gro_file_name = args[0] # gro input file
        residues_to_be_removed = identify_000_xyz_residues(input_gro_file_name)
        remove_these_residues (input_gro_file_name, residues_to_be_removed)
