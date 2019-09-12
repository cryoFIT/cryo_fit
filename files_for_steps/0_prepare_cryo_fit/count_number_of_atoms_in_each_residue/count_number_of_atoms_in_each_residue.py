# author: Doonam

import glob, os, sys

def count_number_of_atoms_in_each_residue(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_file_name = input_pdb_file_name[:-4] + "_w_atom_numbers.txt"
  
  if os.path.isfile(output_file_name) == True:
    cmd = "rm " + output_file_name
    os.system(cmd)
  
  f_out = open(output_file_name, "a")
  
  accumulated_atom_num = 0
  old_RES = ''
  old_chain = ''
  old_resnum = ''
  first_ATOM_HETATM_line = True
  line_num = 0
  for line in f_in:
    RES = line[17:20]
    chain = line[21:22]
    resnum = line[23:26]
    
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      #print line_num
      line_num += 1
      
      if (first_ATOM_HETATM_line == True):
        first_ATOM_HETATM_line = False
        old_RES = RES
        old_chain = chain
        old_resnum = resnum
        accumulated_atom_num += 1
        continue
      
      if ((old_RES == RES) and (old_chain == chain) and (old_resnum == resnum)):
        accumulated_atom_num += 1
        continue
      else:
        #write_this = str(line_num) + " " + str(old_RES) + " " + str(old_chain) + " " + str(old_resnum) + " " + str(accumulated_atom_num-1) + "\n"
        write_this = str(old_RES) + " " + str(old_chain) + " " + str(old_resnum) + " " + str(accumulated_atom_num) + "\n"
        f_out.write(write_this)
        accumulated_atom_num = 0
      old_RES = RES
      old_chain = chain
      old_resnum = resnum
    
      
  f_in.close()
  f_out.close()
################ end of count_number_of_atoms_in_each_residue function

if (__name__ == "__main__") :
  args=sys.argv[1:]
  if len(args)<1:
    print "Please specify one input PDB file"
    print "example usage: python divide_pdb_file_by_chains.py input.pdb"
    sys.exit("exit now")
    
  else: # user provided one pdb input file
    input_pdb_file_name=args[0] # pdb input file
    count_number_of_atoms_in_each_residue(input_pdb_file_name)
