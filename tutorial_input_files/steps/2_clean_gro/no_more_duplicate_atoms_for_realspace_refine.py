# author: Doonam
# purpose: fully run this cleaning code (no more duplicates) for ribosome
# In total, it took 3 hours with 16GM RAM macbookpro

import glob, os, sys

def check_duplicate_atoms(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  bool_same_atom = False
  atom_list = []
  res_name_list = []
  chain_list = []
  res_num_list = []
  for line in f_in:
    print "line during checking:", line
    atom = line[12:16]
    res_name = line[17:20]
    chain = line[21:22]
    res_num = line[22:27]
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      for i in range(len(chain_list)):
        if atom == atom_list[i] and res_name == res_name_list[i] and chain == chain_list[i] \
            and res_num == res_num_list[i]:
          f_in.close()
          bool_same_atom = True
          return bool_same_atom
      atom_list.append(atom)
      res_name_list.append(res_name)
      chain_list.append(chain)
      res_num_list.append(res_num)
      
  f_in.close()
  return bool_same_atom
# end of check_duplicate_atoms function

def no_more_duplicate_atoms(input_pdb_file_name):
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_more_duplicate_atoms.pdb"
  f_in = open(input_pdb_file_name)
  f_out = open(output_pdb_file_name, "w")
  atom_list = []
  res_name_list = []
  chain_list = []
  res_num_list = []
  for line in f_in:
    print "line during changing:", line
    atom = line[12:16]
    res_name = line[17:20]
    chain = line[21:22]
    res_num = line[22:27]
    bool_same_atom = False
    if line[0:3] == "END" or line[0:3] == "TER":
      f_out.write(line)
    elif line[0:4] == "ATOM" or line[0:6] == "HETATM" or line[0:6] == "SPLICE":
      # print "atom: ", atom
      # print "res_name: ", res_name
      # print "chain: ", chain
      # print "res_num: ", res_num
      
      for i in range(len(chain_list)):
        if atom == atom_list[i] and res_name == res_name_list[i] and chain == chain_list[i] \
            and res_num == res_num_list[i]:
          bool_same_atom = True
          print "same atom"
          print "atom_list[i]: ", atom_list[i]
          print "res_name_list[i]: ", res_name_list[i]
          print "chain_list[i]: ", chain_list[i]
          print "res_num_list[i]: ", res_num_list[i]
          print "\n"
          new_line = line[:22] + str(int(res_num)+1000) + line[26:]
          print "new_line:", new_line
          f_out.write(new_line)
          break
      if (bool_same_atom == False):
        f_out.write(line)
      bool_same_atom = False # reinitialization

      atom_list.append(atom)
      res_name_list.append(res_name)
      chain_list.append(chain)
      res_num_list.append(res_num)
      
  f_in.close()
  f_out.close()
  return output_pdb_file_name
# end of no_more_duplicate_atoms function

if (__name__ == "__main__") :
  args=sys.argv[1:]
  input_pdb_file_name=args[0] # pdb input file
  
  # one time execution
  bool_same_atom = check_duplicate_atoms(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name
  
  while (bool_same_atom == True):
    output_pdb_file_name = no_more_duplicate_atoms(output_pdb_file_name)
    bool_same_atom = check_duplicate_atoms(output_pdb_file_name)
  
  cmd = "mv " + output_pdb_file_name + " " + input_pdb_file_name[:-4] + "_no_more_duplicate_atoms.pdb"
  os.system(cmd)
