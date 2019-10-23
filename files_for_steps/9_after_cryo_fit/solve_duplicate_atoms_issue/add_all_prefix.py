# author: Doonam

import glob, os, random, string, sys

   
def change_chain_if_duplicated(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_duplicate_atoms_w_prefix_to_chain_ID.pdb"
  f_out = open(output_pdb_file_name, "w")
  
  two_char_chain_array = []
  
  bool_very_first_line = True
  old_chain       = '' # just initial
  old_atom_number = '' # just initial
  new_two_char_chain = '' # just initial
  
  all_single_char = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"#$%&\'()*+,-./:;<=>?@[]^_`{|}~ '
  all_single_char_list = list(all_single_char)
  new_two_char_written = False # initial
  
  for line in f_in:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      #print line
      chain = line[21:22]
      two_char_chain = line[20:22]
      atom_number = line[22:26]
      
      if bool_very_first_line == True:
        bool_very_first_line = False
        f_out.write(line)
        two_char_chain_array.append(two_char_chain)
      else: # bool_very_first_line = False
        
        if (atom_number < old_atom_number):
          two_char_chain_exists_in_array = two_char_chain in two_char_chain_array
          if (two_char_chain_exists_in_array == False):
            two_char_chain_array.append(two_char_chain)
            f_out.write(line)
            new_two_char_written = False
          else:
            while (1):
              single_prefix = random.choice(all_single_char_list)
              new_two_char_chain = str(single_prefix) + str(chain)
              new_two_char_chain_exists_in_array = str(new_two_char_chain) in two_char_chain_array
              if (new_two_char_chain_exists_in_array == False):
                break
              
            two_char_chain_array.append(new_two_char_chain)
            new_line = line[:20] + str(single_prefix) + line[21:]
            f_out.write(new_line)
            new_two_char_written = True
            
        elif old_chain != chain:
          two_char_chain_exists_in_array = str(two_char_chain) in two_char_chain_array
          if (two_char_chain_exists_in_array == False):
            two_char_chain_array.append(two_char_chain)
            f_out.write(line)
            new_two_char_written = False
          else:
            while (1):
              single_prefix = random.choice(all_single_char_list)
              new_two_char_chain = str(single_prefix) + str(chain)
              new_two_char_chain_exists_in_array = str(new_two_char_chain) in two_char_chain_array
              if (new_two_char_chain_exists_in_array == False):
                break
           
            two_char_chain_array.append(new_two_char_chain)
            new_line = line[:20] + str(single_prefix) + line[21:]
            f_out.write(new_line)
            new_two_char_written = True
        
        elif (new_two_char_written == True):
          new_line = line[:20] + str(new_two_char_chain) + line[22:]
          f_out.write(new_line)
        else:
          f_out.write(line)
      old_chain = chain
      old_atom_number = atom_number
      
    else:
      f_out.write(line)

  f_in.close()
  f_out.close()
########## end of change_chain_if_duplicated function


if (__name__ == "__main__") :
   if len (sys.argv) < 2:
      print "How to use: python add_more_chain.py <user>.pdb "
      print "Example:    python add_more_chain.py extracted_0.4_ps.pdb"
      sys.exit(0)
      
   args=sys.argv[1:]
   input_pdb_file_name=args[0] # pdb input file
   
   change_chain_if_duplicated(input_pdb_file_name)
