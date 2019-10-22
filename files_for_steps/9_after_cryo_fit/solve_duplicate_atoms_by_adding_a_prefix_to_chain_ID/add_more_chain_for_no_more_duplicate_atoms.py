# author: Doonam

import glob, os, sys
   
def change_chain_if_duplicated(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_more_duplicate_atoms_w_a_prefix_to_chain_name.pdb"
  f_out = open(output_pdb_file_name, "w")
  chain_array = []
  chain_big_array = [] # to know whether a chain starts newly
  bool_first_line = True
  for line in f_in:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      #print "line:", line
      bool_chain_just_added_to_array = False
      chain = line[21:22]
      chain_exists_in_chain_array = chain in chain_array
      #if (chain_exists_in_chain_array == True):
      #    print "This chain already exists in chain_array"
      if bool_first_line == True:
        old_chain = chain
        bool_first_line = False
        f_out.write(line)
      else: #bool_first_line == False
        #print "old_chain:", old_chain
        #print "chain:", chain
        chain_exists_in_big_array = chain in chain_big_array
        if (chain_exists_in_big_array == True):
          #print "This chain already exists in chain_big_array"
          #print "chain_exists_in_big_array:", chain_exists_in_big_array
        #if chain_exists_in_big_array == True:
          new_line = line[:20] + "a" + line[21:]
          #print "line[:20]:", line[:20]
          #print "line[21:]:", line[21:]
          #print "new_line:", new_line
          f_out.write(new_line)
        elif old_chain != chain:
          #print "old_chain != chain"
          f_out.write(line)
          chain_big_array.append(old_chain) # chain 2 starts after many chain 1 atoms
          #print "chain_big_array:", chain_big_array
        else:
          f_out.write(line)
        old_chain = chain
      #print "chain:", chain
      try:
          #print "tried"
          chain_array.index(chain)
      except: # a new chain met
          chain_array.append(chain)
    else:
      f_out.write(line)
    
    #print "\n\n"  
  f_in.close()
  f_out.close()
# end of change_chain_if_duplicated function

if (__name__ == "__main__") :
   if len (sys.argv) < 2:
      print "How to use: python add_more_chain.py file.pdb "
      print "Example:    python add_more_chain.py 80S_fit_by_chi_to_1_manual_fixed_trp.pdb"
      sys.exit(0)
      
   args=sys.argv[1:]
   input_pdb_file_name=args[0] # pdb input file
   
   change_chain_if_duplicated(input_pdb_file_name)
   
