# author: Doonam


import glob, os, sys

def clean_main(input_pdb_file_name):
  print "\n\t\tFor realspace_refine, remove CRYST line to avoid error \"Crystal symmetry mismatch between different files\""
  print "\t\t\t\t      U instead of RU to avoid error"
  print "\t\t\t\t      add element at the end of each line"
  print "\n\t\tFor molprobity, change OC1 and OC2 to avoid error"
  print "\t\t\t\t      Omit MODEL and ENDMDL lines to avoid error"
  print "\t\t\tOmit hydrogen since it is not essential to run real_space_refine and molprobity"

  output_pdb_file_name = clean (input_pdb_file_name) 
  
  final_output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_real_space_refine_molprobity.pdb"
  cmd = "mv " + output_pdb_file_name + " " + final_output_pdb_file_name
  os.system(cmd)
  
  return final_output_pdb_file_name
##################################### end of clean_main function


def deal_OC1(atom, new_line):
  if atom == "OC1":
    new_line = new_line[:13] + "OXT" + new_line[16:]
  return new_line 
##################################### end of def deal_OC1(atom)


def clean(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    CRYST_MODEL = line[:5]
    ENDMDL = line[:6]
    element = line[13:14]
    atom = line[13:16]
    
    if ((CRYST_MODEL == "CRYST") or (CRYST_MODEL == "MODEL") or (ENDMDL == "ENDMDL")): # remove CRYST for real_space_refine, remove MODEL for molprobity
      print "\t\t\t omitted ", line
      continue
    elif atom == "OC2":
      print "\t\t\t omitted OC2 for molprobity running"
      continue
    elif (line[0:4] != "ATOM") and (line[0:6] != "HETATM"):
      f_out.write(line)
    else: # most cases
      if (element == "H"):
        print "\t\t\t omiited hydrogen since it is not essential to run real_space_refine and molprobity"
        continue
      new_line = line[:75] + "  " + element + line[79:] + "\n"
      new_line = new_line[:17] + new_line[17:20].replace(' RA', ' A ') + new_line[20:] #residue = line[17:20]
      new_line = new_line[:17] + new_line[17:20].replace(' RU', ' U ') + new_line[20:]
      new_line = new_line[:17] + new_line[17:20].replace(' RG', ' G ') + new_line[20:]
      new_line = new_line[:17] + new_line[17:20].replace(' RC', ' C ') + new_line[20:]
      new_line = deal_OC1(atom, new_line)
      f_out.write(new_line)
      
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
##################################### end of clean function

if (__name__ == "__main__") :
  args=sys.argv[1:]
  if (len(args) >= 1):
    input_pdb_file_name = args[0] # pdb input file
    if (input_pdb_file_name[-4:] != ".pdb"):
      print "Entered input_pdb_file is not in pdb format."
      exit(1)
  else:
    print "Please provide a pdb file."
    exit(1)
  
  clean_main(input_pdb_file_name)