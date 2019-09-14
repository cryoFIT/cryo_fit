# author: Doonam

import glob, os, sys

def clean_main(input_pdb_file_name):
  print "\n\tFor realspace_refine, remove CRYST line to avoid error \"Crystal symmetry mismatch between different files\""
  print "\t\t  U instead of RU to avoid error"
  print "\t\t  add element at the end of each line"
  print "\n\tFor molprobity, change OC1 and OC2 to avoid error"
  print "\t\t  Omit MODEL and ENDMDL lines to avoid error"

  output_pdb_file_name = clean (input_pdb_file_name) 
  
  final_output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_real_space_refine_molprobity.pdb"
  cmd = "mv " + output_pdb_file_name + " " + final_output_pdb_file_name
  os.system(cmd)
  
  #output_pdb_file_name = clean_for_chimeraX (input_pdb_file_name) 
##################################### end of clean_main function


def ILE_CD_to_CD1(residue, atom, new_line):
  if residue == "ILE":
    if atom == "CD ":
      new_line = new_line[:13] + "CD1" + new_line[16:]
  return new_line 
##################################### end of ILE_CD_to_CD1(residue, atom, new_line)


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
    residue = line[17:20]
    
    if ((CRYST_MODEL == "CRYST") or (CRYST_MODEL == "MODEL") or (ENDMDL == "ENDMDL")): # remove CRYST for real_space_refine, remove MODEL for molprobity
      #print "\t\t\t omitted ", line # print unnecessary \n at the end
      print "\t\t\t omitted ", line[:len(line)-2]
      continue
    elif atom == "OC2":
      print "\t\t\t omitted OC2 for molprobity running"
      continue
    elif (line[0:4] != "ATOM") and (line[0:6] != "HETATM"):
      f_out.write(line)
    else: # most cases
      #8/15/2018, I confirmed that having hydrogen from cryo_fit is totally OK for following real_space_refine and molprobity
      # keeping hydrogen is essential for movie showing in ChimeraX
      #if (element == "H"):
      #  print "\t\t\t omiited hydrogen since it is not essential to run real_space_refine and molprobity"
      #  continue
      new_line = line[:75] + "  " + element + line[79:] + "\n"
      new_line = new_line[:17] + new_line[17:20].replace(' RA', ' A ') + new_line[20:] #residue = line[17:20]
      new_line = new_line[:17] + new_line[17:20].replace(' RU', ' U ') + new_line[20:]
      new_line = new_line[:17] + new_line[17:20].replace(' RG', ' G ') + new_line[20:]
      new_line = new_line[:17] + new_line[17:20].replace(' RC', ' C ') + new_line[20:]
      
      new_line = ILE_CD_to_CD1(residue, atom, new_line) 
      # http://www.phenix-online.org/pipermail/phenixbb/2015-October/022450.html
      # For RSR of cryo_fitted adenylate kinase, I need this function
      
      new_line = deal_OC1(atom, new_line)
      f_out.write(new_line)
      
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
##################################### end of clean function


'''
def clean_for_chimeraX(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_chimeraX.pdb"
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
      #8/15/2018, I confirmed that having hydrogen from cryo_fit is totally OK for following real_space_refine and molprobity
      # keeping hydrogen is essential for movie showing in ChimeraX
      #if (element == "H"):
      #  print "\t\t\t omiited hydrogen since it is not essential to run real_space_refine and molprobity"
      #  continue
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
##################################### end of clean_for_chimeraX function
'''


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
