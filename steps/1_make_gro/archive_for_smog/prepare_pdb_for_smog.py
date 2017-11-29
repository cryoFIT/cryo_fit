# author: Doonam, Serdal (idea of having TER)

# no more 100k atomnumber
# no more HETATM

import glob, os, sys

def clean_main(input_pdb_file_name):
  output_pdb_file_name = leave_ATOM_END_HETATM_TER (input_pdb_file_name)
  output_pdb_file_name = HETATM_to_ATOM (output_pdb_file_name)
  output_pdb_file_name = no_more_100k_atom_num(output_pdb_file_name)
  output_pdb_file_name = leave_TER_only_if_TER(output_pdb_file_name)
    
  final_output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_smog.pdb"
  cmd = "mv " + output_pdb_file_name + " " + final_output_pdb_file_name
  os.system(cmd)
# end of clean_main function

def leave_ATOM_END_HETATM_TER(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_left_ATOM_HETATM_TER.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    if line[0:4] == "ATOM" or line[0:3] == "END" or line[0:6] == "HETATM" or line[0:3] == "TER":
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name    
# end of leave_ATOM_END_HETATM_TER function


def HETATM_to_ATOM(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_more_HETATM.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    if line[0:6] == "HETATM":
      new_line = "ATOM  " + line[6:]
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of HETATM_to_ATOM function

def leave_TER_only_if_TER(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_left_TER_only_if_TER.pdb"
  f_out = open(output_pdb_file_name, 'wt') 
  for line in f_in:
    if line[0:3] == "TER":
      new_line = "TER\n"
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of leave_TER_only_if_TER function


def remove_former_files():
  for to_be_removed_file in glob.glob("*_cleaned_for_smog*"):
    cmd = "rm " + to_be_removed_file
    os.system(cmd)
# end of remove_former_files function

def no_more_100k_atom_num(input_pdb_file_name):
  # to solve truncation problem of 80S.mdfit_fitted_by_chimera_cleaned_for_gromacs_chain_2.pdb
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_100k_AtomNum.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    if line[0:4] == "ATOM":
      atom_number = line[5:11]
      if (atom_number == " *****"):
        f_out.write(line)
      else:
        if int(atom_number) >= 200000:
          put_void = 6-len(str(int(atom_number)-200000))
          new_line = line[:5] + " "*put_void + str(int(atom_number)-200000) + line[11:]
          f_out.write(new_line)
        elif int(atom_number) >= 100000:
          put_void = 6-len(str(int(atom_number)-100000))
          new_line = line[:5] + " "*put_void + str(int(atom_number)-100000) + line[11:]
          f_out.write(new_line)
        else:
          f_out.write(line)
    else: # line[0:6] == "HETATM" or line[0:3] == "TER"
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of no_more_100k_atom_num function

if (__name__ == "__main__") :
  remove_former_files()

  print "\tInput format: python clean_pdb_for_gromacs.py input_pdb_file_name"
  print "\tExample usage: python clean_pdb_for_gromacs.py 80S.CS.nomag_MIA.pdb"
  args=sys.argv[1:]
  
  input_pdb_file_name = args[0] # pdb input file
  clean_main(input_pdb_file_name)
  '''
  
  if (len(args) >= 1):
    input_pdb_file_name = args[0] # pdb input file
    if (input_pdb_file_name[-4:] != ".pdb"):
      print "Entered input_pdb_file is not in pdb format"
      exit(1)
  else:
    print "provide pdb_file"
    exit(1)
  '''
  