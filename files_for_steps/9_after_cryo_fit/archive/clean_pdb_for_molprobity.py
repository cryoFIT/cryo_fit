# author: Doonam
# purpose: confirmed to be essential for molprobity to run with Dieter's pdb file

import glob, os, sys

def clean_main(input_pdb_file_name):
  output_pdb_file_name = change_OC12 (input_pdb_file_name) # I remove OXT in user's pdb for a pdb file like Dieter's that has OXT before "TER"
  
  final_output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_molprobity.pdb"
  cmd = "mv " + output_pdb_file_name + " " + final_output_pdb_file_name
  os.system(cmd)
  
  return final_output_pdb_file_name
# end of clean_main function

def change_OC12(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_wo_OC12.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    atom = line[13:16]
    if atom == "OC1":
      new_line = line[:13] + "OXT" + line[16:]
      f_out.write(new_line)
    elif atom == "OC2":
      print "\t\t\tOC2 omitted"
      continue
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of change_OC12 function

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
