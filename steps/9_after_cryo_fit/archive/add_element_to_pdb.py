# author: Doonam
# purpose: confirmed to be essential for molprobity to run with Dieter's pdb file
import glob, os, sys

def add_new_element(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_new_element_added.pdb"
  if os.path.isfile(output_pdb_file_name) == True:
    cmd = "rm " + output_pdb_file_name
    os.system(cmd)
  f_out = open(output_pdb_file_name, "a")
  for line in f_in:
    element = line[13:14]
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      new_line = line[:75] + "  " + element + line[79:] + "\n"
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
# end of add_new_element function

if (__name__ == "__main__") :
  #print "input format: python divide_pdb_file_by_chains.py input_pdb_file_name \n"
  #print "example usage: python divide_pdb_file_by_chains.py 80S.CS.nomag_MIA.pdb \n"
  args=sys.argv[1:]
  if len(args)<1:
    count = 0
    for pdb_file in glob.glob("*.pdb"):
      input_pdb_file_name = pdb_file # if there is only 1 pdb file in this folder, use it
      count +=1
      if count == 2:
        print "Please specify one input PDB file"
        print "example usage: python divide_pdb_file_by_chains.py input.pdb"
        sys.exit("divide_pdb_file_by_chains exits now (expecting a pdb file at next run)")
    add_new_element(input_pdb_file_name)
  else: # user provided one pdb input file
    input_pdb_file_name=args[0] # pdb input file
    add_new_element(input_pdb_file_name)
