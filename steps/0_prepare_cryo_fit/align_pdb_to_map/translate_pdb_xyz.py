import glob, os, sys

def translate_xyz(input_pdb_file_name, move_x_by, move_y_by, move_z_by):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_translated" + ".pdb"
  f_out = open(output_pdb_file_name, "w")
  for line in f_in:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      x_coor_former = line[30:38]
      new_x_coor = str(float(x_coor_former) + float(move_x_by))
      splited = new_x_coor.split(".")
      multi_before_period = 4-len(splited[0])
      multi_after_period = 3-len(splited[1])
      new_line = line[:30] + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" "
      
      y_coor_former = line[38:46]
      new_y_coor = str(float(y_coor_former) + float(move_y_by))
      splited = new_y_coor.split(".")
      multi_before_period = 4-len(splited[0])
      multi_after_period = 3-len(splited[1])
      new_line = new_line + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" "
      
      z_coor_former = line[46:54]
      new_z_coor = str(float(z_coor_former) + float(move_z_by))
      splited = new_z_coor.split(".")
      multi_before_period = 4-len(splited[0])
      multi_after_period = 3-len(splited[1])
      new_line = new_line + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" " \
            + line[54:]
      f_out.write(new_line)
      
    elif line[0:3] == "TER":
      f_out.write(line)
  f_in.close()
  f_out.close()
        

if (__name__ == "__main__") :
  print "Input format: python translate_pdb_xyz.py input_pdb_file_name x y z\n"
  print "Example usage: python translate_pdb_xyz.py 80S.CS.nomag_MIA.pdb 52.090 -75.090 -15.360\n"
  args=sys.argv[1:]
  if len(args)<1:
    count = 0
    for pdb_file in glob.glob("*.pdb"):
      input_pdb_file_name = pdb_file # if there is only 1 pdb file in this folder, use it
      count +=1
      if count == 2:
        print "Please specify one input PDB file"
        print "example usage: python translate_pdb_xyz.py input.pdb"
        sys.exit("translate_pdb_xyz exits now (expecting a pdb file at next run)")
    divide_by_chains(input_pdb_file_name)
  else: # user provided one pdb input file
    input_pdb_file_name=args[0] # pdb input file
    move_x_by=args[1]
    move_y_by=args[2]
    move_z_by=args[3] 
    translate_xyz(input_pdb_file_name, move_x_by, move_y_by, move_z_by)
    
