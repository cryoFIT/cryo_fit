# this is needed to prevent initial cell size smaller error
import glob, os, sys

def slightly_change(line, old_line, old_old_line, perturb_xyz_by):
    x_coor_former = old_line[20:28]
    if x_coor_former == "   0.000":
        x_coor_former = old_old_line[20:28]
        perturb_xyz_by = float(perturb_xyz_by) + 0.04 # based on other normal double hydrogens, \
                                        # adding 0.05 seems appropriate (not subtracting),
                                        # 0.05 worked fine at lanl laptop, but not enough at kaguya
    new_x_coor = str(float(x_coor_former) + float(perturb_xyz_by))
    splited = new_x_coor.split(".")
    multi_before_period = 4-len(splited[0])
    multi_after_period = 3-len(splited[1])
    new_line = line[:20] + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" "
    
    y_coor_former = old_line[28:36]
    if y_coor_former == "   0.000":
        y_coor_former = old_old_line[28:36]
    new_y_coor = str(float(y_coor_former) + float(perturb_xyz_by))
    splited = new_y_coor.split(".")
    multi_before_period = 4-len(splited[0])
    multi_after_period = 3-len(splited[1])
    new_line = new_line + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" "
    
    z_coor_former = old_line[36:44]
    if z_coor_former == "   0.000":
        z_coor_former = old_old_line[36:44]
    
    new_z_coor = str(float(z_coor_former) + float(perturb_xyz_by))
    splited = new_z_coor.split(".")
    multi_before_period = 4-len(splited[0])
    multi_after_period = 3-len(splited[1])
    new_line = new_line + multi_before_period*" "+splited[0] + "." + splited [1]+multi_after_period*" " \
          + line[44:]
    return new_line
############## end of slightly_change function

def no_more_000_xyz(input_gro_file_name, perturb_xyz_by):
    f_in = open(input_gro_file_name)
    output_gro_file_name = input_gro_file_name[:-4] + "_no_more_000_xyz.gro"
    f_out = open(output_gro_file_name, 'wt')
    old_old_line = '' # initial value
    old_line = '' # initial value
    for line in f_in:
        x = line[20:28]
        y = line[28:36]
        z = line[36:44]
        if (x == "   0.000" and y == "   0.000" and z == "   0.000"):
            new_line = slightly_change(line, old_line, old_old_line, perturb_xyz_by)
            f_out.write(new_line)
        else:
            f_out.write(line)
        old_old_line = old_line
        old_line = line
    f_in.close()
    f_out.close()
    cmd = "rm " + input_gro_file_name
    os.system(cmd)
    return output_gro_file_name
####################### end of no_more_000_xyz function

if (__name__ == "__main__") :
    args=sys.argv[1:]
    if len(args)<1:
        count = 0
        for gro_file in glob.glob("*.gro"):
            input_gro_file_name = gro_file # if there is only 1 gro file in this folder, use it
            count +=1
            if count == 2:
                print "Please specify one input gro file"
                print "example usage: rename.py input.gro"
                sys.exit("rename_term_res_to_Cres_by_resnum.py exits now (expecting a gro file at next run)")
        no_more_000_xyz(input_gro_file_name)
    else: # is used by cryo_fit python code
        input_gro_file_name = args[0] # gro input file
        perturb_xyz_by = args[1]
        no_more_000_xyz(input_gro_file_name, perturb_xyz_by)
