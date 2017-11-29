import os, subprocess, sys

if len (sys.argv) < 3:
   print "How to use: replace_xyz_from_gro.py file.pdb file.gro"
   print "Example: replace_xyz_from_gro.py ../1_make_top/80S_fit_by_chi_to_1_manual_fixed_trp.pdb dump_10000_steps_5.0_ps.gro"
   sys.exit(0)

input_pdb_file = sys.argv[1]
input_gro_file = sys.argv[2]

if os.path.isfile("pdb_indel.pl") == False:
   cmd = "cp /Users/doonam/Dropbox/research/lanl/computation/run/script/perl/replace_xyz_from_gro/pdb_indel.pl ."
   os.system(cmd)

chain_recovered_needs_marshal = input_gro_file[:-4] + "_chain_recovered_needs_marshal.pdb"

if os.path.isfile(chain_recovered_needs_marshal) == True:
   cmd = "rm " + chain_recovered_needs_marshal
   os.system(cmd)
   
cmd = "perl gro_to_pdb_dn.pl " + input_pdb_file + " " + input_gro_file + " > " + chain_recovered_needs_marshal
os.system(cmd)

chain_recovered = input_gro_file[:-4] + "_chain_recovered.pdb"

if os.path.isfile(chain_recovered) == True:
   cmd = "rm " + chain_recovered
   os.system(cmd)
   
f_in = open(chain_recovered_needs_marshal)
f_out = open(chain_recovered, "a")
for line in f_in:
   if line[:4] == "ATOM" or line[:6] == "HETATM":
      new_line = line[:30] + line[31:]
      #print line
      #print new_line
      f_out.write(new_line)
   else:
      f_out.write(line)
f_in.close()
f_out.close()

cmd = "rm " + chain_recovered_needs_marshal
os.system(cmd)
