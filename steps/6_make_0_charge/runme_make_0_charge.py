import glob, os, sys

def make_0_charge(input_top_file_name):
    output_top_name = input_top_file_name[:-4] + "_0_charge.top"
    cmd = "awk -f changetop.awk < " + input_top_file_name + " > " + output_top_name
    os.system(cmd)

if (__name__ == "__main__") :
  args=sys.argv[1:]
  if len(args)<1:
    count = 0
    for top_file in glob.glob("*.top"):
      input_top_file_name = top_file # if there is only 1 top file in this folder, use it
      count +=1
      if count == 2:
        print "Please specify one input top file"
        print "example usage: runme_make_0_charge.py input.top"
        sys.exit("runme_make_0_charge exits now (expecting a top file at next run)")
    make_0_charge(input_top_file_name)
  else:
    input_top_file_name=args[0] # pdb input file
    make_0_charge(input_top_file_name)
