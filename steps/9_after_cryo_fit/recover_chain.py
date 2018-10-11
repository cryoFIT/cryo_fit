#author:doonam
import os, subprocess, sys, time

if len (sys.argv) < 3:
   print "How to use: recover_chain.py pdb_in_step1.pdb cryo_fitted.pdb"
   print "Example:    recover_chain.py ../1_make_top/cleaned_for_gromacs.pdb cryo_fitted.pdb"
   sys.exit(0)

pdb_before_cryo_fit = sys.argv[1]
pdb_after_cryo_fit = sys.argv[2]

chain_recovered = pdb_after_cryo_fit[:-4] + "_chain_recovered.pdb"

if os.path.isfile(chain_recovered) == True:
   cmd = "rm " + chain_recovered
   os.system(cmd)

pdb_after_cryo_fit_in = open(pdb_after_cryo_fit)

f_out = open(chain_recovered, 'w')

def retrieve_chain(seeking_res, seeking_res_num, search_after_this_line_num_in_ori_pdb):
   #print "after_this_line_num_in_ori_pdb:"+ str(after_this_line_num_in_ori_pdb) + " in retrieve_line_matching_pair()"
   
   pdb_before_cryo_fit_in = open(pdb_before_cryo_fit)
   line_num = 0
   for line_before_cryo_fit in pdb_before_cryo_fit_in:
      line_num = line_num + 1
      if line_num > search_after_this_line_num_in_ori_pdb:
         if (line_before_cryo_fit[:4] == "ATOM") or (line_before_cryo_fit[:6] == "HETATM"):
            res = line_before_cryo_fit[17:20]
            res_num = line_before_cryo_fit[23:27]
            #print "\tres in ori:", res, "res_num in ori:",res_num
            if str(res.strip()) == str(seeking_res.strip()):
               if int(res_num) == int(seeking_res_num):
                  #print "matched"
                  #print "\tres in ori:", res, "res_num in ori:",res_num
                  chain = line_before_cryo_fit[20:22]
                  search_after_this_line_num_in_ori_pdb = line_num
                 # print "\tsearch_after_this_line_num_in_ori_pdb:", search_after_this_line_num_in_ori_pdb
                  pdb_before_cryo_fit_in.close()
                  return chain, search_after_this_line_num_in_ori_pdb
   print "match not found"
   print "\t\tseeking_res:", seeking_res, "seeking_res_num:",seeking_res_num
   print "\t\tsearch_after_this_line_num_in_ori_pdb:", search_after_this_line_num_in_ori_pdb
   exit(1)
   pdb_before_cryo_fit_in.close()
############################# end of retrieve_chain

start = time.time()
search_after_this_line_num_in_ori_pdb = 0 # initial
former_res = '' # initial
former_res_num = '' # initial
for line in pdb_after_cryo_fit_in:
   if line[:5] == "CRYST" or line[:3] == "END" or line[:6] == "HETATM" or \
      line[:5] == "MODEL" or line[:6] == "REMARK" or line[:3] == "TER" or line[:5] == "TITLE":
      f_out.write(line)
   elif line[:4] == "ATOM":
      res = line[17:20]
      res_num = line[23:27]
      new_line = ''
      if ((res == former_res) and (res_num == former_res_num)):
         new_line = line[:20] + former_retrieved_chain + line[22:]
      else:
         retrieved_chain, search_after_this_line_num_in_ori_pdb = retrieve_chain(res, res_num, search_after_this_line_num_in_ori_pdb)
         new_line = line[:20] + retrieved_chain + line[22:]
      f_out.write(new_line)
      former_res = res
      former_res_num = res_num
      former_retrieved_chain = retrieved_chain
end = time.time()
time_took = "chain recovery finished in " + str(round((end-start), 2)) + " seconds (wallclock)."
print time_took
pdb_after_cryo_fit_in.close()
f_out.close()
