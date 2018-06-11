#author:doonam

import os, subprocess, sys

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

def count_TER_before_this(pdb_before_cryo_fit, count_TER_until_this_line):
   pdb_before_cryo_fit_in = open(pdb_before_cryo_fit)
   line_num = 0
   TER_number = 0
   for line_before_cryo_fit in pdb_before_cryo_fit_in:
      line_num = line_num + 1
      #print "line_num:",line_num, " in count_TER_before_this()"
      #print "count_TER_until_this_line:",count_TER_until_this_line
      if int(line_num) < int(count_TER_until_this_line):
         if line_before_cryo_fit[:3] == "TER":
            TER_number = TER_number + 1
   pdb_before_cryo_fit_in.close()
   return TER_number
#end of count_TER_before_this(pdb_before_cryo_fit, count_TER_until_this_line)

def retrieve_1st_line_matching_pair(pdb_before_cryo_fit, seeking_res, seeking_res_num):
   pdb_before_cryo_fit_in = open(pdb_before_cryo_fit)
   line_num = 0
   for line_before_cryo_fit in pdb_before_cryo_fit_in:
      line_num = line_num + 1
      if (line_before_cryo_fit[:4] == "ATOM") or (line_before_cryo_fit[:6] == "HETATM"):
         res = line_before_cryo_fit[17:20]
         res_num = line_before_cryo_fit[23:27]
         #print "res:",res," res_num:",res_num, " in retrieve_1st_line_matching_pair()"
         if str(res.strip()) == str(seeking_res.strip()):
            if int(res_num) == int(seeking_res_num):
               first_line_num_matching_pair = line_num
               pdb_before_cryo_fit_in.close()
               return first_line_num_matching_pair
# end of retrieve_1st_line_matching_pair(pdb_before_cryo_fit, seeking_res, seeking_res_num):

def retrieve_line_matching_pair(pdb_before_cryo_fit, seeking_res, seeking_res_num, after_this_line_num_in_ori_pdb):
   #print "after_this_line_num_in_ori_pdb:"+ str(after_this_line_num_in_ori_pdb) + " in retrieve_line_matching_pair()"
   pdb_before_cryo_fit_in = open(pdb_before_cryo_fit)
   line_num = 0
   for line_before_cryo_fit in pdb_before_cryo_fit_in:
      line_num = line_num + 1
      if line_num > after_this_line_num_in_ori_pdb:
         if (line_before_cryo_fit[:4] == "ATOM") or (line_before_cryo_fit[:6] == "HETATM"):
            res = line_before_cryo_fit[17:20]
            res_num = line_before_cryo_fit[23:27]
            if str(res.strip()) == str(seeking_res.strip()):
               if int(res_num) == int(seeking_res_num):
                  first_line_num_matching_pair = line_num
                  pdb_before_cryo_fit_in.close()
                  return first_line_num_matching_pair
# end of retrieve_line_matching_pair(pdb_before_cryo_fit, seeking_res, seeking_res_num):

def retrieve_chain(seeking_res, seeking_res_num, TER_number):
   pdb_before_cryo_fit_in = open(pdb_before_cryo_fit)
   TER_counted = 0
   line_num = 0
   for line_before_cryo_fit in pdb_before_cryo_fit_in:
      line_num = line_num + 1
      TER_can = line_before_cryo_fit[:3]
      if TER_can == "TER":
         TER_counted = TER_counted + 1
         pass
      #if line_before_cryo_fit[:4] == "ATOM" # not works for chimera fitted pdb
      if (line_before_cryo_fit[:4] == "ATOM") or (line_before_cryo_fit[:6] == "HETATM"):
         chain = line_before_cryo_fit[20:22]
         res_num = line_before_cryo_fit[23:27]
         if str(res) == str(seeking_res):
            if int(res_num) == int(seeking_res_num):
               if int(TER_counted) == int(TER_number):
                  pdb_before_cryo_fit_in.close()
                  after_this_line_num_in_ori_pdb = line_num
                  return chain, after_this_line_num_in_ori_pdb #retrieved_chain
   pdb_before_cryo_fit_in.close()
   after_this_line_num_in_ori_pdb = line_num
   return "not_retrieved", after_this_line_num_in_ori_pdb
# end of retrieve_chain(seeking_res, seeking_res_num, TER_number):

line_num = 0
first_chain_retrieved = False
for line in pdb_after_cryo_fit_in:
   line_num = line_num + 1
   if line[:5] == "CRYST" or line[:3] == "END" or line[:6] == "HETATM" or \
      line[:5] == "MODEL" or line[:6] == "REMARK" or line[:3] == "TER" or line[:5] == "TITLE":
      f_out.write(line)
   elif line[:4] == "ATOM":
      res = line[17:20]
      res_num = line[23:27]
      line_num_matching_pair = ''
      #print pdb_after_cryo_fit + "'s " + " res:",res," res_num:",res_num, " in main function"
      if first_chain_retrieved == False:
         line_num_matching_pair = retrieve_1st_line_matching_pair(pdb_before_cryo_fit, res, res_num)
      else:
         line_num_matching_pair = retrieve_line_matching_pair(pdb_before_cryo_fit, res, res_num, after_this_line_num_in_ori_pdb)
      #print "line_num_matching_pair:",line_num_matching_pair, " in main function"
      TER_number = count_TER_before_this(pdb_before_cryo_fit, line_num_matching_pair)
      retrieved_chain, after_this_line_num_in_ori_pdb = retrieve_chain(res, res_num, TER_number)
      first_chain_retrieved = True
      new_line = ''
      if retrieved_chain == "not_retrieved":
         new_line = line[:20] + "  " + line[22:]
      else:
         new_line = line[:20] + retrieved_chain + line[22:]
      f_out.write(new_line)
      
pdb_after_cryo_fit_in.close()
f_out.close()