# author: Doonam, Serdal (idea of having TER)

# essential
# for amber03.ff
# purpose 1:
## remove entire pdb lines
#### if it has P, O1P, O2P, O3P atoms in 5' nucleic acids (identified by resid 1, previous TER sign and new chain)
#### if it has O3P atoms in non-5' (e.g. middle or 3') nucleic acids

# purpose 2 (G->RG)

# purpose 3 (AH->RA) # for chimera derived ribosomes

# optional
# purpose 4 (write ATOM, HETATM, TER only)

# purpose 5 (if RESIDUE NAME is 4 characters):
### make it 3 character (for example SERA -> SER in 5t2a.pdb (cit->pdb by phenix.cif_as_pdb))

# purpose 6 (rename OP1/OP2 to O1P/O2P if MIA)

# purpose 7 (add TER between chains)

# purpose 8-1 (if one wants):
#### change MIA to adenine (both rename and actual replacement from SCH3 to H)

# purpose 8-2 (if a developer wants):
#### remove MIA


import glob, os, sys

def clean_main(input_pdb_file_name, bool_rna_name_reposition, bool_remove_MIA, bool_MIA_to_A, bool_remove_metals):
  output_pdb_file_name = leave_ATOM_END_HETATM_TER (input_pdb_file_name)
  
  if (bool_remove_metals == "True"):
   output_pdb_file_name = remove_metals (output_pdb_file_name)

  output_pdb_file_name = start_resnum_at_1_at_each_chain (output_pdb_file_name) # to deal with nucleosome pdb file
  #output_pdb_file_name = start_atom_num_at_1_at_each_chain (output_pdb_file_name) # to minimize any atom number related error
  
  output_pdb_file_name = clean_RNA_for_chimera_derived_ribosome (output_pdb_file_name)  
  output_pdb_file_name = clean_RNA_residues_for_amber03 (output_pdb_file_name)
  output_pdb_file_name = clean_RNA_atoms_for_amber03 (output_pdb_file_name)
  output_pdb_file_name = remove_the_fourth_character(output_pdb_file_name)
  
  output_pdb_file_name = no_more_100k_atom_num (output_pdb_file_name)
  
  output_pdb_file_name = put_TER_between_chains (output_pdb_file_name)
  output_pdb_file_name = clean_HETATM_7C4 (output_pdb_file_name)
  output_pdb_file_name = clean_RNA_OP1 (output_pdb_file_name, bool_remove_MIA, bool_MIA_to_A)
  
  output_pdb_file_name = remove_water (output_pdb_file_name) # gromacs cannot handle water
  
  output_pdb_file_name = remove_OXT (output_pdb_file_name) # for a pdb file like Dieter's that has OXT before "TER"
  
  final_output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_gromacs.pdb"
  cmd = "mv " + output_pdb_file_name + " " + final_output_pdb_file_name
  os.system(cmd)
  
  return final_output_pdb_file_name
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

def remove_metals(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_metal.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    if (line[17:19] == "MG") or (line[17:20] == " MG"): # MGI is possible
      # _MG is for original pdb file, MG_ is for swisspdb viewer saved pdb file
      print "MG removed"
      continue
    elif (line[17:19] == "ZN") or (line[17:20] == " ZN"):
      # _ZN is for original pdb file, ZN_ is for swisspdb viewer saved pdb file
      print "ZN removed"
      continue
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of remove_metals


def start_resnum_at_1_at_each_chain(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_resnum_starts_at_1.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  former_chain_name = ''
  subtract_this_number_from_res_num = ''
  for line in f_in:
    if (line[:3] == "END" or line[:3] == "TER"):
      f_out.write(line)
    else:
      chain = line[21:22]
      res_num = line[22:26]
      if (chain != former_chain_name): # A new chain starts
        former_chain_name = chain
        if (res_num != 1):
          subtract_this_number_from_res_num = int(res_num) - 1
      new_res_num = int(res_num) - subtract_this_number_from_res_num
      new_line = ''
      if (new_res_num < 10):
        new_line = line[:22] + "   " + str(int(res_num) - subtract_this_number_from_res_num) + line[26:]
      elif (new_res_num < 100):
        new_line = line[:22] + "  " + str(int(res_num) - subtract_this_number_from_res_num) + line[26:]
      elif (new_res_num < 1000):
        new_line = line[:22] + " " + str(int(res_num) - subtract_this_number_from_res_num) + line[26:]
      else:
        new_line = line[:22] + str(int(res_num) - subtract_this_number_from_res_num) + line[26:]
      f_out.write(new_line)
    
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of start_resnum_at_1_at_each_chain


def start_atom_num_at_1_at_each_chain(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_atomnum_starts_at_1.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  former_chain_name = ''
  new_atom_num = ''
  for line in f_in:
    if (line[:3] == "END" or line[:3] == "TER"):
      f_out.write(line)
    else:
      chain = line[21:22]
      atom_num = line[6:11]
      if (chain != former_chain_name): # A new chain starts
        former_chain_name = chain
        new_atom_num = 0
      new_atom_num = new_atom_num + 1
      new_line = ''
      if (new_atom_num < 10):
        new_line = line[:6] + "    " + str(new_atom_num) + line[11:]
      elif (new_atom_num < 100):
        new_line = line[:6] + "   " + str(new_atom_num) + line[11:]
      elif (new_atom_num < 1000):
        new_line = line[:6] + "  " + str(new_atom_num) + line[11:]
      elif (new_atom_num < 10000):
        new_line = line[:6] + " " + str(new_atom_num) + line[11:]
      else:
        new_line = line[:6] + str(new_atom_num) + line[11:]
      f_out.write(new_line)
    
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of start_atom_num_at_1_at_each_chain


def clean_RNA_for_chimera_derived_ribosome(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_AH.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  new_line = '' #essential initial value
  for line in f_in:
    residue = line[17:20]
    if (residue == "AH ") or (residue == "AI "):
      new_line = line[:17] + line[17:20].replace('AH ', ' RA') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('AI ', ' RA') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "TH ") or (residue == "TI "):
      new_line = line[:17] + line[17:20].replace('TH ', ' RT') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('TI ', ' RT') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "GH ") or (residue == "GI "):
      new_line = line[:17] + line[17:20].replace('GH ', ' RG') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('GI ', ' RG') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "CH ") or (residue == "CI "):
      new_line = line[:17] + line[17:20].replace('CH ', ' RC') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('CI ', ' RC') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "UH ") or (residue == "UI "):
      new_line = line[:17] + line[17:20].replace('UH ', ' RU') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('UI ', ' RU') + new_line[20:]
      f_out.write(new_line)
    else:
      f_out.write(line)
    
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of clean_RNA_for_chimera_derived_ribosome function


def clean_RNA_residues_for_amber03(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_amber03.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    residue = line[17:20]    
    if (residue == "  A") or (residue == "A  "): # residue == "A  " is for chimera
      new_line = line[:17] + line[17:20].replace('  A', ' RA') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('A  ', ' RA') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "  T") or (residue == "T  "):
      new_line = line[:17] + line[17:20].replace('  T', ' RT') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('T  ', ' RT') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "  G") or (residue == "G  "):
      new_line = line[:17] + line[17:20].replace('  G', ' RG') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('G  ', ' RG') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "  C") or (residue == "C  "):
      new_line = line[:17] + line[17:20].replace('  C', ' RC') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('C  ', ' RC') + new_line[20:]
      f_out.write(new_line)
    elif (residue == "  U") or (residue == "U  "):
      new_line = line[:17] + line[17:20].replace('  U', ' RU') + line[20:]
      new_line = new_line[:17] + new_line[17:20].replace('U  ', ' RU') + new_line[20:]
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of clean_RNA_residues_for_amber03 function


def clean_RNA_atoms_for_amber03(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_amber03.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    third_character_of_atom = line[15:16]
    if (third_character_of_atom == "*"):
      new_line = line[:15] + line[15:16].replace('*', '\'') + line[16:]
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of clean_RNA_atoms_for_amber03 function


def remove_the_fourth_character(input_pdb_file_name):
  # remove the fourth character in residue part, it should be either two characters like RA or \
  # three characters like ASP
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_4th_char.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    TER_candidate = line[0:3]
    if (TER_candidate == "HET"):
      fourth_character = line[20:21]
      new_line = line[:20] + " " + line[21:]
      f_out.write(new_line)
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of remove_the_fourth_character function

def put_TER_between_chains(input_pdb_file_name):
  f_in = open(input_pdb_file_name, 'r')
  output_pdb_file_name = input_pdb_file_name[:-4] + "_TER_added.pdb"
  f_out = open(output_pdb_file_name, 'w')
  old_residue_number = -1 # just initial value of course
  for line in f_in:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      residue_number = line[22:26]
      #print "\nold_residue_number: ", old_residue_number
      #print "residue_number: ", residue_number
      if (int(residue_number) < int(old_residue_number)):
        #print "a new chain starts, although chain name may not have changed"
        if old_TER_candidate != "TER":
          f_out.write("TER\n")
    f_out.write(line)
    old_residue_number = residue_number
    old_TER_candidate = line[0:3]
    #print "old_TER_candidate:", old_TER_candidate
  f_out.close()
  f_in.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of put_TER_between_chains function

def clean_HETATM_7C4(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_HETATM.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    residue = line[17:20]
    if residue == "7C4":
      print residue, " removed"
      continue
    elif residue == "BMA":
      print residue, " removed"
      continue
    elif residue == "CSX":
      print residue, " removed"
      continue
    elif residue == "GDP":
      print residue, " removed"
      continue
    elif residue == "HYP": # needed to run emd_3981
      print residue, " removed"
      continue
    elif residue == "ILX":
      print residue, " removed"
      continue
    elif residue == "NAG":
      print residue, " removed"
      continue
    elif residue == "SEP":
      print residue, " removed"
      continue
    elif residue == "TRX":
      print residue, " removed"
      continue
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of clean_HETATM_7C4 function

def remove_OXT(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_wo_HOH.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    atom = line[13:16]
    if atom == "OXT":
      print line, " removed due to OXT"
      continue
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of remove_OXT function

def remove_water(input_pdb_file_name):
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_wo_HOH.pdb"
  f_out = open(output_pdb_file_name, 'wt')
          
  for line in f_in:
    residue = line[17:20]
    if residue == "HOH":
      print residue, " removed"
      continue
    else:
      f_out.write(line)
  f_in.close()
  f_out.close()
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)
  return output_pdb_file_name
# end of remove_water function

def clean_RNA_OP1(input_pdb_file_name, bool_remove_MIA, bool_MIA_to_A):
  bool_MIA_to_A_done = 0
  bool_MIA_found = 0
  output_pdb_file_name = input_pdb_file_name[:-4] + "_cleaned_for_gromacs.pdb"
  print "\toutput_pdb_file_name: ", output_pdb_file_name
  with open(output_pdb_file_name, "w") as f_out:
    with open(input_pdb_file_name) as f_in:
      first_res_num_encountered = 0 # initial
      first_chain_encountered = 0 # initial
      new_chain_start = 0 # just initial
      consider_this_resnum_as_the_first = 0 # just initial
      deal_this_res_num_as_first = -999 # I deliberately put an impossible res_num
      for line in f_in:
        TER_candidate = line[0:3]
        atom = line[13:16]
        residue = line[17:20]
        chain = line[21:22]
        res_num = line[23:26]
        if first_chain_encountered == 0:
          first_chain_encountered = 1
          new_chain_start = 0
        else:
          if old_chain != chain:
            new_chain_start = 1
            consider_this_resnum_as_the_first = res_num
          else:
            new_chain_start = 0
        old_chain = chain
        trimmed_residue = residue.replace(" ", "")
          
        if trimmed_residue == "RA" or trimmed_residue == "RT" or trimmed_residue == "RU" \
          or trimmed_residue == "RG" or trimmed_residue == "RC" \
          or trimmed_residue == "DA" or trimmed_residue == "DT" or trimmed_residue == "DU" \
          or trimmed_residue == "DG" or trimmed_residue == "DC":
          if first_res_num_encountered == 1:
            deal_this_res_num_as_first = res_num
            first_res_num_encountered = 0
          trimmed_res_num = res_num.replace(" ", "")
          #print "trimmed_res_num:", trimmed_res_num, "."
          # clean the "first" nucleic acid
          if trimmed_res_num == "1" or deal_this_res_num_as_first == res_num or \
                    res_num == consider_this_resnum_as_the_first : 
            if (atom != "P  ") and (atom != "OP1") and (atom != "O1P") and (atom != "OP2") and (atom != "O2P") \
              and (atom != "OP3") :
              f_out.write(line)
            else:
              if res_num == "   1":
                print "removed atom (", atom, ") in the very first nucleic acid (",residue,res_num,") in chain (",chain,")"
              elif deal_this_res_num_as_first == res_num:
                print "removed atom (", atom, ") in the first nucleic acid (",residue,res_num,") in chain (",chain,") identified by previous TER sign"
              else:
                print "removed atom (", atom, ") in the first nucleic acid (",residue,res_num,") in this chain (",chain,")"

          # clean the "non-first" nucleic acid
          else:
            if atom != "OP3" :
              f_out.write(line)
            else:
              print "removed OP3 in the non-first nucleic acid"
        elif residue == "MIA": 
          bool_MIA_found = 1
          OP1_OP2_changed = 0
          if (bool_remove_MIA == 0):
            if atom == "OP1":
              OP1_OP2_changed = 1
              new_line = line[:13] + line[13:16].replace('OP1', 'O1P') + line[16:]
              f_out.write(new_line)
            if atom == "OP2":
              OP1_OP2_changed = 1
              new_line = line[:13] + line[13:16].replace('OP2', 'O2P') + line[16:]
              f_out.write(new_line)
            if (bool_MIA_to_A == 1): # not ideal, this is just temporary for development purpose
              if atom == "S10" or atom == "C11" or atom == "C12" or atom == "C13" or atom == "C14" or atom == "C15" or atom == "C16": 
                continue
              new_line = line[:17] + line[17:20].replace('MIA', 'A  ') + line[20:]
              f_out.write(new_line)
              print "changed MIA to A"
              bool_MIA_to_A_done = 1
            else:
              if (OP1_OP2_changed != 1):
                f_out.write(line)
          #else: # not ideal, this is just temporary for development purpose
            #print "MIA removed"
        else:
          f_out.write(line)
        if TER_candidate == "TER":
          first_res_num_encountered= 1
  f_out.close()
  f_in.close()
  
  cmd = "rm " + input_pdb_file_name
  os.system(cmd)  
  return output_pdb_file_name
# end of clean_RNA_OP1 function

def remove_former_files():
  for to_be_removed_file in glob.glob("*_cleaned_for_gromacs*"):
    cmd = "rm " + to_be_removed_file
    os.system(cmd)
  for to_be_removed_file in glob.glob("*_removed.pdb"):
    cmd= "rm " + to_be_removed_file
    os.system(cmd)
# end of remove_former_files function

def know_the_biggest_atom_num(input_pdb_file_name):
  # to solve truncation problem of 80S.mdfit_fitted_by_chimera_cleaned_for_gromacs_chain_2.pdb
  f_in = open(input_pdb_file_name)
  max_atom_numer = -999
  for line in f_in:
    atom_number = line[5:11]
    #print "atom_number:", atom_number
    if atom_number != " *****":
      if (atom_number > max_atom_numer):
        max_atom_numer = atom_number
  f_in.close()
  return max_atom_numer
# end of know_the_biggest_atom_num function

def no_more_100k_atom_num(input_pdb_file_name):
  # to solve truncation problem of 80S.mdfit_fitted_by_chimera_cleaned_for_gromacs_chain_2.pdb
  f_in = open(input_pdb_file_name)
  output_pdb_file_name = input_pdb_file_name[:-4] + "_no_100k_AtomNum.pdb"
  f_out = open(output_pdb_file_name, 'wt')
  for line in f_in:
    if line[0:4] == "ATOM":
      atom_number = line[5:11]
      print "atom_number:", atom_number
      
      try :
        int(atom_number)
      except :
        print "atom_number is not integer like A0000 in ala_fixed....pdb, so just write this line without reducing atom_number"
        f_out.write(line)
        continue
        
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

  print "\tInput format: python clean_pdb_for_gromacs.py input_pdb_file_name bool_rna_name_reposition bool_remove_MIA bool_MIA_to_A"
  print "\tExample usage: python clean_pdb_for_gromacs.py 80S.CS.nomag_MIA.pdb 0 0 0"

  print "\tIf a user does \"python clean_pdb_for_gromacs.py <input_pdb_file_name>\" only, then by default 0 0 0"
  args=sys.argv[1:]
  if (len(args) >= 1):
    input_pdb_file_name = args[0] # pdb input file
    if (input_pdb_file_name[-4:] != ".pdb"):
      print "Entered input_pdb_file is not in pdb format."
      exit(1)
  else:
    print "Please provide a pdb file."
    exit(1)
  
  bool_rna_name_reposition = 0 # default_value
  bool_remove_MIA = 0 # default_value
  bool_MIA_to_A = 0 # default_value
  bool_remove_metals = 1 # default_value
  if (len(args) == 2):
    bool_rna_name_reposition = args[1]
  elif (len(args) == 3):
    bool_rna_name_reposition = args[1]
    bool_remove_MIA = args[2]
  elif (len(args) == 4):
    bool_rna_name_reposition = args[1]
    bool_remove_MIA = args[2]
    bool_MIA_to_A = args[3]
  else: # len(args) = 5 # python cryo_fit code uses this else
    bool_rna_name_reposition = args[1]
    bool_remove_MIA = args[2]
    bool_MIA_to_A = args[3]
    bool_remove_metals = args[4]
  clean_main(input_pdb_file_name, bool_rna_name_reposition, bool_remove_MIA, bool_MIA_to_A, \
                               bool_remove_metals)

  ''' # for dealing many pdb files
  args=sys.argv[1:]
  if len(args)<1:
    count = 0
    for pdb_file in glob.glob("*.pdb"):
      input_pdb_file_name = pdb_file # if there is only 1 pdb file in this folder, use it
      count +=1
      if count == 2:
        print "Please specify one input PDB file"
        print "example usage: clean_pdb_for_gromacs.py input.pdb"
        sys.exit("clean_pdb_for_gromacs exits now (expecting a pdb file at next run)")
    clean(input_pdb_file_name)
  else:
    input_pdb_file_name=args[0] # pdb input file
    clean(input_pdb_file_name)
  '''
  
