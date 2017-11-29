# Sanbonmatsu  12/11/01
# reads pdb
# to call 
# push(@INC, "/Users/kys/Perl");                       # Add my library
# require('pdb_indel.pl');
# &PDB($Line1); 
# for reference, here is print statement
#0         1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456
#|     |    |    |   | |    |  |       |       |       |     |     |         |
#ATOM      0  O   ILE    43      47.695 215.020 259.022  1.00  0.00      P48
# printf ('%6s%5s%5s%4s%2s%4s%4s%8.3f%8.3f%8.3f%6.2f%6.2f%s',$at,$ano,$ana,$rna,$ch,$rno,$spa_new,$x,$y,$z,$q1,$q2,$p);

sub PDB {
#Check pdb file
# at = atomtype
# ano = atomnumber
# ana = atomname
# rna = res name
# ch  = chain
# rno = res number
# spa = spaces
# ag  = general type of atom
# atom names defined without spaces
        $at=substr($_,0,6);
        $ano=substr($_,6,5);
        $ana=substr($_,11,5);
        $rna=substr($_,16,4);
        $ch=substr($_,20,2);
        $rno=substr($_,22,5);
        $spa=substr($_,28,3);
        $x=substr($_,30,8);
        $y=substr($_,38,8);
        $z=substr($_,46,8);
        $q1=substr($_,54,6);
        $q2=substr($_,60,6);
	$p=substr($_,66,10);
}
# return true
1;
