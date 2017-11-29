#!/bin/awk -f

BEGIN{
   FLAG=0;
}

/\[ atoms \]/{FLAG=1}
/\[ pairs \]/{FLAG=2}
/\[ bonds \]/{FLAG=3}
/\[ exclusions \]/{FLAG=4}
/\[ angles \]/{FLAG=5}
/\[ dihedrals \]/{FLAG=6}
/\[ system \]/{FLAG=7}

{
#  if (FLAG !=1 || substr($0,1,1) = "[" || substr($0,1,1)==";" ) 
   if (FLAG !=1 || substr($0,1,1) == "[" || substr($0,1,1)==";")  
     print $0 

   else 
#     printf "%6i %10s %6i %6s %6s %6i %10.4f      %-11.3f \n", $1,$2,$3,$4,$5,$6,$7,$8
      printf "%6i %10s %6i %6s %6s %6i %10.4f      %-11.3f \n", $1,$2,$3,$4,$5,$6,0,$8
     

} 
