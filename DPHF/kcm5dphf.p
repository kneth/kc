(* New HP435 Pascal Version of Derpar and Hopf continuation
   program. Version containing both the derpar-algorithm
   and the Hopf-bifurcation calculation. Version 25/10 1994 *)


PROGRAM kcm5dphf(INPUT, OUTPUT); 

 IMPORT
  ARG;

 $include 'kcm5const.p'$
 $include 'dphf.def'$
 $include 'kcm5proc.p'$
 $include 'dphf.proc'$

 BEGIN
  detnumparam;
  IF (numparam=1) THEN sp_derpar_driver
  ELSE IF (numparam=2) THEN hf_derpar_driver;
 END.
