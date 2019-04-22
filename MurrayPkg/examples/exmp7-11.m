(*
 * Example 7.6: Fingertip rolling on an object
 *
 * RMM, 6 Dec 93
 *
 *)

<<exmp7-6.m			(* initial calculations *)

g3 = Simplify[ Lie[g1, g2, q] ];
g4 = Simplify[ Lie[g1, g3, q] ];
g5 = Simplify[ Lie[g2, g3, q], Trig->False ];
