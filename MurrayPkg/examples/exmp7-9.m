(*
 * Example 7.9: Disk rolling on a plane
 *
 * RMM, 6 Dec 93
 *
 *)

<<exmp7-4.m			(* get initial calculations *)

g3 = Simplify[ Lie[g1, g2, q] ];
g4 = Simplify[ Lie[g2, g3, q] ];
