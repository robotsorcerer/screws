(*
 * Example 7.10: kinematic car
 *
 * RMM, 6 Dec 93
 *
 *)

<<exmp7-5.m			(* get initial computations *)

g3 = Simplify[ Lie[g1, g2, q] ];
g4 = Simplify[ Lie[g1, g3, q] ];
