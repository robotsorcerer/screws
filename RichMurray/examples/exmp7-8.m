(*
 * Example 7.8:  Planar space robot
 *
 * RMM, 6 Dec 93
 *
 *)

<<exmp7-3.m			(* get initial calculations *)

(* This term has to be simplified in an odd way to avoid strange terms *)
g3 = Simplify[Simplify[ Lie[g1, g2, q], Trig->False ]]
