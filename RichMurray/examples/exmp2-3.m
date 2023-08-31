(*
 * Example 2.3:  (no title)
 *
 * RMM, 4 Nov 93
 *
 *)

<<Screws.m

(* Problem setup *)
w = {0,0,1};
q = {0,l1,0};
xi = Join[-Skew[w]. q, w]

(* Compute the exponential of the twist *)
expxi = TwistExp[xi, th]

(* Now compute the transformation between frames *)
gab0 = RPToHomogeneous[IdentityMatrix[3], {0,l1,0}]
gab = trigSimplify[expxi . gab0]
