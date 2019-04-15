(*
 * Example 2.4: Rotational motion of a one degree of freedom manipulator
 *
 * RMM, 4 Nov 93
 *
 *)

<<Screws.m

(* Define some trig simplification rules *)
trigexpand = { Sin[x_]^2 -> (1 - Cos[x]^2) };
trigreduce = { (1 - Cos[x_]^2) -> Sin[x]^2 };
trigSimplify[expr_] := Simplify[expr /. trigexpand, Trig->False] /. trigreduce;

(* Build the rotation matrix *)
R = SkewExp[{0,0,1}, th];

(* Calculate velocities *)
ws = UnSkew[ trigSimplify[ Dt[R, t] . Transpose[R] ]]
wb = UnSkew[ trigSimplify[ Transpose[R] . Dt[R, t] ]]
