(*
 * Example 6.1:  Dynamics of an idealized planar pendulum
 *
 * RMM, 1 Dec 93
 *
 *)

<<Jac.m

(* Description of the constraints *)
h = (x^2 + y^2 - l^2)/2;
A = Jac[{h}, {x,y}];

(* Dynamics *)
M = {{m, 0}, {0, m}};
G = {0, m g};

(* Calculate the Lagrange multiplier *)
lambda = Simplify[
  Inverse[A . Inverse[M] . Transpose[A]] . (
    A . Inverse[M] . (-G) - (Jac[A.{xdot,ydot}, {x,y}] . {xdot,ydot})
  ) /. {x^2 -> l^2 - y^2}
]

(* Figure out the tension in the rod *)
norm[x_] := Sqrt[x . x];
tension = Simplify[ norm[{{x}, {y}} . lambda] /. {x^2 -> l^2 - y^2} ];
