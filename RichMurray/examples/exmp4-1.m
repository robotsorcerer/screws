(*
 * Example 4.1:  Dynamics of a spherical pendulum
 *
 * RMM, 5 Nov 93
 *
 *)

<<Jac.m				(* define Jacobian function *)

(* Define some trig simplification rules *)
trigexpand = { Sin[x_]^2 -> (1 - Cos[x]^2) };
trigreduce = { (1 - Cos[x_]^2) -> Sin[x]^2 };
trigSimplify[expr_] := Simplify[expr /. trigexpand, Trig->False] /. trigreduce;

(* Define variables for the configuration, velocity, acceleration *)
q = {th, phi};
v = {dth, dphi};
a = {ddth, ddphi};

(* Compute the position and velocity of the mass *)
r = {l Sin[th] Cos[phi], l Sin[th] Sin[phi], -l Cos[th]}
dr = Jac[r, q] . v;

(* Calculate the Lagrangian *)
L = trigSimplify[Simplify[
      m/2 dr.dr + m g l Cos[th] /. Sin[x_]^2 -> (1-Cos[x]^2),
      Trig->False
    ]]

(* Calulate the different terms in Lagrange's equations *)
tmp1 = Simplify[ Jac[L, v] ];
tmp2 = Jac[tmp1, q] . v + Jac[tmp1, v] . a;
tmp3 = Simplify[Jac[L, q]];

(* Put together the equations of motion *)
eqs = tmp2 - tmp3
