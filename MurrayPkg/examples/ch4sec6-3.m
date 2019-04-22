(*
 * 6.3 Example: A planar manipulator moving in a slot
 *
 * RMM, 6 Aug 93
 *
 *)

<<Jac.m				(* include Jac[] function *)

(* Normal vector and end-effector location *)
n = {Sin[alpha], -Cos[alpha]}
p = {l1 Cos[th1] + l2 Cos[th1+th2], l1 Sin[th1] + l2 Sin[th1+th2]}

(* Compute the constraint *)
h = Simplify[(p - {l,0}) . n];

(* Find the gradient of the constraint *)
dh = Jac[h, {th1, th2}]

(* Store the inverse kinematics function *)
f = {
  ArcTan[Sin[alpha] s / (l + Cos[alpha] s)] +
    ArcCos[(s^2 + 2 l Cos[alpha] s + l^2 + l1^2 - l2^2)/
           (2 l1 Sqrt[s^2 + 2 l Cos[alpha] s + l^2])],
  Pi + ArcCos[(l1^2 + l2^2 - l^2 - 2 l Cos[alpha] s - s^2)/(2 l1 l2)]
}

(* Compute the Jacobian of this mapping *)
J = Jac[f, {s}]
