(*
 * Example 7.4: Disk rolling on a plane
 *
 * RMM, 6 Dec 93
 *
 *)

<<Jac.m

q = {x, y, theta, phi};

(* rolling constraints *)
w1 = {1, 0, 0, -rho Cos[theta]};
w2 = {0, 1, 0, -rho Sin[theta]};

(* verify the vector fields are correct *)
g1 = {rho Cos[theta], rho Sin[theta], 0, 1};
g2 = {0,0,1,0};

z11 = w1.g1;
z12 = w1.g2;
z21 = w2.g1;
z22 = w2.g2;
