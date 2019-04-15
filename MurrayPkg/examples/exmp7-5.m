(*
 * Example 7.5: kinematic car
 *
 * RMM, 6 Dec 93
 *
 *)

q = {x, y, theta, phi};

(* rolling constraints *)
w1 = {Sin[theta+phi], -Cos[theta+phi], 0, -l Cos[phi]};
w2 = {Sin[theta], -Cos[theta], 0, 0};

(* verify the vector fields are correct *)
g1 = {Cos[theta], Sin[theta], 1/l Sin[phi]/Cos[phi], 0};
g2 = {0,0,0,1};

z11 = Simplify[ w1.g1 ];
z12 = w1.g2;
z21 = w2.g1;
z22 = w2.g2;
