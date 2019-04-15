(*
 * Example 7.6: Fingertip rolling on an object
 *
 * RMM, 6 Dec 93
 *
 *)

q = {q1, q2, q3, q4, q5};
v = {dq1, dq2, dq3, dq4, dq5};
daf = {dq1, dq2};	dao = {dq3, dq4};	dpsi = dq5;
Rpsi = {{Cos[q5], -Sin[q5]}, {-Sin[q5], -Cos[q5]}};

(* geometric parameters for the finger and object *)
Ko = {{0,0}, {0,0}};		Kf = {{1,0}, {0,1}};
Mo = {{1,0}, {0,1}};		Mf = {{1,0}, {0, Cos[q1]}};
To = {0, 0};			Tf = {0, -Tan[q1]};

(* Write down the rolling kinematics *)
{w1, w2} = Jac[Mf.daf - Rpsi.Mo.dao, v];
     w3  = Jac[Tf.Mf.daf + To.Mo.dao - dpsi, v];

(* verify the vector fields are correct *)
g1 = {0, 1/Cos[q1], -Sin[q5], -Cos[q5], -Sin[q1]/Cos[q1]};
g2 = {-1, 0, -Cos[q5], Sin[q5], 0};

z11 = Simplify[ w1.g1 ];
z12 = Simplify[ w1.g2 ];
z21 = Simplify[ w2.g1 ];
z22 = Simplify[ w2.g2 ];
z31 = Simplify[ w3.g1 ];
z32 = Simplify[ w3.g2 ];
