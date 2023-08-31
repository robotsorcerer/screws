(*
 * 2.3 Example: dynamics of a two-link planar robot
 *
 * Richard M. Murray
 * 8 July 91
 *
 * This file uses Jac.m, a short package written by John Hauser
 # for computing Jacobians.
 *)

<<RobotLinks.m
<<Jac.m				(* package for evaluating derivatives *)
Needs["Algebra`Trigonometry`"];

(* Generalized coordinates for the manipulator *)
q = {th1, th2}
dq = {dth1, dth2}

(* Forward kinematics for centers of mass *)
x1b = R1 Cos[th1];			dx1b = Jac[x1b, q] . dq;
y1b = R1 Sin[th1];			dy1b = Jac[y1b, q] . dq;

x2b = L1 Cos[th1] + R2 Cos[th1 + th2];	dx2b = Jac[x2b, q] . dq;
y2b = L1 Sin[th1] + R2 Sin[th1 + th2];	dy2b = Jac[y2b, q] . dq;

xend = L1 Cos[th1] + L2 Cos[th1 + th2];
yend = L1 Sin[th1] + L2 Sin[th1 + th2];

(* Compute the jacobian of the end-effector position map *)
J = Simplify[Jac[{xend, yend}, q], Trig->False];
Jinv = Simplify[Inverse[J]];
Jinvdot = dq . Jac[Jinv, q]

(* Construct the lagrangian *)
L = Together[Expand[
    1/2 m1 (dx1b^2 + dy1b^2) + 1/2 I1 dth1^2 +
    1/2 m2 (dx2b^2 + dy2b^2) + 1/2 I2 (dth1 + dth2)^2
]];

(* Now figure out the equations of motion *)
Inertia = Expand[Together[ Jac[Jac[L, dq], dq] ], Trig->True];
Coriolis = InertiaToCoriolis[Inertia, q, dq];

(* Now figure out the coefficients for the entries in the inertia matrix *)
alpha = Inertia[[1,1]] /. th2 -> Pi/2;
beta = Together[(Inertia[[1,1]] - alpha) /. th2 -> 0] / 2;
delta = Inertia[[2,2]];
