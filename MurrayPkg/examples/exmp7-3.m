(*
 * Example 7.3:  Planar space robot
 *
 * RMM, 6 Dec 93
 *
 *)

<<Jac.m

(* Generalized coordinates and velocities *)
q = {psi1, psi2, theta}
v = {dpsi1, dpsi2, dtheta}

(*
 * Define the kinetic energy of the system
 *
 * This only includes the rotational contributions since we have
 * already gotten rid of the center of mass (which remains constant).
 *
 *)

(* compute the velocity of the ends of the links *)
x = r Cos[theta] + l Cos[theta+psi];
y = r Sin[theta] + l Sin[theta+psi];
vel = Simplify[
  (Jac[x, {theta,psi}] . {dtheta,dpsi})^2 +
  (Jac[y, {theta,psi}] . {dtheta,dpsi})^2
];

v1 = vel /. {psi->psi1, dpsi->dpsi1};
v2 = Simplify[ vel /. {psi->psi2-Pi, dpsi->dpsi2, r->-r} ];

(* now compute the kinetic energy *)
K = 1/2 I dtheta^2 + 1/2 m v1 + 1/2 m v2;
M = Simplify[ Jac[Jac[K, v], v] ];

(* Now compute the momentum map due to the symmetry in theta *)
momen = Simplify[ D[K, dtheta] ];
omega = Jac[momen, v];

(* Compute the null space vector fields *)
g1 = {1, 0, -omega[[1]] / omega[[3]]};
g2 = {0, 1, -omega[[2]] / omega[[3]]};

z1 = omega . g1;
z2 = omega . g2;
