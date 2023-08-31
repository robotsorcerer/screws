(*
 * Example 7.2:  Hopping robot in flight
 *
 * RMM, 6 Dec 93
 *
 * Kinematics for a hopping robot in free flight.
 *)

<<Jac.m

(* Generalized coordinates and velocities *)
q = {psi, l, theta}
v = {dpsi, dl, dtheta}

(* Define the kinetic energy of the system *)
K = 1/2 I dtheta^2 + 1/2 m (l+d)^2 (dtheta+dpsi)^2 + 1/2 m dl^2

(* Now compute the momentum map due to the symmetry in theta *)
momen = D[K, dtheta];
omega = Jac[momen, v];

(* Verify the the give vector fields are in the null space of omega *)
g1 = {1, 0, -m(l+d)^2 / (I + m(l+d)^2)};
g2 = {0, 1, 0};

z1 = omega . g1;
z2 = omega . g2;
