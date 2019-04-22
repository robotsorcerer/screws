(*
 * Example 4.8:  Comparison of joint space and workspace controllers
 *
 * RMM, 4 Aug 93
 *
 * This file uses Simulate.m, a simulation program for Mathematica.
 * Simulate.m is an undocumented, unsupported package written by
 * John Hauser and Richard Murray.  Contact one of them if you want to
 * try getting it to work on your system.
 *)

<<Simulate.m			(* Mathematica simulator *)
<<ch4sec2-3.m			(* get dynamics for planar manipulator *)

(* Joint space controller *)
BuildSystem[joint,
  States[{q, dq}];
  Derivs[{th1dot, th2dot, dth1dot, dth2dot}];
  Inputs[{tau1, tau2}];
  Outputs[{x,y}];
  Params[m1->2, m2->1, I1->0, I2->0];
  Params[th1d->1, th2d->1, L1->1, L2->1, R1->1, R2->1];
  Params[kv->10,kp->100];

  {tau1, tau2} =
    Inertia . (-kv {dth1,dth2} - kp {th1-th1d, th2-th2d}) + Coriolis.dq;

  (* Finger dynamics *)
  {th1dot,th2dot} = dq;
  {dth1dot,dth2dot} = Inverse[Inertia] . ({tau1, tau2} - Coriolis . dq);

  (* Now generate the output of the end-effector *)
  x = L1 Cos[th1] + L2 Cos[th1+th2];
  y = L1 Sin[th1] + L2 Sin[th1+th2];
];

(* Workspace controller *)
BuildSystem[work,
  States[{q, dq}];
  Derivs[{th1dot, th2dot, dth1dot, dth2dot}];
  Inputs[{tau1, tau2}];
  Outputs[{x,y}];
  Local[{dx,dy,xd,yd}];
  Params[m1->2, m2->1, I1->0, I2->0, th1d->1, th2d->1];
  Params[L1->1, L2->1, R1->1, R2->1];
  Params[kv->10,kp->100];

  (* Calculate the position and speed of the end-effector *)
  x = L1 Cos[th1] + L2 Cos[th1+th2];
  y = L1 Sin[th1] + L2 Sin[th1+th2];

  {dx, dy} = J . {dth1,dth2};

  (* Figure out the desired position of the end-effector *)
  xd = L1 Cos[th1d] + L2 Cos[th1d+th2d];
  yd = L1 Sin[th1d] + L2 Sin[th1d+th2d];

  (* Computed torque control law *)
  {tau1, tau2} = 
    Inertia . Jinv . (-kv {dx,dy} - kp {x-xd, y-yd}) +
    Coriolis.dq + Inertia . Jinvdot . {dx,dy};

  (* Finger dynamics *)
  {th1dot,th2dot} = dq;
  {dth1dot,dth2dot} = Inverse[Inertia] . ({tau1, tau2} - Coriolis . dq);

];

(* Run the simulations and generate some data *)
Simu[joint, {0,2}, Initial->{0,0.5,0,0}, Output->"joint.dat", Method->lsoda]
Simu[work,  {0,2}, Initial->{0,0.5,0,0}, Output->"work.dat", Method->lsoda]
