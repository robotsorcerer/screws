(*
 * Example: Two SCARA fingers grasping a box
 *
 * Richard M. Murray
 * 27 Jan 92
 *
 *)

<<scara.m			(* read scara kinematics *)

(* Construct the change of coordinates from contact to spatial coordinates *)
(* Only do the construction at the identity object configuration *)

Rc1s1 = {{0,1,0}, {0,0,1}, {1,0,0}};  pc1s1 = {0,b-r,a};
gc1s1 = RPToHomogeneous[Rc1s1, pc1s1];
Adgs1c1inv = RigidAdjoint[RigidInverse[gc1s1]];

Rc2s2 = {{1,0,0}, {0,0,-1}, {0,1,0}};  pc2s2 = {0,-b+r,a};
gc2s2 = RPToHomogeneous[Rc2s2, pc2s2];
Adgs2c2inv = RigidAdjoint[RigidInverse[gc2s2]];

(* Contact wrench basis *)
B = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,1}};

(* Calculate the hand Jacobian *)
Jh1 = Transpose[B] . Adgs1c1inv . Js;
Jh2 = Transpose[B] . Adgs2c2inv . Js /. {l1->l3, l2->l4, th1->th3, th2->th4};

(* Now verify the null space vectors *)
z1 = Jh1 . {0, 0, 1, 0} /. {
  l1 Sin[th1] + l2 Sin[th1+th2] -> 0,
  l1 Cos[th1] + l2 Cos[th1+th2] -> b-r
};

z2 = Jh2 . {0, 0, 1, 0} /. {
  l3 Sin[th3] -> -l4 Sin[th3+th4],
  l3 Cos[th3] + l4 Cos[th3+th4] -> r-b
}
