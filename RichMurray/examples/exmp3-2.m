(*
 * Example 3.2:  Elbow forward kinematics
 *
 * Richard M. Murray
 * 22 January 1992
 *
 * This file contains a portion of the file elbow.m
 *)

<<RobotLinks.m

(* twist axes for elbow robot, reference frame at base *)
xi1 = RevoluteTwist[{0,0,l0}, {0,0,1}];		(* base *)
xi2 = RevoluteTwist[{0,0,l0}, {-1,0,0}];
xi3 = RevoluteTwist[{0,l1,l0}, {-1,0,0}];	(* elbow *)
xi4 = RevoluteTwist[{0,l1+l2,l0}, {0,0,1}];	(* wrist *)
xi5 = RevoluteTwist[{0,l1+l2,l0}, {-1,0,0}];
xi6 = RevoluteTwist[{0,l1+l2,l0}, {0,1,0}];

gst0 = RPToHomogeneous[IdentityMatrix[3], {0,l1+l2,l0}];

(* Forward Kinematics (plus tuned simplification) *)
gst = Simplify[
  Simplify[ TwistExp[xi1,th1] . TwistExp[xi2,th2] . TwistExp[xi3,th3] ] .
  TwistExp[xi4,th4] . TwistExp[xi5,th5] . TwistExp[xi6,th6] . gst0,
  Trig->False
];

p = RigidPosition[gst]
R = RigidOrientation[gst]
