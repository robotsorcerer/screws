(*
 * Example 3.3:  SCARA kinematics with alternate base frame
 *
 * RMM, 4 Nov 93
 *
 *)

<<RobotLinks.m

(* twist axes for elbow robot, reference frame at end-effector *)
xi1 = RevoluteTwist[{0,-l1-l2,0}, {0,0,1}];	(* base *)
xi2 = RevoluteTwist[{0,-l2,0}, {0,0,1}];	(* elbow *)
xi3 = RevoluteTwist[{0,0,0}, {0,0,1}];		(* wrist *)
xi4 = PrismaticTwist[{0,0,0}, {0,0,1}];	

gst0 = RPToHomogeneous[IdentityMatrix[3], {0,0,0}];

(* Forward Kinematics *)
gst = Simplify[
  TwistExp[xi1,th1] . TwistExp[xi2,th2] .
  TwistExp[xi3,th3] . TwistExp[xi4,th4] . gst0
];
