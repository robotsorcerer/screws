(*
 * Example 3.9:  Jacobian for the Stanford arm
 *
 * RMM, 4 Nov 93
 *
 *)

<<RobotLinks.m

(* twist axes for SCARA robot, reference frame at base *)
xi1 = RevoluteTwist[{0,0,0}, {0,0,1}];		(* base *)
xi2 = RevoluteTwist[{0,0,l0}, {-1,0,0}];
xi3 = PrismaticTwist[{0,0,l0}, {0,1,0}];	(* shoulder *)
xi4 = RevoluteTwist[{0,l1,l0}, {0,0,1}];	(* wrist *)
xi5 = RevoluteTwist[{0,l1,l0}, {-1,0,0}];
xi6 = RevoluteTwist[{0,l1,l0}, {0,1,0}];

gst0 = RPToHomogeneous[IdentityMatrix[3], {0,l1,0}];

(* Forward Kinematics (not really needed) *)
gst = Simplify[
  TwistExp[xi1,th1] . TwistExp[xi2,th2] .
  TwistExp[xi3,th3] . TwistExp[xi4,th4] .
  TwistExp[xi5,th5] . TwistExp[xi6,th6] . g0
];

(*
 * Calculate the Jacobian "by inspection"
 *
 *)

(* First joint (not very exciting) *)
w1 = {0,0,1};
q1 = {0,0,l0};
xi1 = Join[-Skew[w1] . q1, w1]

(* Second joint *)
w2 = {-1,0,0};
w2p = SkewExp[w1,th1] . w2;
xi2p = Join[-Skew[w2p] . q1, w2p]

(* Third (prismatic) joint *)
v3 = {0,1,0};
v3p = SkewExp[w1,th1] . SkewExp[w2,th2] . v3;
xi3p = Join[v3p, {0,0,0}]

(* Figure out location of the wrist *)
qwp = {0,0,l0} + SkewExp[w1,th1] . SkewExp[w2,th2] . {0, l1+th3,0}

(* Now figure out wrist axes *)
w4 = {0,0,1};
w4p = SkewExp[w1,th1] . SkewExp[w2,th2] . w4

w5 = {-1,0,0};
w5p = SkewExp[w1,th1] . SkewExp[w2,th2] . SkewExp[w4,th4] . w5

w6 = {0,1,0};	
w5p = Simplify[
  SkewExp[w1,th1] . SkewExp[w2,th2] . SkewExp[w4,th4] .
  SkewExp[w5,th5] . w6,
  Trig->False
]
