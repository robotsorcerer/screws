(*
 * Example 3.1:  SCARA forward kinematics
 *
 * Richard M. Murray
 * 22 January 1992
 *
 * This file contains a portion of the file scara.m
 *)

<<RobotLinks.m

(* twist axes for SCARA robot, reference frame at base *)
xi1 = {0,0,0, 0,0,1};		(* base *)
xi2 = {l1,0,0, 0,0,1};		(* elbow *)
xi3 = {l1+l2,0,0, 0,0,1};	(* wrist revolute *)
xi4 = {0,0,1, 0,0,0};		(* wrist prismatic *)

gst0 = RPToHomogeneous[IdentityMatrix[3], {0,l1+l2,0}];

(* Calculate out the individual exponentials *)
expxi1 = TwistExp[xi1,th1];
expxi2 = TwistExp[xi2,th2];
expxi3 = TwistExp[xi3,th3];
expxi4 = TwistExp[xi4,th4];

(* Forward Kinematics *)
gst = Simplify[ expxi1 . expxi2 . expxi3 . expxi4 . gst0 ];

