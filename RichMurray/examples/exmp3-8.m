(*
 * Example 3.8:  Jacobian for a SCARA robot
 *
 * Richard M. Murray
 * 22 January 1992
 *
 * This file contains a portion of the file scara.m
 *)

<<exmp3-1.m			(* forward kinematics for SCARA *)

(* Spatial Jacobian *)
Js = Simplify[StackCols[
  xi1,
  RigidAdjoint[expxi1] . xi2,
  Simplify[ RigidAdjoint[expxi1 . expxi2] . xi3 ],
  Simplify[ RigidAdjoint[expxi1 . expxi2 . expxi3] . xi4 ]
]];
