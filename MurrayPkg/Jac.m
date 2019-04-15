(*
 * Jac.m - jacobians and Lie derivatives
 *
 * John Hauser
 * 1989
 *
 *)

BeginPackage["Jac`"]

Jac::usage = "Jac[f,x] computes the derivative of f with respect to x.";
Lie::usage = "Lie[f,g,x] computes the Lie bracket of f and g wrt x.";
LieD::usage = "LieD[f,h,x] compute the Lie derivitive of h wrt f.";
Adj::usage = "Adj[v1,v2,x,k] calculate the kth bracket of v2 wrt v1.";

Begin["Private`"]

Jac[f_, x_] := 
        If[
                VectorQ[ f ],
                Table[ D[ f[[i]], x[[j]] ], {i, Length[f]} , {j, Length[x]} ],
                Table[ D[ f, x[[j]] ], {j, Length[x]} ]
          ]

Lie[v1_, v2_, x_] := Jac[v2,x].v1 - Jac[v1,x].v2

Adj[v1_, v2_, x_, k_] :=
        If[ k==0, v2, Lie[ v1, Adj[ v1, v2, x, k-1 ], x] ]

LieD[f_, h_, x_] := Jac[h,x].f

End[]
EndPackage[]
