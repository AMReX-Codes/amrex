(* definition of tau *)

tauxxdef = tauxx[i_ + half, j_] -> 2 mu[i + half, j] dudx[i + half, j]

tauyydef = tauyy[i_, j_ + half] -> 2 mu[i, j + half] dvdy[i, j + half]

tauxydef = tauxy[i_, j_] -> mu[i, j] (dudy[i, j] + dvdx[i, j])

(* definitions of derivatives *)

dudydef1 = dudy[i_, j_ + half] -> (u[i, j + 1] - u[i, j])/hy

dvdxdef1 = dvdx[i_ + half, j_] -> (v[i + 1, j] - v[i, j])/hx

dudydef2 = dudy[i_ + half, j_] -> (u[i, j + 1] - u[i, j - 1] + u[i + 1, j + 1] - u[i + 1, j - 1])/(4 hy)

dvdxdef2 = dvdx[i_, j_ + half] -> (v[i + 1, j + 1] - v[i - 1, j + 1] + v[i + 1, j] - v[i - 1, j])/(4 hx)

dudxdef  = dudx[i_ + half, j_] -> (u[i + 1, j] - u[i, j])/hx

dvdydef  = dvdy[i_, j_ + half] -> (v[i, j + 1] - v[i, j])/hy

(* definitions used for fortran output *)

murepl1 = mu[i_, j_ + half] -> muY[i, j + 1]

murepl2 = mu[i_ + half, j_] -> muX[i + 1, j]

urepl = u[i_, j_] -> U[i, j, 1]

vrepl = v[i_, j_] -> U[i, j, 2]

(* 
dependentCellsNotCovered is a function which returns a logical
expression suitable for inclusion in fortran.  Give an expression,
exp, we wish to determine which mesh locations are accessed by the
expression.  However, we do not wish to examine all possible
locations, only those outside the grid patch region.  So we provide a
second argument, which is a boolean function taking two arguments.
The combination will give logical expressions testing the mask for
cells utilized by the expression and for which the boolean function,
logfunc[il,jl], evaluates as true. The third argument is the name of
the mask array.
*)

Clear[dependentCellsNotCovered]
dependentCellsNotCovered[exp_, logfunc_, maskfun_] := 
  Module[{cond, lexp}, cond = False; lexp = exp; 
    For[il = -2, il <= +2, il++, 
      For[jl = -2, jl <= +2, jl++, 
        For[kl = 1, kl <= 2, kl++, 
          If[logfunc[il, jl] && 
              Coefficient[Expand[exp], U[i + il, j + jl, kl]] =!= 0, 
            cond = cond || maskfun[i + il, j + jl] > 0]]]]; cond]

(*
dependentCellsCovered is more or less the converse of
dependentCellsNotCovered.  dependentCellsCovered will return true if
all the cells in the expression (properly restricted by the mask
function) ARE covered.  dependentCellsNotCovered returned true if
one or more of the cells were not covered.
*)

Clear[dependentCellsCovered]
dependentCellsCovered[exp_, logfunc_, maskfun_] := 
  Module[{cond, lexp}, cond = True; lexp = exp; 
    For[il = -2, il <= +2, il++, 
      For[jl = -2, jl <= +2, jl++, 
        For[kl = 1, kl <= 2, kl++, 
          If[logfunc[il, jl] && 
              Coefficient[Expand[exp], U[i + il, j + jl, kl]] =!= 0, 
            cond = cond && maskfun[i + il, j + jl] == 0]]]]; cond]
(*
this is an alternate definition which treats cells with mask .eq. 1 as
untrusted, rather than mask .gt. 0
*)

dependentCellsNotCoveredalt[exp_, logfunc_, maskfun_] := 
  Module[{cond, lexp}, cond = False; lexp = exp; 
    For[il = -2, il <= +2, il++, 
      For[jl = -2, jl <= +2, jl++, 
        For[kl = 1, kl <= 2, kl++, 
          If[logfunc[il, jl] && 
              Coefficient[Expand[exp], U[i + il, j + jl, kl]] =!= 0, 
            cond = cond || maskfun[i + il, j + jl] == 1]]]]; cond]

(* 
definitions for Do One-sided Derivative in X direction.  if sign is
positive, it means extend the stencil in the positivie x direction.
if negative, extend in other direction
*)

DOneX[u_, i_, j_, sign_] := (-u[i + 2, j] + 4*u[i + 1, j] - 3*u[i, j])/(2*hx) /; sign == 1

DOneX[u_, i_, j_, sign_] := (u[i - 2, j] - 4*u[i - 1, j] + 3*u[i, j])/(2*hx) /; sign == -1

DOneX[u_, i_, j_, k_, sign_] := (-u[i + 2, j, k] + 4*u[i + 1, j, k] - 3*u[i, j, k])/(2*hx) /; sign == 1

DOneX[u_, i_, j_, k_, sign_] := (u[i - 2, j, k] - 4*u[i - 1, j, k] + 3*u[i, j, k])/(2*hx) /; sign == -1

(*
definitions for Do One-sided Derivative in Y direction.  if sign is
positive, it means extend the stencil in the positivie y direction.
if negative, extend in other direction
*)

DOneY[u_, i_, j_, sign_] := (-u[i, j + 2] + 4*u[i, j + 1] - 3*u[i, j])/(2*hy) /; sign == 1
DOneY[u_, i_, j_, sign_] := (u[i, j - 2] - 4*u[i, j - 1] + 3*u[i, j])/(2*hy) /;  sign == -1
DOneY[u_, i_, j_, k_, sign_] := (-u[i, j + 2, k] + 4*u[i, j + 1, k] - 3*u[i, j, k])/(2*hy) /; sign == 1
DOneY[u_, i_, j_, k_, sign_] := (u[i, j - 2, k] - 4*u[i, j - 1, k] + 3*u[i, j, k])/(2*hy) /; sign == -1

(* definitions for two sided deriv in x direction. *)

DTwoX[u_, i_, j_] := (u[i + 1, j] - u[i - 1, j])/(2 hx)
DTwoX[u_, i_, j_, k_] := (u[i + 1, j, k] - u[i - 1, j, k])/(2 hx)
DTwoY[u_, i_, j_] := (u[i, j + 1] - u[i, j - 1])/(2 hy)
DTwoY[u_, i_, j_, k_] := (u[i, j + 1, k] - u[i, j - 1, k])/(2 hy)

(*
DeleteFile["dog.mf"];
CopyFile["DV_2D.mF", "dog.mf"];
Splice["dog.mf"];
DeleteFile["DV_2D.F"];
CopyFile["dog.f", "DV_2D.F"];
SetFileDate["DV_2D.F"];
*)
