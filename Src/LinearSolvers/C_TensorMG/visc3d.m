Off[General::spell,General::spell1];

SetOptions[$Output, PageWidth->73];
Get["Format.m"]

(*
definition of tau
*)

tauxxdef = tauxx[i_ + half, j_ ,k_] -> 
		2*mu[i+half,j,k]*dudx[i+half,j,k]
tauyydef = tauyy[i_, j_ + half,k_] ->
		2*mu[i,j+half,k] * dvdy[i,j+half,k]
tauzzdef = tauzz[i_,j_,k_+half] ->
		2*mu[i,j,k+half]*dwdz[i,j,k+half]
tauxydef = tauxy[i_,j_,k_] -> mu[i,j,k]*(dudy[i,j,k]+dvdx[i,j,k])
tauxzdef = tauxz[i_,j_,k_] -> mu[i,j,k]*(dudz[i,j,k]+dwdx[i,j,k])
tauyzdef = tauyz[i_,j_,k_] -> mu[i,j,k]*(dvdz[i,j,k]+dwdy[i,j,k])
(*
definitions of derivatives
*)
(*
diagonal derivatives
*)
dudxdef = dudx[i_+half,j_,k_] -> 
		(u[i+1,j,k]-u[i,j,k])/hx
dvdydef = dvdy[i_,j_+half,k_] ->
		(v[i,j+1,k]-v[i,j,k])/hy
dwdzdef = dwdz[i_,j_,k_+half] ->
		(w[i,j,k+1]-w[i,j,k])/hz
(*
dudy
*)
dudydef1 = dudy[i_,j_+half,k_] -> (u[i,j+1,k]-u[i,j,k])/hy
dudydef2 = dudy[i_+half,j_,k_] -> 
		(u[i,j+1,k]-u[i,j-1,k]+u[i+1,j+1,k]-u[i+1,j-1,k])/(4*hy)
(*
dudz
*)
dudzdef1 = dudz[i_,j_,k_+half]->(u[i,j,k+1]-u[i,j,k])/hz
dudzdef2 = dudz[i_+half,j_,k_] ->
		(u[i,j,k+1]-u[i,j,k-1]+u[i+1,j,k+1]-u[i+1,j,k-1])/ (4*hz)
(*
dvdx
*)
dvdxdef1 = dvdx[i_+half,j_,k_] -> (v[i+1,j,k]-v[i,j,k])/hx
dvdxdef2 = dvdx[i_,j_+half,k_] ->
		(v[i+1,j+1,k]-v[i-1,j+1,k]+v[i+1,j,k]-v[i-1,j,k])/(4*hx)
(*
dvdz
*)
dvdzdef1 = dvdz[i_,j_,k_+half]->(v[i,j,k+1]-v[i,j,k])/hz
dvdzdef2 = dvdz[i_,j_+half,k_]->
		(v[i,j,k+1]-v[i,j,k-1]+v[i,j+1,k+1]-v[i,j+1,k-1])/(4*hz)
(*
dwdx
*)
dwdxdef1 = dwdx[i_+half,j_,k_]->(w[i+1,j,k]-w[i,j,k])/hx
dwdxdef2 = dwdx[i_,j_,k_+half]->
		(w[i+1,j,k]-w[i-1,j,k]+w[i+1,j,k+1]-w[i-1,j,k+1])/(4*hx)
(*
dwdy
*)
dwdydef1 = dwdy[i_,j_+half,k_] ->
		(w[i,j+1,k]-w[i,j,k])/hy
dwdydef2 = dwdy[i_,j_,k_+half] ->
	(w[i,j+1,k]-w[i,j-1,k]+w[i,j+1,k+1]-w[i,j-1,k+1])/(4*hy)
(*
definitions used for fortran output
*)
murepl1 = mu[i_,j_+half,k_] -> muY[i,j+1,k]
murepl2 = mu[i_+half,j_,k_] -> muX[i+1,j,k]
murepl3 = mu[i_,j_,k_+half] -> muZ[i,j,k+1]
urepl = u[i_,j_,k_] -> U[i,j,k,1]
vrepl = v[i_,j_,k_] -> U[i,j,k,2]
wrepl = w[i_,j_,k_] -> U[i,j,k,3]
(*
dependentCellsNotCovered is a function which returns a logical expression 
suitable for inclusion in fortran.  Give an expression, exp, we wish to 
determine which mesh locations are accessed by the expression.  However, we 
do not wish to examine all possible locations, only those outside the grid 
patch region.  So we provide a second argument, which is a boolean function 
taking two arguments.  The combination will give logical expressions testing 
the mask for cells utilized by the
expression and for which the boolean function, logfunc[il,jl], evaluates as 
true. The third argument is the name of the mask array
*)

Clear[ dependentCellsNotCovered ]
dependentCellsNotCovered[exp_ , logfunc_ ,maskfun_] :=
  Module[{cond,lexp,il,jl,kl,ml},
	cond = False;
	lexp = exp;
	For[il=-1,il<=+1,il++,
      For[jl=-1,jl<=+1,jl++,
        For[kl=-1,kl<=+1,kl++,
          For[ml=1,ml<=3,ml++,
            If[ (logfunc[il,jl,kl]) &&
	          (Coefficient[
		        Expand[ 
		          exp
		        ] ,
		        U[i+il,j+jl,k+kl,ml]
	          ] =!= 0), cond = cond || (maskfun[i+il,j+jl,k+kl]>0)
	        ]
	      ]
        ]
      ]
    ];
    cond
  ]

(*
dependentCellsCovered is the logical inverse of dependentCellsNotCovered
*)

Clear[ dependentCellsCovered ]
dependentCellsCovered[exp_ , logfunc_ ,maskfun_] :=
  Module[{cond,lexp,il,jl,kl,ml},
	cond = True;
	lexp = exp;
	For[il=-1,il<=+1,il++,
      For[jl=-1,jl<=+1,jl++,
        For[kl=-1,kl<=+1,kl++,
          For[ml=1,ml<=3,ml++,
            If[ (logfunc[il,jl,kl]) &&
	          (Coefficient[
		        Expand[ 
		          exp
		        ] ,
		        U[i+il,j+jl,k+kl,ml]
	          ] =!= 0), cond = cond && (maskfun[i+il,j+jl,k+kl]==0)
	        ]
	      ]
        ]
      ]
    ];
    cond
  ]
(*
definitions for two sided derivs
*)

DTwoX[u_,i_,j_,k_] := (u[i+1,j,k]-u[i-1,j,k])/(2*hx)
DTwoX[u_,i_,j_,k_,n_] := (u[i+1,j,k,n]-u[i-1,j,k,n])/(2*hx)
DTwoY[u_,i_,j_,k_] := (u[i,j+1,k]-u[i,j-1,k])/(2*hy)
DTwoY[u_,i_,j_,k_,n_] := (u[i,j+1,k,n]-u[i,j-1,k,n])/(2*hy)
DTwoZ[u_,i_,j_,k_] := (u[i,j,k+1]-u[i,j,k-1])/(2*hz)
DTwoZ[u_,i_,j_,k_,n_] := (u[i,j,k+1,n]-u[i,j,k-1,n])/(2*hz)

(*
definitions for Do One-sided Derivative in X direction.  if sign is positive, 
it means extend the stencil in the positivie x direction.  if negative, 
extend in other direction
*)

Clear[DOneX]
DOneX[u_, i_, j_, k_, sign_] := 
		(-u[i + 2, j, k] + 4*u[i + 1, j, k] - 3*u[i, j, k])/
									(2*hx)  /; sign == 1
DOneX[u_, i_, j_, k_, sign_] := 
		(u[i - 2, j, k] - 4*u[i - 1, j, k] + 3*u[i, j, k])/
									(2*hx) /; sign == -1
DOneX[u_, i_, j_, k_, n_, sign_] := 
		(-u[i + 2, j, k, n] + 4*u[i + 1, j, k, n] - 3*u[i, j, k, n])/
									(2*hx) /; sign == 1
DOneX[u_, i_, j_, k_, n_, sign_] := 
		(u[i - 2, j, k, n] - 4*u[i - 1, j, k, n] + 3*u[i, j, k, n])/
									(2*hx) /; sign == -1
(*
definitions for Do One-sided Derivative in Y direction.  if sign is positive, 
it means extend the stencil in the positivie y direction.  if negative, 
extend
in other direction
*)

Clear[DOneY]
DOneY[u_, i_, j_, k_, sign_] := 
		(-u[i, j + 2, k] + 4*u[i, j + 1, k] - 3*u[i, j, k])/
									(2*hy) /; sign == 1
DOneY[u_, i_, j_, k_, sign_] := 
		(u[i, j - 2, k] - 4*u[i, j - 1, k] + 3*u[i, j, k])/
									(2*hy) /; sign == -1
DOneY[u_, i_, j_, k_, n_, sign_] := 
		(-u[i, j + 2, k, n] + 4*u[i, j + 1, k, n] - 3*u[i, j, k, n])/
									(2*hy) /; sign == 1
DOneY[u_, i_, j_, k_, n_, sign_] := 
		(u[i, j - 2, k, n] - 4*u[i, j - 1, k, n] + 3*u[i, j, k, n])/
									(2*hy) /; sign == -1
(*
definitions for Do One-sided Derivative in Z direction.  if sign is positive, 
it means extend the stencil in the positivie z direction.  if negative, 
extend
in other direction
*)
Clear[DOneZ]
DOneZ[u_, i_, j_, k_, sign_] := 
		(-u[i, j, k + 2] + 4*u[i, j, k + 1] - 3*u[i, j, k])/
									(2*hz) /; sign == 1
DOneZ[u_, i_, j_, k_, sign_] := 
		(u[i, j, k - 2] - 4*u[i, j, k - 1] + 3*u[i, j, k])/
									(2*hz)  /; sign == -1
DOneZ[u_, i_, j_, k_, n_, sign_] := 
		(-u[i, j, k + 2, n] + 4*u[i, j, k + 1, n] - 3*u[i, j, k, n])/
									(2*hz) /; sign == 1
DOneZ[u_, i_, j_, k_, n_, sign_] := 
		(u[i, j, k - 2, n] - 4*u[i, j, k - 1, n] + 3*u[i, j, k, n])/
									(2*hz)  /; sign == -1
(*
useful one-sided derivatives
*)
Clear[dvdxalt, dudyalt, dvdzalt, dwdyalt, dudzalt, dwdxalt]
dvdxalt[i_, j_ + half, k_, 
    sign_] := (DOneX[v, i, j, k, sign] + DOneX[v, i, j + 1, k, sign])/2
dudyalt[i_+half,j_,k_,sign_] := (DOneY[u,i  ,j,k,sign]+
							     DOneY[u,i+1,j,k,sign])/2
dvdzalt[i_,j_+half,k_,sign_] := (DOneZ[v,i,j  ,k,sign]+
							     DOneZ[v,i,j+1,k,sign])/2
dwdyalt[i_,j_,k_+half,sign_] := (DOneY[w,i,j,k  ,sign]+
								 DOneY[w,i,j,k+1,sign])/2
dudzalt[i_+half,j_,k_,sign_] := (DOneZ[u,i  ,j,k,sign]+
								 DOneZ[u,i+1,j,k,sign])/2
dwdxalt[i_,j_,k_+half,sign_] := (DOneX[w,i,j,k  ,sign]+
								 DOneX[w,i,j,k+1,sign])/2
(* setup to use the Format.m package  *)

(*

SetOptions[$Output,PageWidth->73];

*)

(*
substitution for all derivatives and variables
*)

allDerivAllUV = {dudxdef,
				 dvdydef,
				 dwdzdef,
				 dudydef1,dudydef2,
				 dudzdef1,dudzdef2,
				 dvdxdef1,dvdxdef2,
				 dvdzdef1,dvdzdef2,
				 dwdxdef1,dwdxdef2,
				 dwdydef1,dwdydef2,
				 urepl,
				 vrepl,
				 wrepl};
allUVW = {       urepl,
				 vrepl,
				 wrepl};
(*
transverse u derivs
*)

tduext[half+i,j,k] = trandere[i+1,j,k,1];
tduext[half+i-1,j,k] = tranderw[i-1,j,k,1];
tduext[i,j+half,k] = trandern[i,j+1,k,1];
tduext[i,j-1+half,k] = tranders[i,j-1,k,1];
tduext[i,j,k+half] = trandert[i,j,k+1,1];
tduext[i,j,k-1+half] = tranderb[i,j,k-1,1];
(*
transverse v derivs
*)

tdvext[half+i,j,k] = trandere[i+1,j,k,2];
tdvext[half+i-1,j,k] = tranderw[i-1,j,k,2];
tdvext[i,j+half,k] = trandern[i,j+1,k,2];
tdvext[i,j-1+half,k] = tranders[i,j-1,k,2];
tdvext[i,j,k+half] = trandert[i,j,k+1,2];
tdvext[i,j,k-1+half] = tranderb[i,j,k-1,2];
(*
transverse w derivs
*)
tdwext[half+i,j,k] = trandere[i+1,j,k,3];
tdwext[half+i-1,j,k] = tranderw[i-1,j,k,3];
tdwext[i,j+half,k] = trandern[i,j+1,k,3];
tdwext[i,j-1+half,k] = tranders[i,j-1,k,3];
tdwext[i,j,k+half] = trandert[i,j,k+1,3];
tdwext[i,j,k-1+half] = tranderb[i,j,k-1,3];

(*
an alternate approach which seeks to automatically determine which 
direction to use for one sided deriv
*)

altgen[lhs_,indx_,indy_,indz_,
       exp_,expalt_,expext_,
       varindx_,derivindx_,
       indexcond_,mask_] :=
	Block[
		{tmpcond,tmp,tmpalt,
		 depplus,depminus,sign,
		 line1,line2,line3,line4,line5,
		 icase,jcase,kcase},
		(* conditions are False if expression is safe to use *)
		tmpcond = dependentCellsNotCovered[
					exp[indx,indy,indz] //. allDerivAllUV,
					indexcond,mask];
		depplus = dependentCellsNotCovered[
					expalt[indx,indy,indz,+1] //.allDerivAllUV,
					indexcond,mask];
		depminus = dependentCellsNotCovered[
					expalt[indx,indy,indz,-1] //.allDerivAllUV,
					indexcond,mask];
		(* temporary *)
		If[ depplus === False , 
			sign = 1,
			If[ depminus === False ,
				sign = -1,
				sign = 0  (* means neither of one 
				             sided is safe *)
			]
		];
		(* treat 3 different cases *)
		Which[
			tmpcond === False,
			(* exp does not extend into masked region *)
			FortranAssign[
				lhs,
				exp[indx,indy,indz] //.allDerivAllUV,
				AssignToArray->{U},
				AssignPrecision->Infinity
			],
			
			tmpcond =!= False && sign != 0,
			(* exp extends outside, output conditional mask *)

			tmp = FortranAssign[
				tmpcond,
				AssignToArray->{mask},
				AssignIndent->"",
				AssignPrecision->Infinity
			];
			tmpalt =dependentCellsNotCovered[
					expalt[indx,indx,indz,sign]//.allDerivAllUV,
						indexcond,mask];
			line1 = StringForm["      if(``) then ", tmp];
			line2 = FortranAssign[
				lhs,
				expalt[indx,indy,indz,sign]//.allDerivAllUV,
				AssignToArray->{U},
				AssignPrecision->Infinity
			];
			line3 = StringForm["      else"];
			line4 = FortranAssign[
				lhs,
				exp[indx,indy,indz] //.allDerivAllUV,
				AssignToArray->{U},
				AssignPrecision->Infinity
			];
			If[tmpalt =!= False ,
				(* this is an error, the alternate form should
				   be specified to be inside safe region *)
		   		line5 = StringForm[" error in tandergen"]
		   	,
		   		line5 = StringForm["      endif"]
			];
			ColumnForm[{line1,line2,line3,line4,line5}],
			
			True,
			(* cannot use exp or expalt, must be
			   externally supplied derivative *)

			(* need to figure out where to evaluate
			   the transverse derivative *)
			icase = Which[ indx === i-1+half, -1,
						   indx === i       ,  0,
						   indx === i  +half, +1,
						   _, Print["error in icase"]];
			jcase = Which[ indy === j-1+half, -1,
						   indy === j       ,  0,
						   indy === j  +half, +1,
						   _, Print["error in jcase"]];

			kcase = Which[ indz === k-1+half, -1,
						   indz === k       ,  0,
						   indz === k  +half, +1,
						   _, Print["error in kcase"]];
		
			FortranAssign[
					lhs, expext[i+icase,
								j+jcase,
								k+kcase,
								varindx,derivindx],
					AssignToArray->{trandere,tranderw,
									trandern,tranders,
									trandert,tranderb},
					AssignPrecision->Infinity]	
		]	
	]
(*
a short-hand function
*)

FA[x_] := FortranAssign[x, AssignToArray->{U,muX,muY,muZ,a,u, maskn,maske,maskw, masks,maskt,maskb},
						AssignIndent->"", AssignPrecision->Infinity];		

(*
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D1.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D1.F"];
CopyFile[ "dog.f", "DV_3D1.F" ];
SetFileDate["DV_3D1.F"];
*)

(*
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D2.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D2.F"];
CopyFile[ "dog.f", "DV_3D2.F" ];
SetFileDate["DV_3D2.F";
*)

(*
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D3.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D3.F"];
CopyFile[ "dog.f", "DV_3D3.F" ];
SetFileDate["DV_3D3.F"];
*)

(*
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D4.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D4.F"];
CopyFile[ "dog.f", "DV_3D4.F" ];
SetFileDate["DV_3D4.F"];
*)
