(*^

::[	Information =

	"This is a Mathematica Notebook file.  It contains ASCII text, and can be
	transferred by email, ftp, or other text-file transfer utility.  It should
	be read or edited using a copy of Mathematica or MathReader.  If you 
	received this as email, use your mail application or copy/paste to save 
	everything from the line containing (*^ down to the line containing ^*)
	into a plain text file.  On some systems you may have to give the file a 
	name ending with ".ma" to allow Mathematica to recognize it as a Notebook.
	The line below identifies what version of Mathematica created this file,
	but it can be opened using any other version as well.";

	FrontEndVersion = "X Window System Mathematica Notebook Front End Version 2.2";

	X11StandardFontEncoding; 
	
	fontset = title, inactive, noPageBreakBelow, noPageBreakInGroup, nohscroll, preserveAspect, groupLikeTitle, center, M7, bold, e8,  24, fontName, "times";
	fontset = subtitle, inactive, noPageBreakBelow, noPageBreakInGroup, nohscroll, preserveAspect, groupLikeTitle, center, M7, bold, e6,  18, fontName, "times";
	fontset = subsubtitle, inactive, noPageBreakBelow, noPageBreakInGroup, nohscroll, preserveAspect, groupLikeTitle, center, M7, italic, e6,  14, fontName, "times";
	fontset = section, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, grayBox, M22, bold, a20,  18, fontName, "times";
	fontset = subsection, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, blackBox, M19, bold, a15,  14, fontName, "times";
	fontset = subsubsection, inactive, noPageBreakBelow, nohscroll, preserveAspect, groupLikeSection, whiteBox, M18, bold, a12,  12, fontName, "times";
	fontset = text, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = smalltext, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  10, fontName, "times";
	fontset = input, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeInput, M42, N23, bold,  12, fontName, "courier";
	fontset = output, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23, L-5,  12, fontName, "courier";
	fontset = message, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23,  12, fontName, "courier";
	fontset = print, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23,  12, fontName, "courier";
	fontset = info, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeOutput, M42, N23,  12, fontName, "courier";
	fontset = postscript, PostScript, formatAsPostScript, output, inactive, noPageBreakInGroup, nowordwrap, preserveAspect, groupLikeGraphics, M7, l34, w282, h287,  12, fontName, "courier";
	fontset = name, inactive, noPageBreakInGroup, nohscroll, preserveAspect, M7, italic, B65535,  10, fontName, "times";
	fontset = header, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7, italic,  12, fontName, "times";
	fontset = leftheader,  12, fontName, "times";
	fontset = footer, inactive, nohscroll, noKeepOnOnePage, preserveAspect, center, M7, italic,  12, fontName, "times";
	fontset = leftfooter,  12, fontName, "times";
	fontset = help, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = clipboard, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = completions, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "courier";
	fontset = special1, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = special2, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = special3, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = special4, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";
	fontset = special5, inactive, nohscroll, noKeepOnOnePage, preserveAspect, M7,  12, fontName, "times";paletteColors = 128; automaticGrouping; currentKernel; 
]
:[font = section; inactive; preserveAspect; startGroup]
definition of tau
:[font = input; preserveAspect; startGroup]
tauxxdef = tauxx[i_ + half, j_ ] -> 
		2*mu[i+half,j]*dudx[i+half,j]
:[font = output; output; inactive; preserveAspect; endGroup]
tauxx[half + (i_), j_] -> 2*dudx[half + i, j]*mu[half + i, j]
;[o]
tauxx[half + (i_), j_] -> 2 dudx[half + i, j] mu[half + i, j]
:[font = input; preserveAspect; startGroup]
tauyydef = tauyy[i_, j_ + half] ->
		2*mu[i,j+half] * dvdy[i,j+half]
:[font = output; output; inactive; preserveAspect; endGroup]
tauyy[i_, half + (j_)] -> 2*dvdy[i, half + j]*mu[i, half + j]
;[o]
tauyy[i_, half + (j_)] -> 2 dvdy[i, half + j] mu[i, half + j]
:[font = input; preserveAspect; startGroup]
tauxydef = tauxy[i_,j_] -> mu[i,j]*(dudy[i,j]+dvdx[i,j])
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "tauxydef"
     is similar to existing symbols {tauxxdef, tauyydef}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "tauxy"
     is similar to existing symbols {tauxx, tauyy}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dudy"
     is similar to existing symbols {dudx, dvdy}.
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell
     will be suppressed during this calculation.
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
tauxy[i_, j_] -> (dudy[i, j] + dvdx[i, j])*mu[i, j]
;[o]
tauxy[i_, j_] -> (dudy[i, j] + dvdx[i, j]) mu[i, j]
:[font = section; inactive; preserveAspect; startGroup]
definitions of derivatives
:[font = input; preserveAspect; startGroup]
dudydef1 = dudy[i_,j_+half] -> (u[i,j+1]-u[i,j])/hy
:[font = output; output; inactive; preserveAspect; endGroup]
dudy[i_, half + (j_)] -> (-u[i, j] + u[i, 1 + j])/hy
;[o]
                         -u[i, j] + u[i, 1 + j]
dudy[i_, half + (j_)] -> ----------------------
                                   hy
:[font = input; preserveAspect; startGroup]
dvdxdef1 = dvdx[i_+half,j_] -> (v[i+1,j]-v[i,j])/hx
:[font = output; output; inactive; preserveAspect; endGroup]
dvdx[half + (i_), j_] -> (-v[i, j] + v[1 + i, j])/hx
;[o]
                         -v[i, j] + v[1 + i, j]
dvdx[half + (i_), j_] -> ----------------------
                                   hx
:[font = input; preserveAspect; startGroup]
dudydef2 = dudy[i_+half,j_] -> 
		(u[i,j+1]-u[i,j-1]+u[i+1,j+1]-u[i+1,j-1])/(4*hy)
:[font = output; output; inactive; preserveAspect; endGroup]
dudy[half + (i_), j_] -> 
 
  (-u[i, -1 + j] + u[i, 1 + j] - u[1 + i, -1 + j] + 
 
     u[1 + i, 1 + j])/(4*hy)
;[o]
dudy[half + (i_), j_] -> 
 
  (-u[i, -1 + j] + u[i, 1 + j] - u[1 + i, -1 + j] + 
 
     u[1 + i, 1 + j]) / (4 hy)
:[font = input; preserveAspect; startGroup]
dvdxdef2 = dvdx[i_,j_+half] ->
		(v[i+1,j+1]-v[i-1,j+1]+v[i+1,j]-v[i-1,j])/(4*hx)
:[font = output; output; inactive; preserveAspect; endGroup]
dvdx[i_, half + (j_)] -> 
 
  (-v[-1 + i, j] - v[-1 + i, 1 + j] + v[1 + i, j] + 
 
     v[1 + i, 1 + j])/(4*hx)
;[o]
dvdx[i_, half + (j_)] -> 
 
  (-v[-1 + i, j] - v[-1 + i, 1 + j] + v[1 + i, j] + 
 
     v[1 + i, 1 + j]) / (4 hx)
:[font = input; preserveAspect; startGroup]
dudxdef = dudx[i_+half,j_] -> 
		(u[i+1,j]-u[i,j])/hx
:[font = output; output; inactive; preserveAspect; endGroup]
dudx[half + (i_), j_] -> (-u[i, j] + u[1 + i, j])/hx
;[o]
                         -u[i, j] + u[1 + i, j]
dudx[half + (i_), j_] -> ----------------------
                                   hx
:[font = input; preserveAspect; startGroup]
dvdydef = dvdy[i_,j_+half] ->
		(v[i,j+1]-v[i,j])/hy
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dvdy[i_, half + (j_)] -> (-v[i, j] + v[i, 1 + j])/hy
;[o]
                         -v[i, j] + v[i, 1 + j]
dvdy[i_, half + (j_)] -> ----------------------
                                   hy
:[font = section; inactive; preserveAspect; startGroup]
definitions used to test taylor expansions
:[font = input; preserveAspect; startGroup]
taylorudef = u[i_,j_]->
			U[x0,y0]+
			DuDx[x0,y0]*((i+1/2)*hx-x0)+
			DuDy[x0,y0]*((j+1/2)*hy-y0)+
			D2uDy2[x0,y0]/2*((i+1/2)*hx-x0)^2+
			D2uDx2[x0,y0]/2*((j+1/2)*hy-y0)^2+
			D2uDxDy[x0,y0]*((i+1/2)*hx-x0)*((j+1/2)*hy-y0)
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "DuDy"
     is similar to existing symbol "DuDx".
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "D2uDx2"
     is similar to existing symbol "D2uDy2".
:[font = output; output; inactive; preserveAspect; endGroup]
u[i_, j_] -> 
 
  (hx*(1/2 + i) - x0)*DuDx[x0, y0] + 
 
   (hy*(1/2 + j) - y0)*DuDy[x0, y0] + 
 
   (hx*(1/2 + i) - x0)*(hy*(1/2 + j) - y0)*D2uDxDy[x0, y0] + 
 
   ((hy*(1/2 + j) - y0)^2*D2uDx2[x0, y0])/2 + 
 
   ((hx*(1/2 + i) - x0)^2*D2uDy2[x0, y0])/2 + U[x0, y0]
;[o]
u[i_, j_] -> 
 
       1
  (hx (- + i) - x0) DuDx[x0, y0] + 
       2
 
        1
   (hy (- + j) - y0) DuDy[x0, y0] + 
        2
 
        1                 1
   (hx (- + i) - x0) (hy (- + j) - y0) D2uDxDy[x0, y0] + 
        2                 2
 
        1           2
   (hy (- + j) - y0)  D2uDx2[x0, y0]
        2
   --------------------------------- + 
                   2
 
        1           2
   (hx (- + i) - x0)  D2uDy2[x0, y0]
        2
   --------------------------------- + U[x0, y0]
                   2
:[font = input; preserveAspect]
dog = dudy[i+half,j] //. {dudydef1,dudydef2,taylorudef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,x0->hx,y0->hy/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx, hy/2]
;[o]
         hy
DuDy[hx, --]
         2
:[font = input; preserveAspect]
dog = dudy[i,j+half] //. {dudydef1,dudydef2,taylorudef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,x0->hx/2,y0->hy}
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx/2, hy]
;[o]
     hx
DuDy[--, hy]
     2
:[font = input; preserveAspect]
dog = dudx[i+half,j] //. {dudxdef,taylorudef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,x0->hx,y0->hy/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DuDx[hx, hy/2]
;[o]
         hy
DuDx[hx, --]
         2
:[font = input; preserveAspect; startGroup]
taylorvdef = v[i_,j_]->
			V[x0,y0]+
			DvDx[x0,y0]*((i+1/2)*hx-x0)+
			DvDy[x0,y0]*((j+1/2)*hy-y0)+
			D2vDy2[x0,y0]/2*((i+1/2)*hx-x0)^2+
			D2vDx2[x0,y0]/2*((j+1/2)*hy-y0)^2+
			D2vDxDy[x0,y0]*((i+1/2)*hx-x0)*((j+1/2)*hy-y0)
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "taylorvdef"
     is similar to existing symbol "taylorudef".
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "DvDx"
     is similar to existing symbol "DuDx".
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "DvDy"
     is similar to existing symbols {DuDy, DvDx}.
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "D2vDy2"
     is similar to existing symbol "D2uDy2".
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell1
     will be suppressed during this calculation.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "D2vDx2"
     is similar to existing symbols {D2uDx2, D2vDy2}.
:[font = output; output; inactive; preserveAspect; endGroup]
v[i_, j_] -> 
 
  (hx*(1/2 + i) - x0)*DvDx[x0, y0] + 
 
   (hy*(1/2 + j) - y0)*DvDy[x0, y0] + 
 
   (hx*(1/2 + i) - x0)*(hy*(1/2 + j) - y0)*D2vDxDy[x0, y0] + 
 
   ((hy*(1/2 + j) - y0)^2*D2vDx2[x0, y0])/2 + 
 
   ((hx*(1/2 + i) - x0)^2*D2vDy2[x0, y0])/2 + V[x0, y0]
;[o]
v[i_, j_] -> 
 
       1
  (hx (- + i) - x0) DvDx[x0, y0] + 
       2
 
        1
   (hy (- + j) - y0) DvDy[x0, y0] + 
        2
 
        1                 1
   (hx (- + i) - x0) (hy (- + j) - y0) D2vDxDy[x0, y0] + 
        2                 2
 
        1           2
   (hy (- + j) - y0)  D2vDx2[x0, y0]
        2
   --------------------------------- + 
                   2
 
        1           2
   (hx (- + i) - x0)  D2vDy2[x0, y0]
        2
   --------------------------------- + V[x0, y0]
                   2
:[font = input; preserveAspect]
dog = dvdy[i,j+half] //. {dvdydef, taylorvdef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,x0->hx/2,y0->hy}
:[font = output; output; inactive; preserveAspect; endGroup]
DvDy[hx/2, hy]
;[o]
     hx
DvDy[--, hy]
     2
:[font = input; preserveAspect]
dog = dvdx[i,j+half] //. {dvdxdef1,dvdxdef2,taylorvdef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,x0->hx/2,y0->hy}
:[font = output; output; inactive; preserveAspect; endGroup]
DvDx[hx/2, hy]
;[o]
     hx
DvDx[--, hy]
     2
:[font = input; preserveAspect]
dog = dvdx[i+half,j] //. {dvdxdef1,dvdxdef2,taylorvdef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,x0->hx,y0->hy/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DvDx[hx, hy/2]
;[o]
         hy
DvDx[hx, --]
         2
:[font = input; preserveAspect; startGroup]
taylormudef = mu[i_,j_]->
			MU[x0,y0]+
			DmuDx[x0,y0]*((i+1/2)*hx-x0)+
			DmuDy[x0,y0]*((j+1/2)*hy-y0)

:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "taylormudef"
     is similar to existing symbol "taylorudef".
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "DmuDx"
     is similar to existing symbol "DuDx".
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "DmuDy"
     is similar to existing symbols {DmuDx, DuDy}.
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
mu[i_, j_] -> 
 
  (hx*(1/2 + i) - x0)*DmuDx[x0, y0] + 
 
   (hy*(1/2 + j) - y0)*DmuDy[x0, y0] + MU[x0, y0]
;[o]
mu[i_, j_] -> 
 
       1
  (hx (- + i) - x0) DmuDx[x0, y0] + 
       2
 
        1
   (hy (- + j) - y0) DmuDy[x0, y0] + MU[x0, y0]
        2
:[font = section; inactive; preserveAspect; startGroup]
tests
:[font = input; preserveAspect]
dog = tauxx[i+half,j] //.
		{tauxxdef,dudxdef,taylormudef,taylorudef};
:[font = input; preserveAspect; startGroup]
 dog //. {half->1/2,i->0,j->0,x0->hx,y0 -> hy/2 } 
:[font = output; output; inactive; preserveAspect; endGroup]
2*DuDx[hx, hy/2]*MU[hx, hy/2]
;[o]
           hy         hy
2 DuDx[hx, --] MU[hx, --]
           2          2
:[font = input; preserveAspect]
dog = tauyy[i,j+half] //. 
		{tauyydef,dvdydef,taylormudef,taylorvdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,x0->hx/2,y0->hy}
:[font = output; output; inactive; preserveAspect; endGroup]
2*DvDy[hx/2, hy]*MU[hx/2, hy]
;[o]
       hx         hx
2 DvDy[--, hy] MU[--, hy]
       2          2
:[font = input; preserveAspect]
dog = tauxy[i,j+half] //.
		{tauxydef,dudydef1,dudydef2,dvdxdef1,dvdxdef2,
		 taylormudef,taylorudef,taylorvdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,x0->hx/2,y0->hy}
:[font = output; output; inactive; preserveAspect; endGroup]
(DuDy[hx/2, hy] + DvDx[hx/2, hy])*MU[hx/2, hy]
;[o]
      hx             hx          hx
(DuDy[--, hy] + DvDx[--, hy]) MU[--, hy]
      2              2           2
:[font = input; preserveAspect]
dog = tauxy[i+half,j] //.
		{tauxydef,dudydef1,dudydef2,dvdxdef1,dvdxdef2,
		 taylormudef,taylorudef,taylorvdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,x0->hx,y0->hy/2}
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
(DuDy[hx, hy/2] + DvDx[hx, hy/2])*MU[hx, hy/2]
;[o]
          hy             hy          hy
(DuDy[hx, --] + DvDx[hx, --]) MU[hx, --]
          2              2           2
:[font = section; inactive; preserveAspect; startGroup]
definitions used for fortran output
:[font = input; preserveAspect; startGroup]
murepl1 = mu[i_,j_+half] -> muY[i,j+1]
:[font = output; output; inactive; preserveAspect; endGroup]
mu[i_, half + (j_)] -> muY[i, 1 + j]
;[o]
mu[i_, half + (j_)] -> muY[i, 1 + j]
:[font = input; preserveAspect; startGroup]
murepl2 = mu[i_+half,j_] -> muX[i+1,j]
:[font = output; output; inactive; preserveAspect; endGroup]
mu[half + (i_), j_] -> muX[1 + i, j]
;[o]
mu[half + (i_), j_] -> muX[1 + i, j]
:[font = input; preserveAspect; startGroup]
urepl = u[i_,j_] -> U[i,j,1]
:[font = output; output; inactive; preserveAspect; endGroup]
u[i_, j_] -> U[i, j, 1]
;[o]
u[i_, j_] -> U[i, j, 1]
:[font = input; preserveAspect; startGroup]
vrepl = v[i_,j_] -> U[i,j,2]
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "vrepl"
     is similar to existing symbol "urepl".
:[font = output; output; inactive; preserveAspect; endGroup]
v[i_, j_] -> U[i, j, 2]
;[o]
v[i_, j_] -> U[i, j, 2]
:[font = input; preserveAspect; startGroup]
tauxy[i,j+half] //. 
	{tauxydef,dudydef1,dudydef2,dvdxdef1,dvdxdef2,
	 murepl1,murepl2,urepl,vrepl}
:[font = output; output; inactive; preserveAspect; endGroup]
muY[i, 1 + j]*((-U[i, j, 1] + U[i, 1 + j, 1])/hy + 
 
    (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
       U[1 + i, j, 2] + U[1 + i, 1 + j, 2])/(4*hx))
;[o]
               -U[i, j, 1] + U[i, 1 + j, 1]
muY[i, 1 + j] (---------------------------- + 
                            hy
 
    (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
       U[1 + i, j, 2] + U[1 + i, 1 + j, 2]) / (4 hx))
:[font = input; preserveAspect; startGroup]
(tauxy[i,j+half]-tauxy[i,j-1+half]+
 tauxx[i+half,j]-tauxx[i-1+half,j]) //. 
	{tauxxdef,tauxydef,
	 dudydef1,dudydef2,dvdxdef1,dvdxdef2,dudxdef,dvdydef,
	 murepl1,murepl2,urepl,vrepl}
:[font = output; output; inactive; preserveAspect; endGroup]
(-2*muX[i, j]*(-U[-1 + i, j, 1] + U[i, j, 1]))/hx + 
 
  (2*muX[1 + i, j]*(-U[i, j, 1] + U[1 + i, j, 1]))/hx - 
 
  muY[i, j]*((-U[i, -1 + j, 1] + U[i, j, 1])/hy + 
 
     (-U[-1 + i, -1 + j, 2] - U[-1 + i, j, 2] + 
 
        U[1 + i, -1 + j, 2] + U[1 + i, j, 2])/(4*hx)) + 
 
  muY[i, 1 + j]*((-U[i, j, 1] + U[i, 1 + j, 1])/hy + 
 
     (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
        U[1 + i, j, 2] + U[1 + i, 1 + j, 2])/(4*hx))
;[o]
-2 muX[i, j] (-U[-1 + i, j, 1] + U[i, j, 1])
-------------------------------------------- + 
                     hx
 
  2 muX[1 + i, j] (-U[i, j, 1] + U[1 + i, j, 1])
  ---------------------------------------------- - 
                        hx
 
             -U[i, -1 + j, 1] + U[i, j, 1]
  muY[i, j] (----------------------------- + 
                          hy
 
     (-U[-1 + i, -1 + j, 2] - U[-1 + i, j, 2] + 
 
        U[1 + i, -1 + j, 2] + U[1 + i, j, 2]) / (4 hx)) + 
 
                 -U[i, j, 1] + U[i, 1 + j, 1]
  muY[i, 1 + j] (---------------------------- + 
                              hy
 
     (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
        U[1 + i, j, 2] + U[1 + i, 1 + j, 2]) / (4 hx))
:[font = subsection; inactive; preserveAspect; startGroup]
Extract stencil information from operator
:[font = input; preserveAspect; startGroup]
dog = alpha*a[i,j]*U[i,j,1]-
	beta*(hy*(tauxx[i+half,j]-tauxx[i-1+half,j]) +
	      hx*(tauxy[i,j+half]-tauxy[i,j-1+half]))/vol //.
	      {tauxxdef,tauxydef,
	 		dudydef1,dudydef2,dvdxdef1,dvdxdef2,
	 		dudxdef,dvdydef,
	 		murepl1,murepl2,urepl,vrepl,vol->hx*hy}
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "beta"
     is similar to existing symbol "Beta".
:[font = output; output; inactive; preserveAspect; endGroup]
alpha*a[i, j]*U[i, j, 1] - 
 
  (beta*(hy*((-2*muX[i, j]*(-U[-1 + i, j, 1] + U[i, j, 1]))/
 
           hx + (2*muX[1 + i, j]*
 
             (-U[i, j, 1] + U[1 + i, j, 1]))/hx) + 
 
       hx*(-(muY[i, j]*((-U[i, -1 + j, 1] + U[i, j, 1])/
 
                hy + (-U[-1 + i, -1 + j, 2] - 
 
                  U[-1 + i, j, 2] + U[1 + i, -1 + j, 2] + 
 
                  U[1 + i, j, 2])/(4*hx))) + 
 
          muY[i, 1 + j]*((-U[i, j, 1] + U[i, 1 + j, 1])/
 
              hy + (-U[-1 + i, j, 2] - 
 
                U[-1 + i, 1 + j, 2] + U[1 + i, j, 2] + 
 
                U[1 + i, 1 + j, 2])/(4*hx)))))/(hx*hy)
;[o]
alpha a[i, j] U[i, j, 1] - 
 
             -2 muX[i, j] (-U[-1 + i, j, 1] + U[i, j, 1])
  (beta (hy (-------------------------------------------- + 
                                  hx
 
          2 muX[1 + i, j] (-U[i, j, 1] + U[1 + i, j, 1])
          ----------------------------------------------) + 
                                hx
 
                        -U[i, -1 + j, 1] + U[i, j, 1]
       hx (-(muY[i, j] (----------------------------- + 
                                     hy
 
               (-U[-1 + i, -1 + j, 2] - U[-1 + i, j, 2] + 
 
                  U[1 + i, -1 + j, 2] + U[1 + i, j, 2]) / 
 
                (4 hx))) + 
 
                         -U[i, j, 1] + U[i, 1 + j, 1]
          muY[i, 1 + j] (---------------------------- + 
                                      hy
 
             (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
                U[1 + i, j, 2] + U[1 + i, 1 + j, 2]) / 
 
              (4 hx))))) / (hx hy)
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i,j,1]]
:[font = output; output; inactive; preserveAspect; endGroup]
alpha*a[i, j] + (2*beta*muX[i, j])/hx^2 + 
 
  (2*beta*muX[1 + i, j])/hx^2 + (beta*muY[i, j])/hy^2 + 
 
  (beta*muY[i, 1 + j])/hy^2
;[o]
                2 beta muX[i, j]   2 beta muX[1 + i, j]
alpha a[i, j] + ---------------- + -------------------- + 
                        2                    2
                      hx                   hx
 
  beta muY[i, j]   beta muY[i, 1 + j]
  -------------- + ------------------
         2                  2
       hy                 hy
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i+1,j,1]]
:[font = output; output; inactive; preserveAspect; endGroup]
(-2*beta*muX[1 + i, j])/hx^2
;[o]
-2 beta muX[1 + i, j]
---------------------
           2
         hx
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i-1,j,1]]
:[font = output; output; inactive; preserveAspect; endGroup]
(-2*beta*muX[i, j])/hx^2
;[o]
-2 beta muX[i, j]
-----------------
         2
       hx
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i,j+1,1]]
:[font = output; output; inactive; preserveAspect; endGroup]
-((beta*muY[i, 1 + j])/hy^2)
;[o]
  beta muY[i, 1 + j]
-(------------------)
           2
         hy
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i,j-1,1]]
:[font = output; output; inactive; preserveAspect; endGroup]
-((beta*muY[i, j])/hy^2)
;[o]
  beta muY[i, j]
-(--------------)
         2
       hy
:[font = input; preserveAspect; startGroup]
dog = alpha*a[i,j]*U[i,j,2] -
	beta*(hy*(tauxy[i+half,j]-tauxy[i-1+half,j])+
	      hx*(tauyy[i,j+half]-tauyy[i,j-1+half]))/vol //.
	      {tauyydef,tauxydef,
	 		dudydef1,dudydef2,dvdxdef1,dvdxdef2,
	 		dudxdef,dvdydef,
	 		murepl1,murepl2,urepl,vrepl,vol->hx*hy}
:[font = output; output; inactive; preserveAspect; endGroup]
alpha*a[i, j]*U[i, j, 2] - 
 
  (beta*(hx*((-2*muY[i, j]*(-U[i, -1 + j, 2] + U[i, j, 2]))/
 
           hy + (2*muY[i, 1 + j]*
 
             (-U[i, j, 2] + U[i, 1 + j, 2]))/hy) + 
 
       hy*(-(muX[i, j]*((-U[-1 + i, j, 2] + U[i, j, 2])/
 
                hx + (-U[-1 + i, -1 + j, 1] + 
 
                  U[-1 + i, 1 + j, 1] - U[i, -1 + j, 1] + 
 
                  U[i, 1 + j, 1])/(4*hy))) + 
 
          muX[1 + i, j]*((-U[i, j, 2] + U[1 + i, j, 2])/
 
              hx + (-U[i, -1 + j, 1] + U[i, 1 + j, 1] - 
 
                U[1 + i, -1 + j, 1] + U[1 + i, 1 + j, 1])/
 
              (4*hy)))))/(hx*hy)
;[o]
alpha a[i, j] U[i, j, 2] - 
 
             -2 muY[i, j] (-U[i, -1 + j, 2] + U[i, j, 2])
  (beta (hx (-------------------------------------------- + 
                                  hy
 
          2 muY[i, 1 + j] (-U[i, j, 2] + U[i, 1 + j, 2])
          ----------------------------------------------) + 
                                hy
 
                        -U[-1 + i, j, 2] + U[i, j, 2]
       hy (-(muX[i, j] (----------------------------- + 
                                     hx
 
               (-U[-1 + i, -1 + j, 1] + 
 
                  U[-1 + i, 1 + j, 1] - U[i, -1 + j, 1] + 
 
                  U[i, 1 + j, 1]) / (4 hy))) + 
 
                         -U[i, j, 2] + U[1 + i, j, 2]
          muX[1 + i, j] (---------------------------- + 
                                      hx
 
             (-U[i, -1 + j, 1] + U[i, 1 + j, 1] - 
 
                U[1 + i, -1 + j, 1] + U[1 + i, 1 + j, 1]) / 
 
              (4 hy))))) / (hx hy)
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i,j,2]]
:[font = output; output; inactive; preserveAspect; endGroup]
alpha*a[i, j] + (beta*muX[i, j])/hx^2 + 
 
  (beta*muX[1 + i, j])/hx^2 + (2*beta*muY[i, j])/hy^2 + 
 
  (2*beta*muY[i, 1 + j])/hy^2
;[o]
                beta muX[i, j]   beta muX[1 + i, j]
alpha a[i, j] + -------------- + ------------------ + 
                       2                  2
                     hx                 hx
 
  2 beta muY[i, j]   2 beta muY[i, 1 + j]
  ---------------- + --------------------
          2                    2
        hy                   hy
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i+1,j,2]]
:[font = output; output; inactive; preserveAspect; endGroup]
-((beta*muX[1 + i, j])/hx^2)
;[o]
  beta muX[1 + i, j]
-(------------------)
           2
         hx
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i-1,j,2]]
:[font = output; output; inactive; preserveAspect; endGroup]
-((beta*muX[i, j])/hx^2)
;[o]
  beta muX[i, j]
-(--------------)
         2
       hx
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i,j+1,2]]
:[font = output; output; inactive; preserveAspect; endGroup]
(-2*beta*muY[i, 1 + j])/hy^2
;[o]
-2 beta muY[i, 1 + j]
---------------------
           2
         hy
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[dog],U[i,j-1,2]]
:[font = output; output; inactive; preserveAspect; endGroup]
(-2*beta*muY[i, j])/hy^2
;[o]
-2 beta muY[i, j]
-----------------
         2
       hy
:[font = subsubsection; inactive; preserveAspect; startGroup]
dependentCellsNotCovered is a function which returns a logical expression suitable for
inclusion in fortran.  Give an expression, exp, we wish to determine which mesh locations are accessed by the expression.  However, we do not wish to examine
all possible locations, only those outside the grid patch region.  So we provide
a second argument, which is a boolean function taking two arguments.  The combination will give logical expressions testing the mask for cells utilized by the
expression and for which the boolean function, logfunc[il,jl], evaluates as true. The
third argument is the name of the mask array
:[font = input; preserveAspect]
Clear[ dependentCellsNotCovered ]
:[font = input; preserveAspect; endGroup]
dependentCellsNotCovered[exp_ , logfunc_ ,maskfun_] :=
  Module[{cond,lexp},
	cond = False;
	lexp = exp;
	For[il=-2,il<=+2,il++,
      For[jl=-2,jl<=+2,jl++,
        For[kl=1,kl<=2,kl++,
          If[ (logfunc[il,jl]) &&
	        (Coefficient[
		      Expand[ 
		        exp
		      ] ,
		      U[i+il,j+jl,kl]
	        ] =!= 0), cond = cond || (maskfun[i+il,j+jl]>0)
	      ]
        ]
      ]
    ];
    cond
  ]
:[font = subsubsection; inactive; preserveAspect; startGroup]
dependentCellsCovered is more or less the converse of dependentCellsNotCovered.
dependentCellsCovered will return true if all the cells in the expression (properly
restricted by the mask function) ARE covered.  dependentCellsNotCovered returned
true if one or more of the cells were not covered.
:[font = input; preserveAspect]
Clear[ dependentCellsCovered ]
:[font = input; preserveAspect; endGroup]
dependentCellsCovered[exp_ , logfunc_ ,maskfun_] :=
  Module[{cond,lexp},
	cond = True;
	lexp = exp;
	For[il=-2,il<=+2,il++,
      For[jl=-2,jl<=+2,jl++,
        For[kl=1,kl<=2,kl++,
          If[ (logfunc[il,jl]) &&
	        (Coefficient[
		      Expand[ 
		        exp
		      ] ,
		      U[i+il,j+jl,kl]
	        ] =!= 0), cond = cond && (maskfun[i+il,j+jl]==0)
	      ]
        ]
      ]
    ];
    cond
  ]
:[font = subsubsection; inactive; preserveAspect; startGroup]
this is an alternate definition which treats cells with mask .eq. 1 as untrusted,
rather than mask .gt. 0
:[font = input; preserveAspect]
dependentCellsNotCoveredalt[exp_ , logfunc_ ,maskfun_] :=
  Module[{cond,lexp},
	cond = False;
	lexp = exp;
	For[il=-2,il<=+2,il++,
      For[jl=-2,jl<=+2,jl++,
        For[kl=1,kl<=2,kl++,
          If[ (logfunc[il,jl]) &&
	        (Coefficient[
		      Expand[ 
		        exp
		      ] ,
		      U[i+il,j+jl,kl]
	        ] =!= 0), cond = cond || (maskfun[i+il,j+jl]==1)
	      ]
        ]
      ]
    ];
    cond
  ]
:[font = input; preserveAspect; startGroup]
dependentCellsNotCovered[abba*U[i+1,j-1,1]+U[i-1,j-1,2]
							, (#2<0)& , masks  ]
:[font = output; output; inactive; preserveAspect; endGroup]
masks[-1 + i, -1 + j] > 0 || masks[1 + i, -1 + j] > 0
;[o]
masks[-1 + i, -1 + j] > 0 || masks[1 + i, -1 + j] > 0
:[font = input; preserveAspect]
f[x_] := f1*(x-x2)(x-x3)/(x1-x2)/(x1-x3)+
		 f2*(x-x1)(x-x3)/(x2-x1)/(x2-x3)+
		 f3*(x-x1)(x-x2)/(x3-x1)/(x3-x2)
:[font = input; preserveAspect; startGroup]
dog = Simplify[ D[f[x],x] /. {x1->0,x2->-h,x3->-2*h,x->0} ]
:[font = output; output; inactive; preserveAspect; endGroup]
(3*f1 - 4*f2 + f3)/(2*h)
;[o]
3 f1 - 4 f2 + f3
----------------
      2 h
:[font = input; preserveAspect; startGroup]
dog /. {f1->0, f2->h^2, f3->4*h^2 }
:[font = output; output; inactive; preserveAspect; endGroup]
0
;[o]
0
:[font = input; preserveAspect; startGroup]
dog  /. {f1->0,f2->h^3,f3->8*h^3}
:[font = output; output; inactive; preserveAspect; endGroup]
2*h^2
;[o]
   2
2 h
:[font = input; preserveAspect; startGroup]
dog = Simplify[ D[f[x],x,x] /. {x1->0,x2->-h,x3->-2*h,x->0} ]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup; endGroup]
(f1 - 2*f2 + f3)/h^2
;[o]
f1 - 2 f2 + f3
--------------
       2
      h
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for Do One-sided Derivative in X direction.  if sign is positive, 
it means extend the stencil in the positivie x direction.  if negative, extend
in other direction
:[font = input; preserveAspect]
DOneX[u_,i_,j_,sign_] := (-u[i+2,j]+4*u[i+1,j]-3*u[i,j])/
									(2*hx)  /; sign==1
:[font = input; preserveAspect]
DOneX[u_,i_,j_,sign_] := (u[i-2,j]-4*u[i-1,j]+3*u[i,j])/
									(2*hx)  /; sign==-1
:[font = input; preserveAspect]
DOneX[u_,i_,j_,k_,sign_] := (-u[i+2,j,k]+4*u[i+1,j,k]-
							3*u[i,j,k])/(2*hx) /; sign==1
:[font = input; preserveAspect]
DOneX[u_,i_,j_,k_,sign_] := (u[i-2,j,k]-4*u[i-1,j,k]+
							3*u[i,j,k])/(2*hx) /; sign==-1
:[font = input; preserveAspect; startGroup]
Simplify[ DOneX[u,0,0,-1] //. {taylorudef,x0->hx/2,y0->hy/2}]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDx[hx/2, hy/2]
;[o]
     hx  hy
DuDx[--, --]
     2   2
:[font = input; preserveAspect; startGroup]
Simplify[ DOneX[u,0,0,+1] //. {taylorudef,x0->hx/2,y0->hy/2}]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
DuDx[hx/2, hy/2]
;[o]
     hx  hy
DuDx[--, --]
     2   2
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for Do One-sided Derivative in Y direction.  if sign is positive, 
it means extend the stencil in the positivie y direction.  if negative, extend
in other direction
:[font = input; preserveAspect; startGroup]
DOneY[u_,i_,j_,sign_] := (-u[i,j+2]+4*u[i,j+1]-3*u[i,j])/
									(2*hy)  /; sign==1
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "DOneY"
     is similar to existing symbol "DOneX".
:[font = input; preserveAspect]
DOneY[u_,i_,j_,sign_] := (u[i,j-2]-4*u[i,j-1]+3*u[i,j])/
									(2*hy)  /; sign==-1
:[font = input; preserveAspect]
DOneY[u_,i_,j_,k_,sign_] := (-u[i,j+2,k]+4*u[i,j+1,k]-
							3*u[i,j,k])/(2*hy)  /; sign==1
:[font = input; preserveAspect]
DOneY[u_,i_,j_,k_,sign_] := (u[i,j-2,k]-4*u[i,j-1,k]+
							3*u[i,j,k])/(2*hy)  /; sign==-1
:[font = input; preserveAspect; startGroup]
Simplify[ DOneY[u,0,0,-1] //. {taylorudef,x0->hx/2,y0->hy/2}]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx/2, hy/2]
;[o]
     hx  hy
DuDy[--, --]
     2   2
:[font = input; preserveAspect; startGroup]
Simplify[ DOneY[u,0,0,+1] //. {taylorudef,x0->hx/2,y0->hy/2}]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
DuDy[hx/2, hy/2]
;[o]
     hx  hy
DuDy[--, --]
     2   2
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for two sided deriv in x direction.
:[font = input; preserveAspect]
DTwoX[u_,i_,j_] := (u[i+1,j]-u[i-1,j])/(2*hx)
:[font = input; preserveAspect]
DTwoX[u_,i_,j_,k_] := (u[i+1,j,k]-u[i-1,j,k])/(2*hx)
:[font = input; preserveAspect; startGroup]
Simplify[ DTwoX[u,0,0] //. {taylorudef,x0->hx/2,y0->hy/2} ]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDx[hx/2, hy/2]
;[o]
     hx  hy
DuDx[--, --]
     2   2
:[font = input; preserveAspect; startGroup]
DTwoY[u_,i_,j_] := (u[i,j+1]-u[i,j-1])/(2*hy)
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "DTwoY"
     is similar to existing symbol "DTwoX".
:[font = input; preserveAspect]
DTwoY[u_,i_,j_,k_] := (u[i,j+1,k]-u[i,j-1,k])/(2*hy)
:[font = input; preserveAspect; startGroup]
Simplify[ DTwoY[u,0,0] //. {taylorudef,x0->hx/2,y0->hy/2} ]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx/2, hy/2]
;[o]
     hx  hy
DuDy[--, --]
     2   2
:[font = input; preserveAspect; endGroup; endGroup]
DeleteFile[ "dog.mf"];
CopyFile[ "DV_2D.mF" , "dog.mf"];
Splice["dog.mf"];
DeleteFile[ "DV_2D.F"];
CopyFile[ "dog.f", "DV_2D.F" ];
<<"!touch DV_2D.F"

:[font = section; inactive; preserveAspect; startGroup]
Start figuring out how to write stencil
:[font = input; preserveAspect; startGroup]
out1 = alpha*a[i,j]*U[i,j,1]-
	beta*(hy*(tauxx[i+half,j]-tauxx[i-1+half,j]) +
	      hx*(tauxy[i,j+half]-tauxy[i,j-1+half]))/vol //.
	      {tauxxdef,tauxydef,
	 		dudydef1,dudydef2,dvdxdef1,dvdxdef2,
	 		dudxdef,dvdydef,
	 		murepl1,murepl2,urepl,vrepl,vol->hx*hy}
:[font = output; output; inactive; preserveAspect; endGroup]
alpha*a[i, j]*U[i, j, 1] - 
 
  (beta*(hy*((-2*muX[i, j]*(-U[-1 + i, j, 1] + U[i, j, 1]))/
 
           hx + (2*muX[1 + i, j]*
 
             (-U[i, j, 1] + U[1 + i, j, 1]))/hx) + 
 
       hx*(-(muY[i, j]*((-U[i, -1 + j, 1] + U[i, j, 1])/hy + 
 
               (-U[-1 + i, -1 + j, 2] - U[-1 + i, j, 2] + 
 
                  U[1 + i, -1 + j, 2] + U[1 + i, j, 2])/
 
                (4*hx))) + 
 
          muY[i, 1 + j]*((-U[i, j, 1] + U[i, 1 + j, 1])/hy + 
 
             (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
                U[1 + i, j, 2] + U[1 + i, 1 + j, 2])/(4*hx)))
 
       ))/(hx*hy)
;[o]
alpha a[i, j] U[i, j, 1] - 
 
             -2 muX[i, j] (-U[-1 + i, j, 1] + U[i, j, 1])
  (beta (hy (-------------------------------------------- + 
                                  hx
 
          2 muX[1 + i, j] (-U[i, j, 1] + U[1 + i, j, 1])
          ----------------------------------------------) + 
                                hx
 
                        -U[i, -1 + j, 1] + U[i, j, 1]
       hx (-(muY[i, j] (----------------------------- + 
                                     hy
 
               (-U[-1 + i, -1 + j, 2] - U[-1 + i, j, 2] + 
 
                  U[1 + i, -1 + j, 2] + U[1 + i, j, 2]) / 
 
                (4 hx))) + 
 
                         -U[i, j, 1] + U[i, 1 + j, 1]
          muY[i, 1 + j] (---------------------------- + 
                                      hy
 
             (-U[-1 + i, j, 2] - U[-1 + i, 1 + j, 2] + 
 
                U[1 + i, j, 2] + U[1 + i, 1 + j, 2]) / (4 hx)
 
             )))) / (hx hy)
:[font = input; preserveAspect; startGroup]
Coefficient[Expand[out1],U[i,j,1]]
:[font = output; output; inactive; preserveAspect; endGroup]
alpha*a[i, j] + (2*beta*muX[i, j])/hx^2 + 
 
  (2*beta*muX[1 + i, j])/hx^2 + (beta*muY[i, j])/hy^2 + 
 
  (beta*muY[i, 1 + j])/hy^2
;[o]
                2 beta muX[i, j]   2 beta muX[1 + i, j]
alpha a[i, j] + ---------------- + -------------------- + 
                        2                    2
                      hx                   hx
 
  beta muY[i, j]   beta muY[i, 1 + j]
  -------------- + ------------------
         2                  2
       hy                 hy
:[font = subsection; inactive; preserveAspect; startGroup]
setup to use the Format.m package 
:[font = input; preserveAspect]
Off[General::spell,General::spell1];
:[font = input; preserveAspect]
SetOptions[$Output,PageWidth->73];
:[font = input; preserveAspect; startGroup]
<</usr/people/wyc/_math/MathSource/Format/Format.m
:[font = message; inactive; preserveAspect]
exp::shdw: Warning: Symbol exp appears in multiple contexts 
    {Format`, Global`}; definitions in context Format`
     may shadow or be shadowed by other definitions.
:[font = message; inactive; preserveAspect; endGroup]
sign::shdw: Warning: Symbol sign appears in multiple contexts 
    {Format`, Global`}; definitions in context Format`
     may shadow or be shadowed by other definitions.
:[font = input; preserveAspect]
stencilDef[exp_,kin_] := 
	Block[{iind,il,jl,kl,list1,list2,indices,tmp,indtrn},
		indices={
					{ 0, 0,KO},
					{-1, 0,KW},
					{ 0,-1,KS},
					{-1,-1,KSW},
					{-1,+1,KNW}
				};
		list1={};
		list2={};
		For[iind=1,iind<=5,iind++,
			il = indices[[iind,1]];
			jl = indices[[iind,2]];
			indtrn = indices[[iind,3]];
			For[kl=1,kl<=2,kl++,
				tmp=Coefficient[Expand[exp],U[i+il,j+jl,kl]];
				If[tmp=!=0,
					list1=Append[list1,coef[i,j,indtrn,kin,kl]];
					list2=Append[list2,+
					  Coefficient[Expand[exp],U[i+il,j,kl]]] 
				]
			]
		];
		FortranAssign[Evaluate[list1],list2,
				AssignToArray->{muY,muX,a,coef,polycoef}]
	]

							
:[font = input; preserveAspect]
stencilInc[exp_,kin_] := 
	Block[{iind,il,jl,kl,list1,list2,indices,tmp,indtrn},
		indices={
					{ 0, 0,KO},
					{-1, 0,KW},
					{ 0,-1,KS},
					{-1,-1,KSW},
					{-1,+1,KNW}
				};
		list1={};
		list2={};
		For[iind=1,iind<=5,iind++,
			il = indices[[iind,1]];
			jl = indices[[iind,2]];
			indtrn = indices[[iind,3]];
			For[kl=1,kl<=2,kl++,
				tmp=Coefficient[Expand[exp],U[i+il,j+jl,kl]];
				If[tmp=!=0,
					list1=Append[list1,coef[i,j,indtrn,kin,kl]];
					list2=Append[list2,coef[i,j,indtrn,kin,kl]+
					  Coefficient[Expand[exp],U[i+il,j,kl]]] 
				]
			]
		];
		FortranAssign[Evaluate[list1],list2,
				AssignToArray->{muY,muX,a,coef,polycoef}]
	]
	

							
:[font = input; preserveAspect; endGroup]
stencilDec[exp_,kin_] := 
	Block[{iind,il,jl,kl,list1,list2,indices,tmp,indtrn},
		indices={
					{ 0, 0,KO},
					{-1, 0,KW},
					{ 0,-1,KS},
					{-1,-1,KSW},
					{-1,+1,KNW}
				};
		list1={};
		list2={};
		For[iind=1,iind<=5,iind++,
			il = indices[[iind,1]];
			jl = indices[[iind,2]];
			indtrn = indices[[iind,3]];
			For[kl=1,kl<=2,kl++,
				tmp=Coefficient[Expand[exp],U[i+il,j+jl,kl]];
				If[tmp=!=0,
					list1=Append[list1,coef[i,jl,indtrn,kin,kl]];
					list2=Append[list2,coef[i,j,indtrn,kin,kl]-
					  Coefficient[Expand[exp],U[i+il,j,kl]]] 
				]
			]
		];
		FortranAssign[Evaluate[list1],list2,
				AssignToArray->{muY,muX,a,coef,polycoef}]
	]
	

							
:[font = subsection; inactive; preserveAspect; startGroup]
This is slight variation on dependentCellsNotCovered.  Uses Format.m
:[font = input; preserveAspect]
Clear[ depCellNotCovered ]
:[font = input; preserveAspect]
dependentCellIsCovered[exp_ , logfunc_ ,maskfun_] :=
  Module[{cond,lexp},
	cond = False;
	lexp = exp;
	For[il=-2,il<=+2,il++,
      For[jl=-2,jl<=+2,jl++,
        For[kl=1,kl<=2,kl++,
          If[ (logfunc[il,jl]) &&
	        (Coefficient[
		      Expand[ 
		        exp
		      ] ,
		      U[i+il,j+jl,kl]
	        ] =!= 0), cond = cond || (maskfun[i+il,j+jl]>0)
	      ]
        ]
      ]
    ];
    FortranAssign[ cond, AssignToArray->{maskray} ]
  ]
:[font = input; preserveAspect]
Clear[stencoef]
:[font = input; preserveAspect; startGroup]
stencilDef[out1,1]
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
        coef(i,j,KO,1,1)=alpha*a(i,j)+2.d0*beta*muX(i,j)/hx
     &  **2+2.d0*beta*muX(1+i,j)/hx**2+beta*muY(i,j)/hy**2+
     &  beta*muY(i,1+j)/hy**2
        coef(i,j,KW,1,1)=-2.d0*beta*muX(i,j)/hx**2
        coef(i,j,KW,1,2)=-2.5d-1*beta*muY(i,j)/(hx*hy)+2.5d
     &  -1*beta*muY(i,1+j)/(hx*hy)
        coef(i,j,KS,1,1)=alpha*a(i,j)+2.d0*beta*muX(i,j)/hx
     &  **2+2.d0*beta*muX(1+i,j)/hx**2+beta*muY(i,j)/hy**2+
     &  beta*muY(i,1+j)/hy**2
        coef(i,j,KSW,1,2)=-2.5d-1*beta*muY(i,j)/(hx*hy)+2.5
     &  d-1*beta*muY(i,1+j)/(hx*hy)
        coef(i,j,KNW,1,2)=-2.5d-1*beta*muY(i,j)/(hx*hy)+2.5
     &  d-1*beta*muY(i,1+j)/(hx*hy)
;[o]
        coef(i,j,KO,1,1)=alpha*a(i,j)+2.d0*beta*muX(i,j)/hx
     &  **2+2.d0*beta*muX(1+i,j)/hx**2+beta*muY(i,j)/hy**2+
     &  beta*muY(i,1+j)/hy**2
        coef(i,j,KW,1,1)=-2.d0*beta*muX(i,j)/hx**2
        coef(i,j,KW,1,2)=-2.5d-1*beta*muY(i,j)/(hx*hy)+2.5d
     &  -1*beta*muY(i,1+j)/(hx*hy)
        coef(i,j,KS,1,1)=alpha*a(i,j)+2.d0*beta*muX(i,j)/hx
     &  **2+2.d0*beta*muX(1+i,j)/hx**2+beta*muY(i,j)/hy**2+
     &  beta*muY(i,1+j)/hy**2
        coef(i,j,KSW,1,2)=-2.5d-1*beta*muY(i,j)/(hx*hy)+2.5
     &  d-1*beta*muY(i,1+j)/(hx*hy)
        coef(i,j,KNW,1,2)=-2.5d-1*beta*muY(i,j)/(hx*hy)+2.5
     &  d-1*beta*muY(i,1+j)/(hx*hy)
:[font = input; preserveAspect; startGroup]
dog1 = Expand[ vol*(alpha*a[i,j]*U[i,j,1]-
	beta*(hy*(tauxx[i+half,j]-tauxx[i-1+half,j]) +
	      hx*(tauxy[i,j+half]-tauxy[i,j-1+half]))/vol ) //.
	      {tauxxdef,tauxydef,
	 		murepl1,murepl2,urepl,vrepl,vol->hx*hy} ] //InputForm
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
InputForm[2*beta*hy*dudx[-1 + half + i, j]*muX[i, j] - 
 
   2*beta*hy*dudx[half + i, j]*muX[1 + i, j] + 
 
   beta*hx*dudy[i, -1 + half + j]*muY[i, j] + 
 
   beta*hx*dvdx[i, -1 + half + j]*muY[i, j] - 
 
   beta*hx*dudy[i, half + j]*muY[i, 1 + j] - 
 
   beta*hx*dvdx[i, half + j]*muY[i, 1 + j] + 
 
   alpha*hx*hy*a[i, j]*U[i, j, 1]]
;[o]
2*beta*hy*dudx[-1 + half + i, j]*muX[i, j] - 
  2*beta*hy*dudx[half + i, j]*muX[1 + i, j] + 
  beta*hx*dudy[i, -1 + half + j]*muY[i, j] + 
  beta*hx*dvdx[i, -1 + half + j]*muY[i, j] - 
  beta*hx*dudy[i, half + j]*muY[i, 1 + j] - 
  beta*hx*dvdx[i, half + j]*muY[i, 1 + j] + 
  alpha*hx*hy*a[i, j]*U[i, j, 1]
:[font = input; preserveAspect; startGroup]
Splice[ "dog.mf",FormatType->OutputForm]
:[font = message; inactive; preserveAspect]
Splice::splicx: 
   Syntax error in Mathematica input 
    dependentCellsNotCovered[DOneX[u,i,j+1,+1]//..
       allDerivAllUV,
                                                 
       Function[{i,j},(j>0),mask]
:[font = output; output; inactive; preserveAspect; endGroup]
Splice["dog.mf", FormatType -> OutputForm]
;[o]
Splice[dog.mf, FormatType -> OutputForm]
:[font = input; preserveAspect; startGroup]
Splice[ "ratdog.mf",FormatType->OutputForm]
:[font = output; output; inactive; preserveAspect; endGroup]
"ratdog.mf"
;[o]
ratdog.mf
:[font = input; preserveAspect; startGroup]
exntop+exnbot
:[font = output; output; inactive; preserveAspect; endGroup]
(beta*muY[i, j]*(-U[-1 + i, -1 + j, 2] + 
 
       U[1 + i, -1 + j, 2]))/4 + 
 
  (beta*muY[i, j]*(-U[-1 + i, j, 2] + U[1 + i, j, 2]))/4
;[o]
(beta muY[i, j] (-U[-1 + i, -1 + j, 2] + 
 
       U[1 + i, -1 + j, 2])) / 4 + 
 
  beta muY[i, j] (-U[-1 + i, j, 2] + U[1 + i, j, 2])
  --------------------------------------------------
                          4
:[font = input; preserveAspect; startGroup]
dependentCellIsCovered[exntop+exnbot,indexcond,maskray]
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
        (masklox(-1+i,-1+j).gt.0).or.(masklox(-1+i,j).gt.0)
;[o]
        (masklox(-1+i,-1+j).gt.0).or.(masklox(-1+i,j).gt.0)
:[font = input; preserveAspect; startGroup]
Expand[exntop]
:[font = output; output; inactive; preserveAspect; endGroup]
-(beta*muY[i, j]*U[-1 + i, j, 2])/4 + 
 
  (beta*muY[i, j]*U[1 + i, j, 2])/4
;[o]
-(beta muY[i, j] U[-1 + i, j, 2])
--------------------------------- + 
                4
 
  beta muY[i, j] U[1 + i, j, 2]
  -----------------------------
                4
:[font = input; preserveAspect; startGroup]
outdog = 2*beta*hy*dudx[-1+half+i,j] //. 
			{dudxdef,u[i-1,j]->u[i,j]}
:[font = output; output; inactive; preserveAspect; endGroup; endGroup; endGroup]
0
;[o]
0
^*)