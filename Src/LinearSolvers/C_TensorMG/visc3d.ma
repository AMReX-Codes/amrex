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
:[font = subtitle; inactive; preserveAspect; startGroup]
a clone of visc2d.ma  which
should work in 3D
:[font = section; inactive; preserveAspect; startGroup]
definition of tau
:[font = input; preserveAspect; startGroup]
tauxxdef = tauxx[i_ + half, j_ ,k_] -> 
		2*mu[i+half,j,k]*dudx[i+half,j,k]
:[font = output; output; inactive; preserveAspect; endGroup]
tauxx[half + (i_), j_, k_] -> 
 
  2*dudx[half + i, j, k]*mu[half + i, j, k]
;[o]
tauxx[half + (i_), j_, k_] -> 
 
  2 dudx[half + i, j, k] mu[half + i, j, k]
:[font = input; preserveAspect; startGroup]
tauyydef = tauyy[i_, j_ + half,k_] ->
		2*mu[i,j+half,k] * dvdy[i,j+half,k]
:[font = output; output; inactive; preserveAspect; endGroup]
tauyy[i_, half + (j_), k_] -> 
 
  2*dvdy[i, half + j, k]*mu[i, half + j, k]
;[o]
tauyy[i_, half + (j_), k_] -> 
 
  2 dvdy[i, half + j, k] mu[i, half + j, k]
:[font = input; preserveAspect; startGroup]
tauzzdef = tauzz[i_,j_,k_+half] ->
		2*mu[i,j,k+half]*dwdz[i,j,k+half]
:[font = output; output; inactive; preserveAspect; endGroup]
tauzz[i_, j_, half + (k_)] -> 
 
  2*dwdz[i, j, half + k]*mu[i, j, half + k]
;[o]
tauzz[i_, j_, half + (k_)] -> 
 
  2 dwdz[i, j, half + k] mu[i, j, half + k]
:[font = input; preserveAspect; startGroup]
tauxydef = tauxy[i_,j_,k_] -> mu[i,j,k]*
							(dudy[i,j,k]+dvdx[i,j,k])
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
:[font = output; output; inactive; preserveAspect; endGroup]
tauxy[i_, j_, k_] -> 
 
  (dudy[i, j, k] + dvdx[i, j, k])*mu[i, j, k]
;[o]
tauxy[i_, j_, k_] -> 
 
  (dudy[i, j, k] + dvdx[i, j, k]) mu[i, j, k]
:[font = input; preserveAspect; startGroup]
tauxzdef = tauxz[i_,j_,k_] -> mu[i,j,k]*
							(dudz[i,j,k]+dwdx[i,j,k])
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "tauxzdef"
     is similar to existing symbols 
    {tauxxdef, tauxydef, tauzzdef}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "tauxz"
     is similar to existing symbols {tauxx, tauxy, tauzz}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dudz"
     is similar to existing symbols {dudx, dudy, dwdz}.
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell
     will be suppressed during this calculation.
:[font = output; output; inactive; preserveAspect; endGroup]
tauxz[i_, j_, k_] -> 
 
  (dudz[i, j, k] + dwdx[i, j, k])*mu[i, j, k]
;[o]
tauxz[i_, j_, k_] -> 
 
  (dudz[i, j, k] + dwdx[i, j, k]) mu[i, j, k]
:[font = input; preserveAspect; startGroup]
tauyzdef = tauyz[i_,j_,k_] -> mu[i,j,k]*
							(dvdz[i,j,k]+dwdy[i,j,k])
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "tauyzdef"
     is similar to existing symbols 
    {tauxzdef, tauyydef, tauzzdef}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "tauyz"
     is similar to existing symbols {tauxz, tauyy, tauzz}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dvdz"
     is similar to existing symbols {dudz, dvdx, dvdy, 
     dwdz}.
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell
     will be suppressed during this calculation.
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
tauyz[i_, j_, k_] -> 
 
  (dvdz[i, j, k] + dwdy[i, j, k])*mu[i, j, k]
;[o]
tauyz[i_, j_, k_] -> 
 
  (dvdz[i, j, k] + dwdy[i, j, k]) mu[i, j, k]
:[font = section; inactive; preserveAspect; startGroup]
definitions of derivatives
:[font = subsubsection; inactive; preserveAspect; startGroup]
diagonal derivatives
:[font = input; preserveAspect; startGroup]
dudxdef = dudx[i_+half,j_,k_] -> 
		(u[i+1,j,k]-u[i,j,k])/hx
:[font = output; output; inactive; preserveAspect; endGroup]
dudx[half + (i_), j_, k_] -> 
 
  (-u[i, j, k] + u[1 + i, j, k])/hx
;[o]
                             -u[i, j, k] + u[1 + i, j, k]
dudx[half + (i_), j_, k_] -> ----------------------------
                                          hx
:[font = input; preserveAspect; startGroup]
dvdydef = dvdy[i_,j_+half,k_] ->
		(v[i,j+1,k]-v[i,j,k])/hy
:[font = output; output; inactive; preserveAspect; endGroup]
dvdy[i_, half + (j_), k_] -> 
 
  (-v[i, j, k] + v[i, 1 + j, k])/hy
;[o]
                             -v[i, j, k] + v[i, 1 + j, k]
dvdy[i_, half + (j_), k_] -> ----------------------------
                                          hy
:[font = input; preserveAspect; startGroup]
dwdzdef = dwdz[i_,j_,k_+half] ->
		(w[i,j,k+1]-w[i,j,k])/hz
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dwdz[i_, j_, half + (k_)] -> 
 
  (-w[i, j, k] + w[i, j, 1 + k])/hz
;[o]
                             -w[i, j, k] + w[i, j, 1 + k]
dwdz[i_, j_, half + (k_)] -> ----------------------------
                                          hz
:[font = subsubsection; inactive; preserveAspect; startGroup]
dudy
:[font = input; preserveAspect; startGroup]
dudydef1 = dudy[i_,j_+half,k_] -> (u[i,j+1,k]-u[i,j,k])/hy
:[font = output; output; inactive; preserveAspect; endGroup]
dudy[i_, half + (j_), k_] -> 
 
  (-u[i, j, k] + u[i, 1 + j, k])/hy
;[o]
                             -u[i, j, k] + u[i, 1 + j, k]
dudy[i_, half + (j_), k_] -> ----------------------------
                                          hy
:[font = input; preserveAspect; startGroup]
dudydef2 = dudy[i_+half,j_,k_] -> 
		(u[i,j+1,k]-u[i,j-1,k]+u[i+1,j+1,k]-u[i+1,j-1,k])/
													(4*hy)
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dudy[half + (i_), j_, k_] -> 
 
  (-u[i, -1 + j, k] + u[i, 1 + j, k] - 
 
     u[1 + i, -1 + j, k] + u[1 + i, 1 + j, k])/(4*hy)
;[o]
dudy[half + (i_), j_, k_] -> 
 
  (-u[i, -1 + j, k] + u[i, 1 + j, k] - 
 
     u[1 + i, -1 + j, k] + u[1 + i, 1 + j, k]) / (4 hy)
:[font = subsubsection; inactive; preserveAspect; startGroup]
dudz
:[font = input; preserveAspect; startGroup]
dudzdef1 = dudz[i_,j_,k_+half]->(u[i,j,k+1]-u[i,j,k])/hz
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "dudzdef1"
     is similar to existing symbol "dudydef1".
:[font = output; output; inactive; preserveAspect; endGroup]
dudz[i_, j_, half + (k_)] -> 
 
  (-u[i, j, k] + u[i, j, 1 + k])/hz
;[o]
                             -u[i, j, k] + u[i, j, 1 + k]
dudz[i_, j_, half + (k_)] -> ----------------------------
                                          hz
:[font = input; preserveAspect; startGroup]
dudzdef2 = dudz[i_+half,j_,k_] ->
		(u[i,j,k+1]-u[i,j,k-1]+u[i+1,j,k+1]-u[i+1,j,k-1])/
												(4*hz)
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "dudzdef2"
     is similar to existing symbol "dudydef2".
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dudz[half + (i_), j_, k_] -> 
 
  (-u[i, j, -1 + k] + u[i, j, 1 + k] - 
 
     u[1 + i, j, -1 + k] + u[1 + i, j, 1 + k])/(4*hz)
;[o]
dudz[half + (i_), j_, k_] -> 
 
  (-u[i, j, -1 + k] + u[i, j, 1 + k] - 
 
     u[1 + i, j, -1 + k] + u[1 + i, j, 1 + k]) / (4 hz)
:[font = subsubsection; inactive; preserveAspect; startGroup]
dvdx
:[font = input; preserveAspect; startGroup]
dvdxdef1 = dvdx[i_+half,j_,k_] -> (v[i+1,j,k]-v[i,j,k])/hx
:[font = output; output; inactive; preserveAspect; endGroup]
dvdx[half + (i_), j_, k_] -> 
 
  (-v[i, j, k] + v[1 + i, j, k])/hx
;[o]
                             -v[i, j, k] + v[1 + i, j, k]
dvdx[half + (i_), j_, k_] -> ----------------------------
                                          hx
:[font = input; preserveAspect; startGroup]
dvdxdef2 = dvdx[i_,j_+half,k_] ->
		(v[i+1,j+1,k]-v[i-1,j+1,k]+v[i+1,j,k]-v[i-1,j,k])/(4*hx)
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dvdx[i_, half + (j_), k_] -> 
 
  (-v[-1 + i, j, k] - v[-1 + i, 1 + j, k] + 
 
     v[1 + i, j, k] + v[1 + i, 1 + j, k])/(4*hx)
;[o]
dvdx[i_, half + (j_), k_] -> 
 
  (-v[-1 + i, j, k] - v[-1 + i, 1 + j, k] + 
 
     v[1 + i, j, k] + v[1 + i, 1 + j, k]) / (4 hx)
:[font = subsubsection; inactive; preserveAspect; startGroup]
dvdz
:[font = input; preserveAspect; startGroup]
dvdzdef1 = dvdz[i_,j_,k_+half]->(v[i,j,k+1]-v[i,j,k])/hz
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dvdzdef1"
     is similar to existing symbols {dudzdef1, dvdxdef1}.
:[font = output; output; inactive; preserveAspect; endGroup]
dvdz[i_, j_, half + (k_)] -> 
 
  (-v[i, j, k] + v[i, j, 1 + k])/hz
;[o]
                             -v[i, j, k] + v[i, j, 1 + k]
dvdz[i_, j_, half + (k_)] -> ----------------------------
                                          hz
:[font = input; preserveAspect; startGroup]
dvdzdef2 = dvdz[i_,j_+half,k_]->
		(v[i,j,k+1]-v[i,j,k-1]+v[i,j+1,k+1]-v[i,j+1,k-1])/
												(4*hz)
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dvdzdef2"
     is similar to existing symbols {dudzdef2, dvdxdef2}.
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dvdz[i_, half + (j_), k_] -> 
 
  (-v[i, j, -1 + k] + v[i, j, 1 + k] - 
 
     v[i, 1 + j, -1 + k] + v[i, 1 + j, 1 + k])/(4*hz)
;[o]
dvdz[i_, half + (j_), k_] -> 
 
  (-v[i, j, -1 + k] + v[i, j, 1 + k] - 
 
     v[i, 1 + j, -1 + k] + v[i, 1 + j, 1 + k]) / (4 hz)
:[font = subsubsection; inactive; preserveAspect; startGroup]
dwdx
:[font = input; preserveAspect; startGroup]
dwdxdef1 = dwdx[i_+half,j_,k_]->(w[i+1,j,k]-w[i,j,k])/hx
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "dwdxdef1"
     is similar to existing symbol "dvdxdef1".
:[font = output; output; inactive; preserveAspect; endGroup]
dwdx[half + (i_), j_, k_] -> 
 
  (-w[i, j, k] + w[1 + i, j, k])/hx
;[o]
                             -w[i, j, k] + w[1 + i, j, k]
dwdx[half + (i_), j_, k_] -> ----------------------------
                                          hx
:[font = input; preserveAspect; startGroup]
dwdxdef2 = dwdx[i_,j_,k_+half]->
		(w[i+1,j,k]-w[i-1,j,k]+w[i+1,j,k+1]-w[i-1,j,k+1])/
											(4*hx)
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "dwdxdef2"
     is similar to existing symbol "dvdxdef2".
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
dwdx[i_, j_, half + (k_)] -> 
 
  (-w[-1 + i, j, k] - w[-1 + i, j, 1 + k] + 
 
     w[1 + i, j, k] + w[1 + i, j, 1 + k])/(4*hx)
;[o]
dwdx[i_, j_, half + (k_)] -> 
 
  (-w[-1 + i, j, k] - w[-1 + i, j, 1 + k] + 
 
     w[1 + i, j, k] + w[1 + i, j, 1 + k]) / (4 hx)
:[font = subsubsection; inactive; preserveAspect; startGroup]
dwdy
:[font = input; preserveAspect; startGroup]
dwdydef1 = dwdy[i_,j_+half,k_] ->
		(w[i,j+1,k]-w[i,j,k])/hy
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dwdydef1"
     is similar to existing symbols {dudydef1, dwdxdef1}.
:[font = output; output; inactive; preserveAspect; endGroup]
dwdy[i_, half + (j_), k_] -> 
 
  (-w[i, j, k] + w[i, 1 + j, k])/hy
;[o]
                             -w[i, j, k] + w[i, 1 + j, k]
dwdy[i_, half + (j_), k_] -> ----------------------------
                                          hy
:[font = input; preserveAspect; startGroup]
dwdydef2 = dwdy[i_,j_,k_+half] ->
	(w[i,j+1,k]-w[i,j-1,k]+w[i,j+1,k+1]-w[i,j-1,k+1])/(4*hy)
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "dwdydef2"
     is similar to existing symbols {dudydef2, dwdxdef2}.
:[font = output; output; inactive; preserveAspect; endGroup; endGroup; endGroup]
dwdy[i_, j_, half + (k_)] -> 
 
  (-w[i, -1 + j, k] - w[i, -1 + j, 1 + k] + 
 
     w[i, 1 + j, k] + w[i, 1 + j, 1 + k])/(4*hy)
;[o]
dwdy[i_, j_, half + (k_)] -> 
 
  (-w[i, -1 + j, k] - w[i, -1 + j, 1 + k] + 
 
     w[i, 1 + j, k] + w[i, 1 + j, 1 + k]) / (4 hy)
:[font = section; inactive; preserveAspect; startGroup]
definitions used to test taylor expansions
:[font = input; preserveAspect; startGroup]
taylorudef = u[i_,j_,k_]->
			U[x0,y0,z0]+
			DuDx[x0,y0,z0]*((i+1/2)*hx-x0)+
			DuDy[x0,y0,z0]*((j+1/2)*hy-y0)+
			DuDz[x0,y0,z0]*((k+1/2)*hz-z0)+
			D2uDy2[x0,y0,z0]/2*((i+1/2)*hx-x0)^2+
			D2uDx2[x0,y0,z0]/2*((j+1/2)*hy-y0)^2+
			D2uDz2[x0,y0,z0]/2*((k+1/2)*hz-z0)^2+
			D2uDxDy[x0,y0,z0]*((i+1/2)*hx-x0)*((j+1/2)*hy-y0)+
			D2uDxDz[x0,y0,z0]*((i+1/2)*hx-x0)*((k+1/2)*hz-z0)+
			D2uDyDz[x0,y0,z0]*((j+1/2)*hy-y0)*((k+1/2)*hz-z0)
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "DuDy"
     is similar to existing symbol "DuDx".
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "DuDz"
     is similar to existing symbols {DuDx, DuDy}.
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "D2uDx2"
     is similar to existing symbol "D2uDy2".
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "D2uDz2"
     is similar to existing symbols {D2uDx2, D2uDy2}.
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "D2uDxDz"
     is similar to existing symbol "D2uDxDy".
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell1
     will be suppressed during this calculation.
:[font = output; output; inactive; preserveAspect; endGroup]
u[i_, j_, k_] -> 
 
  (hx*(1/2 + i) - x0)*DuDx[x0, y0, z0] + 
 
   (hy*(1/2 + j) - y0)*DuDy[x0, y0, z0] + 
 
   (hz*(1/2 + k) - z0)*DuDz[x0, y0, z0] + 
 
   (hx*(1/2 + i) - x0)*(hy*(1/2 + j) - y0)*
 
    D2uDxDy[x0, y0, z0] + 
 
   (hx*(1/2 + i) - x0)*(hz*(1/2 + k) - z0)*
 
    D2uDxDz[x0, y0, z0] + 
 
   ((hy*(1/2 + j) - y0)^2*D2uDx2[x0, y0, z0])/2 + 
 
   (hy*(1/2 + j) - y0)*(hz*(1/2 + k) - z0)*
 
    D2uDyDz[x0, y0, z0] + 
 
   ((hx*(1/2 + i) - x0)^2*D2uDy2[x0, y0, z0])/2 + 
 
   ((hz*(1/2 + k) - z0)^2*D2uDz2[x0, y0, z0])/2 + 
 
   U[x0, y0, z0]
;[o]
u[i_, j_, k_] -> 
 
       1
  (hx (- + i) - x0) DuDx[x0, y0, z0] + 
       2
 
        1
   (hy (- + j) - y0) DuDy[x0, y0, z0] + 
        2
 
        1
   (hz (- + k) - z0) DuDz[x0, y0, z0] + 
        2
 
        1                 1
   (hx (- + i) - x0) (hy (- + j) - y0) 
        2                 2
 
    D2uDxDy[x0, y0, z0] + 
 
        1                 1
   (hx (- + i) - x0) (hz (- + k) - z0) 
        2                 2
 
    D2uDxDz[x0, y0, z0] + 
 
        1           2
   (hy (- + j) - y0)  D2uDx2[x0, y0, z0]
        2
   ------------------------------------- + 
                     2
 
        1                 1
   (hy (- + j) - y0) (hz (- + k) - z0) 
        2                 2
 
    D2uDyDz[x0, y0, z0] + 
 
        1           2
   (hx (- + i) - x0)  D2uDy2[x0, y0, z0]
        2
   ------------------------------------- + 
                     2
 
        1           2
   (hz (- + k) - z0)  D2uDz2[x0, y0, z0]
        2
   ------------------------------------- + U[x0, y0, z0]
                     2
:[font = input; preserveAspect]
dog = dudy[i+half,j,k] //. {dudydef1,dudydef2,taylorudef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,k->0,x0->hx,y0->hy/2,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx, hy/2, hz/2]
;[o]
         hy  hz
DuDy[hx, --, --]
         2   2
:[font = input; preserveAspect]
dog = dudy[i,j+half,k] //. {dudydef1,dudydef2,taylorudef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,k->0,x0->hx/2,y0->hy,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx/2, hy, hz/2]
;[o]
     hx      hz
DuDy[--, hy, --]
     2       2
:[font = input; preserveAspect; startGroup]
taylorvdef = v[i_,j_,k_]->
			V[x0,y0,z0]+
			DvDx[x0,y0,z0]*((i+1/2)*hx-x0)+
			DvDy[x0,y0,z0]*((j+1/2)*hy-y0)+
			DvDz[x0,y0,z0]*((k+1/2)*hz-z0)+
			D2vDy2[x0,y0,z0]/2*((i+1/2)*hx-x0)^2+
			D2vDx2[x0,y0,z0]/2*((j+1/2)*hy-y0)^2+
			D2vDz2[x0,y0,z0]/2*((k+1/2)*hz-z0)^2+
			D2vDxDy[x0,y0,z0]*((i+1/2)*hx-x0)*((j+1/2)*hy-y0)+
			D2vDxDz[x0,y0,z0]*((i+1/2)*hx-x0)*((k+1/2)*hz-z0)+
			D2vDyDz[x0,y0,z0]*((j+1/2)*hy-y0)*((k+1/2)*hz-z0)
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
General::spell: 
   Possible spelling error: new symbol name "DvDz"
     is similar to existing symbols {DuDz, DvDx, DvDy}.
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
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell
     will be suppressed during this calculation.
:[font = output; output; inactive; preserveAspect; endGroup]
v[i_, j_, k_] -> 
 
  (hx*(1/2 + i) - x0)*DvDx[x0, y0, z0] + 
 
   (hy*(1/2 + j) - y0)*DvDy[x0, y0, z0] + 
 
   (hz*(1/2 + k) - z0)*DvDz[x0, y0, z0] + 
 
   (hx*(1/2 + i) - x0)*(hy*(1/2 + j) - y0)*
 
    D2vDxDy[x0, y0, z0] + 
 
   (hx*(1/2 + i) - x0)*(hz*(1/2 + k) - z0)*
 
    D2vDxDz[x0, y0, z0] + 
 
   ((hy*(1/2 + j) - y0)^2*D2vDx2[x0, y0, z0])/2 + 
 
   (hy*(1/2 + j) - y0)*(hz*(1/2 + k) - z0)*
 
    D2vDyDz[x0, y0, z0] + 
 
   ((hx*(1/2 + i) - x0)^2*D2vDy2[x0, y0, z0])/2 + 
 
   ((hz*(1/2 + k) - z0)^2*D2vDz2[x0, y0, z0])/2 + 
 
   V[x0, y0, z0]
;[o]
v[i_, j_, k_] -> 
 
       1
  (hx (- + i) - x0) DvDx[x0, y0, z0] + 
       2
 
        1
   (hy (- + j) - y0) DvDy[x0, y0, z0] + 
        2
 
        1
   (hz (- + k) - z0) DvDz[x0, y0, z0] + 
        2
 
        1                 1
   (hx (- + i) - x0) (hy (- + j) - y0) 
        2                 2
 
    D2vDxDy[x0, y0, z0] + 
 
        1                 1
   (hx (- + i) - x0) (hz (- + k) - z0) 
        2                 2
 
    D2vDxDz[x0, y0, z0] + 
 
        1           2
   (hy (- + j) - y0)  D2vDx2[x0, y0, z0]
        2
   ------------------------------------- + 
                     2
 
        1                 1
   (hy (- + j) - y0) (hz (- + k) - z0) 
        2                 2
 
    D2vDyDz[x0, y0, z0] + 
 
        1           2
   (hx (- + i) - x0)  D2vDy2[x0, y0, z0]
        2
   ------------------------------------- + 
                     2
 
        1           2
   (hz (- + k) - z0)  D2vDz2[x0, y0, z0]
        2
   ------------------------------------- + V[x0, y0, z0]
                     2
:[font = input; preserveAspect]
dog = dvdy[i,j+half,k] //. {dvdydef, taylorvdef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,k->0,x0->hx/2,y0->hy,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DvDy[hx/2, hy, hz/2]
;[o]
     hx      hz
DvDy[--, hy, --]
     2       2
:[font = input; preserveAspect; startGroup]
taylorwdef = w[i_,j_,k_]->
			W[x0,y0,z0]+
			DwDx[x0,y0,z0]*((i+1/2)*hx-x0)+
			DwDy[x0,y0,z0]*((j+1/2)*hy-y0)+
			DwDz[x0,y0,z0]*((k+1/2)*hz-z0)+
			D2wDy2[x0,y0,z0]/2*((i+1/2)*hx-x0)^2+
			D2wDx2[x0,y0,z0]/2*((j+1/2)*hy-y0)^2+
			D2wDz2[x0,y0,z0]/2*((k+1/2)*hz-z0)^2+
			D2wDxDy[x0,y0,z0]*((i+1/2)*hx-x0)*((j+1/2)*hy-y0)+
			D2wDxDz[x0,y0,z0]*((i+1/2)*hx-x0)*((k+1/2)*hz-z0)+
			D2wDyDz[x0,y0,z0]*((j+1/2)*hy-y0)*((k+1/2)*hz-z0)
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "taylorwdef"
     is similar to existing symbols {taylorudef, 
     taylorvdef}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "DwDx"
     is similar to existing symbols {DuDx, DvDx}.
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "DwDy"
     is similar to existing symbols {DuDy, DvDy, DwDx}.
:[font = message; inactive; preserveAspect]
General::stop: 
   Further output of General::spell
     will be suppressed during this calculation.
:[font = output; output; inactive; preserveAspect; endGroup]
w[i_, j_, k_] -> 
 
  (hx*(1/2 + i) - x0)*DwDx[x0, y0, z0] + 
 
   (hy*(1/2 + j) - y0)*DwDy[x0, y0, z0] + 
 
   (hz*(1/2 + k) - z0)*DwDz[x0, y0, z0] + 
 
   (hx*(1/2 + i) - x0)*(hy*(1/2 + j) - y0)*
 
    D2wDxDy[x0, y0, z0] + 
 
   (hx*(1/2 + i) - x0)*(hz*(1/2 + k) - z0)*
 
    D2wDxDz[x0, y0, z0] + 
 
   ((hy*(1/2 + j) - y0)^2*D2wDx2[x0, y0, z0])/2 + 
 
   (hy*(1/2 + j) - y0)*(hz*(1/2 + k) - z0)*
 
    D2wDyDz[x0, y0, z0] + 
 
   ((hx*(1/2 + i) - x0)^2*D2wDy2[x0, y0, z0])/2 + 
 
   ((hz*(1/2 + k) - z0)^2*D2wDz2[x0, y0, z0])/2 + 
 
   W[x0, y0, z0]
;[o]
w[i_, j_, k_] -> 
 
       1
  (hx (- + i) - x0) DwDx[x0, y0, z0] + 
       2
 
        1
   (hy (- + j) - y0) DwDy[x0, y0, z0] + 
        2
 
        1
   (hz (- + k) - z0) DwDz[x0, y0, z0] + 
        2
 
        1                 1
   (hx (- + i) - x0) (hy (- + j) - y0) 
        2                 2
 
    D2wDxDy[x0, y0, z0] + 
 
        1                 1
   (hx (- + i) - x0) (hz (- + k) - z0) 
        2                 2
 
    D2wDxDz[x0, y0, z0] + 
 
        1           2
   (hy (- + j) - y0)  D2wDx2[x0, y0, z0]
        2
   ------------------------------------- + 
                     2
 
        1                 1
   (hy (- + j) - y0) (hz (- + k) - z0) 
        2                 2
 
    D2wDyDz[x0, y0, z0] + 
 
        1           2
   (hx (- + i) - x0)  D2wDy2[x0, y0, z0]
        2
   ------------------------------------- + 
                     2
 
        1           2
   (hz (- + k) - z0)  D2wDz2[x0, y0, z0]
        2
   ------------------------------------- + W[x0, y0, z0]
                     2
:[font = input; preserveAspect]
dog = dwdy[i,j+half,k] //. {dwdydef1,dwdydef2, taylorwdef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,k->0,x0->hx/2,y0->hy,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
DwDy[hx/2, hy, hz/2]
;[o]
     hx      hz
DwDy[--, hy, --]
     2       2
:[font = input; preserveAspect]
dog = dwdy[i,j,k+half] //. {dwdydef1,dwdydef2, taylorwdef} ;
:[font = input; preserveAspect; startGroup]
dog /. {i->0,j->0,k->0,x0->hx/2,y0->hy/2,z0->hz}
:[font = output; output; inactive; preserveAspect; endGroup]
DwDy[hx/2, hy/2, hz]
;[o]
     hx  hy
DwDy[--, --, hz]
     2   2
:[font = input; preserveAspect; startGroup]
taylormudef = mu[i_,j_]->
			MU[x0,y0,z0]+
			DmuDx[x0,y0,z0]*((i+1/2)*hx-x0)+
			DmuDy[x0,y0,z0]*((j+1/2)*hy-y0)+
			DmuDz[x0,y0,z0]*((k+1/2)*hz-z0)

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
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "DmuDz"
     is similar to existing symbols {DmuDx, DmuDy, DuDz}.
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
mu[i_, j_] -> 
 
  (hx*(1/2 + i) - x0)*DmuDx[x0, y0, z0] + 
 
   (hy*(1/2 + j) - y0)*DmuDy[x0, y0, z0] + 
 
   (hz*(1/2 + k) - z0)*DmuDz[x0, y0, z0] + MU[x0, y0, z0]
;[o]
mu[i_, j_] -> 
 
       1
  (hx (- + i) - x0) DmuDx[x0, y0, z0] + 
       2
 
        1
   (hy (- + j) - y0) DmuDy[x0, y0, z0] + 
        2
 
        1
   (hz (- + k) - z0) DmuDz[x0, y0, z0] + MU[x0, y0, z0]
        2
:[font = section; inactive; preserveAspect; startGroup]
tests
:[font = subsubsection; inactive; preserveAspect; startGroup]
diagonal elements of tau
:[font = input; preserveAspect]
dog = tauxx[i+half,j,k] //.
		{tauxxdef,dudxdef,taylormudef,taylorudef};
:[font = input; preserveAspect; startGroup]
 dog //. {half->1/2,i->0,j->0,k->0,x0->hx,y0 -> hy/2,z0->hz/2 } 
:[font = output; output; inactive; preserveAspect; endGroup]
2*DuDx[hx, hy/2, hz/2]*mu[1/2, 0, 0]
;[o]
           hy  hz     1
2 DuDx[hx, --, --] mu[-, 0, 0]
           2   2      2
:[font = input; preserveAspect]
dog = tauyy[i,j+half,k] //. 
		{tauyydef,dvdydef,taylormudef,taylorvdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx/2,y0->hy,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
2*DvDy[hx/2, hy, hz/2]*mu[0, 1/2, 0]
;[o]
       hx      hz        1
2 DvDy[--, hy, --] mu[0, -, 0]
       2       2         2
:[font = input; preserveAspect]
dog = tauzz[i,j,k+half] //. 
		{tauzzdef,dwdzdef,taylormudef,taylorvdef,taylorwdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx/2,y0->hy/2,z0->hz}
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
2*DwDz[hx/2, hy/2, hz]*mu[0, 0, 1/2]
;[o]
       hx  hy               1
2 DwDz[--, --, hz] mu[0, 0, -]
       2   2                2
:[font = subsubsection; inactive; preserveAspect; startGroup]
 tauxy
:[font = input; preserveAspect]
dog = tauxy[i,j+half,k] //.
		{tauxydef,dudydef1,dudydef2,dvdxdef1,dvdxdef2,
		 taylormudef,taylorudef,taylorvdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx/2,y0->hy,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
(DuDy[hx/2, hy, hz/2] + DvDx[hx/2, hy, hz/2])*mu[0, 1/2, 0]
;[o]
      hx      hz         hx      hz         1
(DuDy[--, hy, --] + DvDx[--, hy, --]) mu[0, -, 0]
      2       2          2       2          2
:[font = input; preserveAspect]
dog = tauxy[i+half,j,k] //.
		{tauxydef,dudydef1,dudydef2,dvdxdef1,dvdxdef2,
		 taylormudef,taylorudef,taylorvdef};
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx,y0->hy/2,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
(DuDy[hx, hy/2, hz/2] + DvDx[hx, hy/2, hz/2])*mu[1/2, 0, 0]
;[o]
          hy  hz             hy  hz      1
(DuDy[hx, --, --] + DvDx[hx, --, --]) mu[-, 0, 0]
          2   2              2   2       2
:[font = subsubsection; inactive; preserveAspect; startGroup]
 tauxz
:[font = input; preserveAspect]
dog = tauxz[i+half,j,k] //.
		{tauxzdef,dudzdef1,dudzdef2,dwdxdef1,dwdxdef2,
		 taylormudef,taylorudef,taylorwdef} ;
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx,y0->hy/2,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
(DuDz[hx, hy/2, hz/2] + DwDx[hx, hy/2, hz/2])*mu[1/2, 0, 0]
;[o]
          hy  hz             hy  hz      1
(DuDz[hx, --, --] + DwDx[hx, --, --]) mu[-, 0, 0]
          2   2              2   2       2
:[font = input; preserveAspect]
dog = tauxz[i,j,k+half] //.
		{tauxzdef,dudzdef1,dudzdef2,dwdxdef1,dwdxdef2,
		 taylormudef,taylorudef,taylorwdef} ;
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx/2,y0->hy/2,z0->hz}
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
(DuDz[hx/2, hy/2, hz] + DwDx[hx/2, hy/2, hz])*mu[0, 0, 1/2]
;[o]
      hx  hy             hx  hy                1
(DuDz[--, --, hz] + DwDx[--, --, hz]) mu[0, 0, -]
      2   2              2   2                 2
:[font = subsubsection; inactive; preserveAspect; startGroup]
 tauyz
:[font = input; preserveAspect]
dog = tauyz[i,j+half,k] //.
		{tauyzdef,dvdzdef1,dvdzdef2,dwdydef1,dwdydef2,
		 taylormudef,taylorvdef,taylorwdef} ;
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx/2,y0->hy,z0->hz/2}
:[font = output; output; inactive; preserveAspect; endGroup]
(DvDz[hx/2, hy, hz/2] + DwDy[hx/2, hy, hz/2])*mu[0, 1/2, 0]
;[o]
      hx      hz         hx      hz         1
(DvDz[--, hy, --] + DwDy[--, hy, --]) mu[0, -, 0]
      2       2          2       2          2
:[font = input; preserveAspect]
dog = tauyz[i,j,k+half] //.
		{tauyzdef,dvdzdef1,dvdzdef2,dwdydef1,dwdydef2,
		 taylormudef,taylorvdef,taylorwdef} ;
:[font = input; preserveAspect; startGroup]
dog //. {half->1/2,i->0,j->0,k->0,x0->hx/2,y0->hy/2,z0->hz}
:[font = output; output; inactive; preserveAspect; endGroup; endGroup; endGroup]
(DvDz[hx/2, hy/2, hz] + DwDy[hx/2, hy/2, hz])*mu[0, 0, 1/2]
;[o]
      hx  hy             hx  hy                1
(DvDz[--, --, hz] + DwDy[--, --, hz]) mu[0, 0, -]
      2   2              2   2                 2
:[font = section; inactive; preserveAspect; startGroup]
definitions used for fortran output
:[font = input; preserveAspect; startGroup]
murepl1 = mu[i_,j_+half,k_] -> muY[i,j+1,k]
:[font = output; output; inactive; preserveAspect; endGroup]
mu[i_, half + (j_), k_] -> muY[i, 1 + j, k]
;[o]
mu[i_, half + (j_), k_] -> muY[i, 1 + j, k]
:[font = input; preserveAspect; startGroup]
murepl2 = mu[i_+half,j_,k_] -> muX[i+1,j,k]
:[font = output; output; inactive; preserveAspect; endGroup]
mu[half + (i_), j_, k_] -> muX[1 + i, j, k]
;[o]
mu[half + (i_), j_, k_] -> muX[1 + i, j, k]
:[font = input; preserveAspect; startGroup]
murepl3 = mu[i_,j_,k_+half] -> muZ[i,j,k+1]
:[font = output; output; inactive; preserveAspect; endGroup]
mu[i_, j_, half + (k_)] -> muZ[i, j, 1 + k]
;[o]
mu[i_, j_, half + (k_)] -> muZ[i, j, 1 + k]
:[font = input; preserveAspect; startGroup]
urepl = u[i_,j_,k_] -> U[i,j,k,1]
:[font = output; output; inactive; preserveAspect; endGroup]
u[i_, j_, k_] -> U[i, j, k, 1]
;[o]
u[i_, j_, k_] -> U[i, j, k, 1]
:[font = input; preserveAspect; startGroup]
vrepl = v[i_,j_,k_] -> U[i,j,k,2]
:[font = message; inactive; preserveAspect]
General::spell1: 
   Possible spelling error: new symbol name "vrepl"
     is similar to existing symbol "urepl".
:[font = output; output; inactive; preserveAspect; endGroup]
v[i_, j_, k_] -> U[i, j, k, 2]
;[o]
v[i_, j_, k_] -> U[i, j, k, 2]
:[font = input; preserveAspect; startGroup]
wrepl = w[i_,j_,k_] -> U[i,j,k,3]
:[font = message; inactive; preserveAspect]
General::spell: 
   Possible spelling error: new symbol name "wrepl"
     is similar to existing symbols {urepl, vrepl}.
:[font = output; output; inactive; preserveAspect; endGroup]
w[i_, j_, k_] -> U[i, j, k, 3]
;[o]
w[i_, j_, k_] -> U[i, j, k, 3]
:[font = subsection; inactive; preserveAspect; startGroup]
dependentCellsNotCovered is a function which returns a logical expression suitable for inclusion in fortran.  Give an expression, exp, we wish to determine which mesh locations are accessed by the expression.  However, we do not wish to examine all possible locations, only those outside the grid patch region.  So we provide a second argument, which is a boolean function taking two arguments.  The combination will give logical expressions testing the mask for cells utilized by the
expression and for which the boolean function, logfunc[il,jl], evaluates as true. The third argument is the name of the mask array
:[font = input; preserveAspect]
Clear[ dependentCellsNotCovered ]
:[font = input; preserveAspect; startGroup]
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
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "lexp"
     is similar to existing symbol "exp".
:[font = input; preserveAspect; startGroup]
dependentCellsNotCovered[abba*U[i+1,j-1,k,1]+U[i-1,j-1,k+1,2]
							, Function[{i,j,k},(k>0)] , 
							  masks  ]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
masks[-1 + i, -1 + j, 1 + k] > 0
;[o]
masks[-1 + i, -1 + j, 1 + k] > 0
:[font = subsection; inactive; preserveAspect; startGroup]
dependentCellsCovered is the logical inverse of dependentCellsNotCovered
:[font = input; preserveAspect]
Clear[ dependentCellsCovered ]
:[font = input; preserveAspect]
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
:[font = input; preserveAspect; startGroup]
dependentCellsCovered[abba*U[i+1,j-1,k,1]+U[i-1,j-1,k+1,2]
							, Function[{i,j,k},(k>0)] , 
							  masks  ]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
masks[-1 + i, -1 + j, 1 + k] == 0
;[o]
masks[-1 + i, -1 + j, 1 + k] == 0
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for two sided derivs
:[font = input; preserveAspect]
DTwoX[u_,i_,j_,k_] := (u[i+1,j,k]-u[i-1,j,k])/(2*hx)
:[font = input; preserveAspect]
DTwoX[u_,i_,j_,k_,n_] := (u[i+1,j,k,n]-u[i-1,j,k,n])/(2*hx)
:[font = input; preserveAspect; startGroup]
Simplify[ DTwoX[u,0,0,0] //. {taylorudef,x0->hx/2,y0->hy/2,
									z0->hz/2} ]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDx[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDx[--, --, --]
     2   2   2
:[font = input; preserveAspect; startGroup]
DTwoY[u_,i_,j_,k_] := (u[i,j+1,k]-u[i,j-1,k])/(2*hy)
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "DTwoY"
     is similar to existing symbol "DTwoX".
:[font = input; preserveAspect]
DTwoY[u_,i_,j_,k_,n_] := (u[i,j+1,k,n]-u[i,j-1,k,n])/(2*hy)
:[font = input; preserveAspect; startGroup]
Simplify[ DTwoY[u,0,0,0] //. {taylorudef,x0->hx/2,y0->hy/2,
									z0->hz/2} ]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDy[--, --, --]
     2   2   2
:[font = input; preserveAspect; startGroup]
DTwoZ[u_,i_,j_,k_] := (u[i,j,k+1]-u[i,j,k-1])/(2*hz)
:[font = message; inactive; preserveAspect; endGroup]
General::spell: 
   Possible spelling error: new symbol name "DTwoZ"
     is similar to existing symbols {DTwoX, DTwoY}.
:[font = input; preserveAspect]
DTwoZ[u_,i_,j_,k_,n_] := (u[i,j,k+1,n]-u[i,j,k-1,n])/(2*hz)
:[font = input; preserveAspect; startGroup]
Simplify[ DTwoZ[u,0,0,0] //. {taylorudef,x0->hx/2,y0->hy/2,
									z0->hz/2} ]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
DuDz[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDz[--, --, --]
     2   2   2
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for Do One-sided Derivative in X direction.  if sign is positive, 
it means extend the stencil in the positivie x direction.  if negative, extend
in other direction
:[font = input; preserveAspect]
DOneX[u_,i_,j_,k_,sign_] := 
		(-u[i+2,j,k]+4*u[i+1,j,k]-3*u[i,j,k])/
									(2*hx)  /; sign==1
:[font = input; preserveAspect]
DOneX[u_,i_,j_,k_,sign_] := 
		(u[i-2,j,k]-4*u[i-1,j,k]+3*u[i,j,k])/
									(2*hx)  /; sign==-1
:[font = input; preserveAspect]
DOneX[u_,i_,j_,k_,n_,sign_] := 
		(-u[i+2,j,k,n]+4*u[i+1,j,k,n]-3*u[i,j,k,n])/
									(2*hx)  /; sign==1
:[font = input; preserveAspect]
DOneX[u_,i_,j_,k_,n_,sign_] := 
		(u[i-2,j,k,n]-4*u[i-1,j,k,n]+3*u[i,j,k,n])/
									(2*hx)  /; sign==-1
:[font = input; preserveAspect; startGroup]
Simplify[ DOneX[u,0,0,0,-1] //. 
		{taylorudef,x0->hx/2,y0->hy/2,z0->hz/2}]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDx[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDx[--, --, --]
     2   2   2
:[font = input; preserveAspect; startGroup]
Simplify[ DOneX[u,0,0,0,+1] //. 
		{taylorudef,x0->hx/2,y0->hy/2,z0->hz/2}]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
DuDx[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDx[--, --, --]
     2   2   2
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for Do One-sided Derivative in Y direction.  if sign is positive, 
it means extend the stencil in the positivie y direction.  if negative, extend
in other direction
:[font = input; preserveAspect; startGroup]
DOneY[u_,i_,j_,k_,sign_] := 
		(-u[i,j+2,k]+4*u[i,j+1,k]-3*u[i,j,k])/
									(2*hy)  /; sign==1
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "DOneY"
     is similar to existing symbol "DOneX".
:[font = input; preserveAspect]
DOneY[u_,i_,j_,k_,sign_] := 
		(u[i,j-2,k]-4*u[i,j-1,k]+3*u[i,j,k])/
									(2*hy)  /; sign==-1
:[font = input; preserveAspect]
DOneY[u_,i_,j_,k_,n_,sign_] := 
		(-u[i,j+2,k,n]+4*u[i,j+1,k,n]-3*u[i,j,k,n])/
									(2*hy)  /; sign==1
:[font = input; preserveAspect]
DOneY[u_,i_,j_,k_,n_,sign_] := 
		(u[i,j-2,k,n]-4*u[i,j-1,k,n]+3*u[i,j,k,n])/
									(2*hy)  /; sign==-1
:[font = input; preserveAspect; startGroup]
Simplify[ DOneY[u,0,0,0,-1] //. 
		{taylorudef,x0->hx/2,y0->hy/2,z0->hz/2}]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDy[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDy[--, --, --]
     2   2   2
:[font = input; preserveAspect; startGroup]
Simplify[ DOneY[u,0,0,0,+1] //. 
		{taylorudef,x0->hx/2,y0->hy/2,z0->hz/2}]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup]
DuDy[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDy[--, --, --]
     2   2   2
:[font = subsection; inactive; preserveAspect; startGroup]
definitions for Do One-sided Derivative in Z direction.  if sign is positive, 
it means extend the stencil in the positivie z direction.  if negative, extend
in other direction
:[font = input; preserveAspect; startGroup]
DOneZ[u_,i_,j_,k_,sign_] := 
		(-u[i,j,k+2]+4*u[i,j,k+1]-3*u[i,j,k])/
									(2*hz)  /; sign==1
:[font = message; inactive; preserveAspect; endGroup]
General::spell: 
   Possible spelling error: new symbol name "DOneZ"
     is similar to existing symbols {DOneX, DOneY}.
:[font = input; preserveAspect]
DOneZ[u_,i_,j_,k_,sign_] := 
		(u[i,j,k-2]-4*u[i,j,k-1]+3*u[i,j,k])/
									(2*hz)  /; sign==-1
:[font = input; preserveAspect]
DOneZ[u_,i_,j_,k_,n_,sign_] := 
		(-u[i,j,k+2,n]+4*u[i,j,k+1,n]-3*u[i,j,k,n])/
									(2*hz)  /; sign==1
:[font = input; preserveAspect]
DOneZ[u_,i_,j_,k_,n_,sign_] := 
		(u[i,j,k-2,n]-4*u[i,j,k-1,n]+3*u[i,j,k,n])/
									(2*hz)  /; sign==-1
:[font = input; preserveAspect; startGroup]
Simplify[ DOneZ[u,0,0,0,-1] //. 
		{taylorudef,x0->hx/2,y0->hy/2,z0->hz/2}]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDz[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDz[--, --, --]
     2   2   2
:[font = input; preserveAspect; startGroup]
Simplify[ DOneZ[u,0,0,0,+1] //. 
		{taylorudef,x0->hx/2,y0->hy/2,z0->hz/2}]
:[font = output; output; inactive; preserveAspect; endGroup]
DuDz[hx/2, hy/2, hz/2]
;[o]
     hx  hy  hz
DuDz[--, --, --]
     2   2   2
:[font = subsubsection; inactive; preserveAspect; startGroup]
useful one-sided derivatives
:[font = input; preserveAspect]
dvdxalt[i_,j_+half,k_,sign_] := (DOneX[v,i,j  ,k,sign]+
							     DOneX[v,i,j+1,k,sign])/2
:[font = input; preserveAspect]
dudyalt[i_+half,j_,k_,sign_] := (DOneY[u,i  ,j,k,sign]+
							     DOneY[u,i+1,j,k,sign])/2
:[font = input; preserveAspect; startGroup]
dvdzalt[i_,j_+half,k_,sign_] := (DOneZ[v,i,j  ,k,sign]+
							     DOneZ[v,i,j+1,k,sign])/2
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "dvdzalt"
     is similar to existing symbol "dvdxalt".
:[font = input; preserveAspect; startGroup]
dwdyalt[i_,j_,k_+half,sign_] := (DOneY[w,i,j,k  ,sign]+
								 DOneY[w,i,j,k+1,sign])/2
:[font = message; inactive; preserveAspect; endGroup]
General::spell1: 
   Possible spelling error: new symbol name "dwdyalt"
     is similar to existing symbol "dudyalt".
:[font = input; preserveAspect; startGroup]
dudzalt[i_+half,j_,k_,sign_] := (DOneZ[u,i  ,j,k,sign]+
								 DOneZ[u,i+1,j,k,sign])/2
:[font = message; inactive; preserveAspect; endGroup]
General::spell: 
   Possible spelling error: new symbol name "dudzalt"
     is similar to existing symbols {dudyalt, dvdzalt}.
:[font = input; preserveAspect; startGroup]
dwdxalt[i_,j_,k_+half,sign_] := (DOneX[w,i,j,k  ,sign]+
								 DOneX[w,i,j,k+1,sign])/2
:[font = message; inactive; preserveAspect; endGroup; endGroup; endGroup]
General::spell: 
   Possible spelling error: new symbol name "dwdxalt"
     is similar to existing symbols {dvdxalt, dwdyalt}.
:[font = subsection; inactive; preserveAspect; startGroup]
setup to use the Format.m package 
:[font = input; preserveAspect]
Off[General::spell,General::spell1];
:[font = input; preserveAspect]
SetOptions[$Output,PageWidth->73];
:[font = input; preserveAspect; startGroup]
<</usr/people/wyc/_math/MathSource/Format/Format.m
:[font = message; inactive; preserveAspect]
exp::shdw: Warning: Symbol exp
     appears in multiple contexts {Format`, Global`};
     definitions in context Format`
     may shadow or be shadowed by other definitions.
:[font = message; inactive; preserveAspect; endGroup]
sign::shdw: 
   Warning: Symbol sign appears in multiple contexts 
    {Format`, Global`}; definitions in context Format`
     may shadow or be shadowed by other definitions.
:[font = subsubsection; inactive; preserveAspect; startGroup]
substitution for all derivatives and variables
:[font = input; preserveAspect]
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
:[font = input; preserveAspect; endGroup]
allUVW = {       urepl,
				 vrepl,
				 wrepl};
:[font = subsubsection; inactive; preserveAspect; startGroup]
transverse u derivs
:[font = input; preserveAspect]
tduext[half+i,j,k] = trandere[i+1,j,k,1];
:[font = input; preserveAspect]
tduext[half+i-1,j,k] = tranderw[i-1,j,k,1];
:[font = input; preserveAspect]
tduext[i,j+half,k] = trandern[i,j+1,k,1];
:[font = input; preserveAspect]
tduext[i,j-1+half,k] = tranders[i,j-1,k,1];
:[font = input; preserveAspect]
tduext[i,j,k+half] = trandert[i,j,k+1,1];
:[font = input; preserveAspect; endGroup]
tduext[i,j,k-1+half] = tranderb[i,j,k-1,1];
:[font = subsubsection; inactive; preserveAspect; startGroup]
transverse v derivs
:[font = input; preserveAspect]
tdvext[half+i,j,k] = trandere[i+1,j,k,2];
:[font = input; preserveAspect]
tdvext[half+i-1,j,k] = tranderw[i-1,j,k,2];
:[font = input; preserveAspect]
tdvext[i,j+half,k] = trandern[i,j+1,k,2];
:[font = input; preserveAspect]
tdvext[i,j-1+half,k] = tranders[i,j-1,k,2];
:[font = input; preserveAspect]
tdvext[i,j,k+half] = trandert[i,j,k+1,2];
:[font = input; preserveAspect; endGroup]
tdvext[i,j,k-1+half] = tranderb[i,j,k-1,2];
:[font = subsubsection; inactive; preserveAspect; startGroup]
transverse w derivs
:[font = input; preserveAspect]
tdwext[half+i,j,k] = trandere[i+1,j,k,3];
:[font = input; preserveAspect]
tdwext[half+i-1,j,k] = tranderw[i-1,j,k,3];
:[font = input; preserveAspect]
tdwext[i,j+half,k] = trandern[i,j+1,k,3];
:[font = input; preserveAspect]
tdwext[i,j-1+half,k] = tranders[i,j-1,k,3];
:[font = input; preserveAspect]
tdwext[i,j,k+half] = trandert[i,j,k+1,3];
:[font = input; preserveAspect; endGroup; endGroup]
tdwext[i,j,k-1+half] = tranderb[i,j,k-1,3];
:[font = subsection; inactive; preserveAspect; startGroup]
an alternate approach which seeks to automatically determine which 
direction to use for one sided deriv
:[font = input; preserveAspect]
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
				AssignToArray->{U}
			],
			
			tmpcond =!= False && sign != 0,
			(* exp extends outside, output conditional mask *)

			tmp = FortranAssign[
				tmpcond,
				AssignToArray->{mask},
				AssignIndent->""
			];
			tmpalt =dependentCellsNotCovered[
					expalt[indx,indx,indz,sign]//.allDerivAllUV,
						indexcond,mask];
			line1 = StringForm["      if(``) then ", tmp];
			line2 = FortranAssign[
				lhs,
				expalt[indx,indy,indz,sign]//.allDerivAllUV,
				AssignToArray->{U}
			];
			line3 = StringForm["      else"];
			line4 = FortranAssign[
				lhs,
				exp[indx,indy,indz] //.allDerivAllUV,
				AssignToArray->{U}
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
									trandert,tranderb}]	
		]	
	]
:[font = input; preserveAspect; startGroup]
altgen[ dwdxt, 
		i,j,k+half,
		dwdx,dwdxalt,
		trandern,3,1,Function[{i,j,k},(k>0)],
		maskn ]
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
        dwdxt=trandern(i,j,1+k,3,1)
;[o]
        dwdxt=trandern(i,j,1+k,3,1)
:[font = input; preserveAspect; startGroup]
altgen[ dudye, 
		i+half,j,k,
		dudy,
		dudyalt,
		tduext,
		1,2,
		Function[{i,j,k},(i<0)],
		maskn ] 
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
;[o]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
:[font = input; preserveAspect; startGroup]
altgen[ dudye, 
		i+half,j,k,
		dudy,
		dudyalt,
		trandere,
		1,2,
		Function[{i,j,k},(i>0)],
		maskn ] 
:[font = output; output; inactive; preserveAspect; endGroup]
        dudye=trandere(1+i,j,k,1,2)
;[o]
        dudye=trandere(1+i,j,k,1,2)
:[font = input; preserveAspect; startGroup]
altgen[ dudye, 
		i+half,j,k,
		dudy,
		dudyalt,
		tranderns,
		1,2,
		Function[{i,j,k},(j<0)],
		maskn ] 
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
StringForm["      if(``) then ", 
 
  (maskn(i,-1+j,k).gt.0).or.(maskn(1+i,-1+j,k).gt.0)]
        dudye=5.d-1*(5.d-1*(-3.d0*U(i,j,k,1)+4.d0*U(i,1+j,
     &  k,1)-U(i,2+j,k,1))/hy+5.d-1*(-3.d0*U(1+i,j,k,1)+4.
     &  d0*U(1+i,1+j,k,1)-U(1+i,2+j,k,1))/hy)
StringForm["      else"]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
StringForm["      endif"]
;[o]
      if((maskn(i,-1+j,k).gt.0).or.(maskn(1+i,-1+j,k).gt.0)
 
  ) then 
        dudye=5.d-1*(5.d-1*(-3.d0*U(i,j,k,1)+4.d0*U(i,1+j,
     &  k,1)-U(i,2+j,k,1))/hy+5.d-1*(-3.d0*U(1+i,j,k,1)+4.
     &  d0*U(1+i,1+j,k,1)-U(1+i,2+j,k,1))/hy)
      else
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
      endif
:[font = input; preserveAspect; startGroup]
altgen[ dudye, 
		i+half,j,k,
		dudy,
		dudyalt,
		trandere,
		1,2,
		Function[{i,j,k},(j>0)],
		maskn ] 
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
StringForm["      if(``) then ", 
 
  (maskn(i,1+j,k).gt.0).or.(maskn(1+i,1+j,k).gt.0)]
        dudye=5.d-1*(5.d-1*(U(i,-2+j,k,1)-4.d0*U(i,-1+j,k,
     &  1)+3.d0*U(i,j,k,1))/hy+5.d-1*(U(1+i,-2+j,k,1)-4.d0
     &  *U(1+i,-1+j,k,1)+3.d0*U(1+i,j,k,1))/hy)
StringForm["      else"]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
StringForm["      endif"]
;[o]
      if((maskn(i,1+j,k).gt.0).or.(maskn(1+i,1+j,k).gt.0)
 
  ) then 
        dudye=5.d-1*(5.d-1*(U(i,-2+j,k,1)-4.d0*U(i,-1+j,k,
     &  1)+3.d0*U(i,j,k,1))/hy+5.d-1*(U(1+i,-2+j,k,1)-4.d0
     &  *U(1+i,-1+j,k,1)+3.d0*U(1+i,j,k,1))/hy)
      else
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
      endif
:[font = input; preserveAspect; startGroup]
altgen[ dudye, 
		i+half,j,k,
		dudy,
		dudyalt,
		tranderb,
		1,2,
		Function[{i,j,k},(k<0)],
		maskn ] 
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
;[o]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
:[font = input; preserveAspect; startGroup]
altgen[ dudye, 
		i+half,j,k,
		dudy,
		dudyalt,
		trandert,
		1,2,
		Function[{i,j,k},(k>0)],
		maskn ] 
:[font = output; output; inactive; preserveAspect; fontLeading = 0; endGroup]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
;[o]
        dudye=2.5d-1*(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1
     &  +j,k,1)+U(1+i,1+j,k,1))/hy
:[font = subsubsection; inactive; preserveAspect; startGroup]
a short-hand function
:[font = input; preserveAspect]
FA[x_] := FortranAssign[x,
						AssignToArray->{U,muX,muY,muZ,a,u,
										maskn,maske,maskw,
										masks,maskt,maskb},
						AssignIndent->"" ];		
:[font = input; preserveAspect]
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D1.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D1.F"];
CopyFile[ "dog.f", "DV_3D1.F" ];
<<"!touch DV_3D1.F"

:[font = input; preserveAspect]
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D2.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D2.F"];
CopyFile[ "dog.f", "DV_3D2.F" ];
<<"!touch DV_3D2.F"

:[font = input; preserveAspect]
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D3.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D3.F"];
CopyFile[ "dog.f", "DV_3D3.F" ];
<<"!touch DV_3D3.F"

:[font = input; preserveAspect]
DeleteFile[ "dog.mf"];
CopyFile[ "DV_3D4.mF" , "dog.mf"];
Splice["dog.mf",FormatType->OutputForm];
DeleteFile[ "DV_3D4.F"];
CopyFile[ "dog.f", "DV_3D4.F" ];
<<"!touch DV_3D4.F"

:[font = input; preserveAspect; startGroup]

 FA[ Coefficient[Expand[
  alpha*a[i,j,k]*u[i,j,k]-beta*(
	hy*hz*(tauxx[i+half  ,j      , k]-
	       tauxx[i-1+half,j      , k])+
	hx*hz*(tauxy[i       ,j+half  ,k]-   
	       tauxy[i       ,j-1+half,k])+
	hx*hy*(tauxz[i       ,j       ,k+half]-
	       tauxz[i       ,j       ,k-1+half]))/vol //. 
	       {murepl1,murepl2,murepl3,
	        tauxxdef,tauyydef,tauxydef,tauxzdef,tauyzdef,
	        tauzzdef,
	        vol->hx*hy*hz,
	        dudxdef,dvdydef,dwdzdef,
	        dudydef1,dudydef2,
	        dudzdef1,dudzdef2,
	        dwdxdef1,dwdxdef2,
	        dvdxdef1,dvdxdef2 }
 ],u[i,j+1,k]]]
	       
:[font = output; output; inactive; preserveAspect; endGroup]
-(beta*muY(i,1+j,k)/hy**2)
;[o]
-(beta*muY(i,1+j,k)/hy**2)
:[font = input; preserveAspect; startGroup]

FA[Coefficient[ Expand[
  alpha*a[i,j,k]*v[i,j,k]-beta*(
	hy*hz*(tauxy[i+half  ,j      , k]-
	       tauxy[i-1+half,j      , k])+
	hx*hz*(tauyy[i       ,j+half  ,k]-   
	       tauyy[i       ,j-1+half,k])+
	hx*hy*(tauyz[i       ,j       ,k+half]-
	       tauyz[i       ,j       ,k-1+half]))/vol //. 
	       {murepl1,murepl2,murepl3,
	        tauxxdef,tauyydef,tauxydef,tauxzdef,tauyzdef,
	        tauzzdef,
	        vol->hx*hy*hz,
	        dudxdef,dvdydef,dwdzdef,
	        dudydef1,dudydef2,
	        dudzdef1,dudzdef2,
	        dwdxdef1,dwdxdef2,
	        dvdxdef1,dvdxdef2,
	        dvdzdef1,dvdzdef2,
	        dwdydef1,dwdydef2 } ] , v[i,j+1,k] ] ]
:[font = output; output; inactive; preserveAspect; endGroup]
-2.d0*beta*muY(i,1+j,k)/hy**2
;[o]
-2.d0*beta*muY(i,1+j,k)/hy**2
:[font = input; preserveAspect; startGroup]

FA[ Coefficient[ Expand[
  alpha*a[i,j,k]*w[i,j,k]-beta*(
	hy*hz*(tauxz[i+half  ,j      , k]-
	       tauxz[i-1+half,j      , k])+
	hx*hz*(tauyz[i       ,j+half  ,k]-   
	       tauyz[i       ,j-1+half,k])+
	hx*hy*(tauzz[i       ,j       ,k+half]-
	       tauzz[i       ,j       ,k-1+half]))/vol //. 
	       {murepl1,murepl2,murepl3,
	        tauxxdef,tauyydef,tauxydef,tauxzdef,tauyzdef,
	        tauzzdef,
	        vol->hx*hy*hz,
	        dudxdef,dvdydef,dwdzdef,
	        dudydef1,dudydef2,
	        dudzdef1,dudzdef2,
	        dwdxdef1,dwdxdef2,
	        dvdxdef1,dvdxdef2,
	        dvdzdef1,dvdzdef2,
	        dwdydef1,dwdydef2 } ] , w[i,j+1,k] ] ]
:[font = output; output; inactive; preserveAspect; endGroup]
-(beta*muY(i,1+j,k)/hy**2)
;[o]
-(beta*muY(i,1+j,k)/hy**2)
:[font = input; preserveAspect; startGroup]

	DTwoY[u,i,j,k] //. allDerivAllUV
:[font = output; output; inactive; preserveAspect; endGroup]
(-U[i, -1 + j, k, 1] + U[i, 1 + j, k, 1])/(2*hy)
;[o]
-U[i, -1 + j, k, 1] + U[i, 1 + j, k, 1]
---------------------------------------
                 2 hy
:[font = input; preserveAspect; startGroup]
FA[ dependentCellsNotCovered[
	DTwoY[u,i,j,k] //. allDerivAllUV,
	Function[{i,j,k},(j<0)],masks] ]
:[font = output; output; inactive; preserveAspect; endGroup]
masks(i,-1+j,k).gt.0
;[o]
masks(i,-1+j,k).gt.0
:[font = input; preserveAspect; startGroup]
FA[ DTwoY[U,i,j,k,n] ]
:[font = output; output; inactive; preserveAspect; endGroup; endGroup; endGroup; endGroup; endGroup]
5.d-1*(-U(i,-1+j,k,n)+U(i,1+j,k,n))/hy
;[o]
5.d-1*(-U(i,-1+j,k,n)+U(i,1+j,k,n))/hy
^*)