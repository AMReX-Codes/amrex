
#include "cache.H"

copy_cache::copy_cache(int Nsets, Real *Dptr, Real *Sptr)
  : nsets(Nsets), dptr(Dptr), sptr(Sptr), bdy_cache(NULL)
{
#if (BL_SPACEDIM == 2)
  dstart = new int[5 * nsets];
  sstart = dstart + nsets;
  dstrid = dstart + 2 * nsets;
  sstrid = dstart + 3 * nsets;
  nvals  = dstart + 4 * nsets;
  for (int i = 0; i < nsets; i++)
    nvals[i] = 0;
#else
  dstart = new int[8 * nsets];
  sstart = dstart + nsets;
  dstrid1 = dstart + 2 * nsets;
  dstrid2 = dstart + 3 * nsets;
  sstrid1 = dstart + 4 * nsets;
  sstrid2 = dstart + 5 * nsets;
  nvals1  = dstart + 6 * nsets;
  nvals2  = dstart + 7 * nsets;
  for (int i = 0; i < nsets; i++) {
    nvals1[i] = 0;
    nvals2[i] = 0;
  }
#endif
}

// sync cache

copy_cache::copy_cache(MultiFab& r, const level_interface& interface,
		       amr_boundary bdy)
{
  assert(r.length() > 0);
  assert(r.nComp() == 1);
  assert(type(r) == IntVect::TheNodeVector());

  int igrid, jgrid, iface, icor, i;

  nsets = 0;
  for (i = 0; i < BL_SPACEDIM; i++) {
    for (igrid = 0; igrid < interface.nboxes(i); igrid++) {
      if (interface.geo(i, igrid) != level_interface::ALL)
	break;
      nsets++;
    }
  }

  // While it is possible for two copies to occur at a corner in 3D,
  // this implies that the edge from that corner in the +x direction
  // exists and does not copy, so the above counts are sufficient.

  Real *baseptr = r[0].dataPtr();
  dptr = sptr = baseptr;

  bdy_cache = NULL;

#if (BL_SPACEDIM == 2)
  dstart = new int[5 * nsets];
  sstart = dstart + nsets;
  dstrid = dstart + 2 * nsets;
  sstrid = dstart + 3 * nsets;
  nvals  = dstart + 4 * nsets;
  for (i = 0; i < nsets; i++)
    nvals[i] = 0;
#else
  dstart = new int[8 * nsets];
  sstart = dstart + nsets;
  dstrid1 = dstart + 2 * nsets;
  dstrid2 = dstart + 3 * nsets;
  sstrid1 = dstart + 4 * nsets;
  sstrid2 = dstart + 5 * nsets;
  nvals1  = dstart + 6 * nsets;
  nvals2  = dstart + 7 * nsets;
  for (i = 0; i < nsets; i++) {
    nvals1[i] = 0;
    nvals2[i] = 0;
  }
#endif

  int iset = 0;
  for (iface = 0; iface < interface.nfaces(); iface++) {
    igrid = interface.fgrid(iface, 0);
    jgrid = interface.fgrid(iface, 1);
    if (igrid < 0 || jgrid < 0 || interface.fgeo(iface) != level_interface::ALL)
      break;
    const Box& b = interface.node_face(iface);
#if (BL_SPACEDIM == 2)
    int dstartj, sstarti, stridi, stridj, nvals;
    stridi = r[igrid].box().length(0);
    stridj = r[jgrid].box().length(0);
    sstarti = r[igrid].dataPtr() - baseptr +
	      b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	      stridi * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
    dstartj = r[jgrid].dataPtr() - baseptr +
	      b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	      stridj * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
    if (interface.fdim(iface) == 0) {
      nvals = b.length(1);
    }
    else {
      nvals = b.length(0);
      stridi = 1;
      stridj = 1;
    }
    set(iset++, dstartj, sstarti, stridj, stridi, nvals);
#else
    int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2, nvals1, nvals2;
    stridi1 = r[igrid].box().length(0);
    stridi2 = stridi1 * r[igrid].box().length(1);
    stridj1 = r[jgrid].box().length(0);
    stridj2 = stridj1 * r[jgrid].box().length(1);
    sstarti = r[igrid].dataPtr() - baseptr +
              b.smallEnd(0) - r[igrid].box().smallEnd(0) +
              stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
	      stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
    dstartj = r[jgrid].dataPtr() - baseptr +
	      b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	      stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
	      stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
    if (interface.fdim(iface) == 0) {
      nvals1 = b.length(1);
      nvals2 = b.length(2);
    }
    else if (interface.fdim(iface) == 1) {
      nvals1 = b.length(0);
      nvals2 = b.length(2);
      stridi1 = 1;
      stridj1 = 1;
    }
    else {
      nvals1 = b.length(0);
      nvals2 = b.length(1);
      stridi2 = stridi1;
      stridj2 = stridj1;
      stridi1 = 1;
      stridj1 = 1;
    }
    set(iset++, dstartj, sstarti, stridj1, stridj2,
	       stridi1, stridi2, nvals1, nvals2);
#endif
  }

#if (BL_SPACEDIM == 2)
  for (icor = 0; icor < interface.ncorners(); icor++) {
    igrid = interface.cgrid(icor, 0);
    jgrid = interface.cgrid(icor, 3);
    // only do interior corners with fine grid on all sides
    if (igrid < 0 || jgrid < 0 || interface.cgeo(icor) != level_interface::ALL)
      break;
    if (jgrid == interface.cgrid(icor, 1)) {
      const Box& b = interface.corner(icor);
      int dstartj, sstarti, stridi, stridj;
      stridi = r[igrid].box().length(0);
      stridj = r[jgrid].box().length(0);
      sstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	        stridi * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	        stridj * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
      set(iset++, dstartj, sstarti, 0, 0, 1);
    }
  }
#else
  int iedge;
  for (iedge = 0; iedge < interface.nedges(); iedge++) {
    igrid = interface.egrid(iedge, 0);
    jgrid = interface.egrid(iedge, 3);
    // only do interior edges with fine grid on all sides
    if (igrid < 0 || jgrid < 0 || interface.egeo(iedge) != level_interface::ALL)
      break;
    if (jgrid == interface.egrid(iedge, 1)) {
      const Box& b = interface.node_edge(iedge);
      int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2, nvals1;
      stridi1 = r[igrid].box().length(0);
      stridi2 = stridi1 * r[igrid].box().length(1);
      stridj1 = r[jgrid].box().length(0);
      stridj2 = stridj1 * r[jgrid].box().length(1);
      sstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	        stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
	        stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	        stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
	        stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
      if ((nvals1 = b.length(0)) > 1) {
	stridi1 = 1;
	stridj1 = 1;
      }
      else if ((nvals1 = b.length(2)) > 1) {
	stridi1 = stridi2;
	stridj1 = stridj2;
      }
      else {
	nvals1 = b.length(1);
      }
      set(iset++, dstartj, sstarti, stridj1, 0,
		 stridi1, 0, nvals1, 1);
    }
  }
  for (icor = 0; icor < interface.ncorners(); icor++) {
    igrid = interface.cgrid(icor, 0);
    jgrid = interface.cgrid(icor, 7);
    // only do interior corners with fine grid on all sides
    if (igrid < 0 || jgrid < 0 || interface.cgeo(icor) != level_interface::ALL)
      break;
    if (interface.cgrid(icor, 3) == interface.cgrid(icor, 1)) {
      if (jgrid != interface.cgrid(icor, 3)) {
	const Box& b = interface.corner(icor);
	int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2;
	stridi1 = r[igrid].box().length(0);
	stridi2 = stridi1 * r[igrid].box().length(1);
	stridj1 = r[jgrid].box().length(0);
	stridj2 = stridj1 * r[jgrid].box().length(1);
	sstarti = r[igrid].dataPtr() - baseptr +
	          b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	          stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
	          stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	dstartj = r[jgrid].dataPtr() - baseptr +
	          b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	          stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
	          stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
	jgrid = interface.cgrid(icor, 5);
	if (jgrid != interface.cgrid(icor, 7)) {
	  stridj1 = r[jgrid].box().length(0);
	  stridj2 = stridj1 * r[jgrid].box().length(1);
	  dstartj = r[jgrid].dataPtr() - baseptr +
	            b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	            stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		    stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	  set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
	}
      }
    }
    else if (interface.cgrid(icor, 5) == interface.cgrid(icor, 1)) {
      if (jgrid != interface.cgrid(icor, 5)) {
	const Box& b = interface.corner(icor);
	int dstartj, sstarti, stridi1, stridi2, stridj1, stridj2;
	stridi1 = r[igrid].box().length(0);
	stridi2 = stridi1 * r[igrid].box().length(1);
	stridj1 = r[jgrid].box().length(0);
	stridj2 = stridj1 * r[jgrid].box().length(1);
	sstarti = r[igrid].dataPtr() - baseptr +
	          b.smallEnd(0) - r[igrid].box().smallEnd(0) +
	          stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
	          stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
	dstartj = r[jgrid].dataPtr() - baseptr +
	          b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	          stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
	          stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
	jgrid = interface.cgrid(icor, 3);
	if (jgrid != interface.cgrid(icor, 7)) {
	  stridj1 = r[jgrid].box().length(0);
	  stridj2 = stridj1 * r[jgrid].box().length(1);
	  dstartj = r[jgrid].dataPtr() - baseptr +
	            b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
	            stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		    stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	  set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
	  if (jgrid == interface.cgrid(icor, 2)) {
	    jgrid = interface.cgrid(icor, 6);
	    if (jgrid != interface.cgrid(icor, 7)) {
	      stridj1 = r[jgrid].box().length(0);
	      stridj2 = stridj1 * r[jgrid].box().length(1);
	      dstartj = r[jgrid].dataPtr() - baseptr +
		        b.smallEnd(0) - r[jgrid].box().smallEnd(0) +
			stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
		        stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
	      set(iset++, dstartj, sstarti, 0, 0, 0, 0, 1, 1);
	    }
	  }
	}
      }
    }
  }
#endif

  bdy.set_sync_cache(this, nsets, iset, r, interface);
  nsets = iset;
}

// border cache

copy_cache::copy_cache(MultiFab& r, const level_interface& interface,
		       amr_boundary bdy, int w)
{
  assert(r.length() > 0);
  assert(r.nComp() == 1);
  assert(type(r) == IntVect::TheNodeVector());

  w = (w < 0 || w > r.nGrow()) ? r.nGrow() : w;

  assert(w == 1);

  int igrid, jgrid, iface;

  nsets = 0;
  for (iface = 0; iface < interface.nfaces(); iface++) {
    if (interface.fgeo(iface) != level_interface::ALL)
      break;
    if (interface.fgrid(iface, 0) >= 0)
      nsets++;
    if (interface.fgrid(iface, 1) >= 0)
      nsets++;
  }

  Real *baseptr = r[0].dataPtr();
  dptr = sptr = baseptr;

  bdy_cache = NULL;

#if (BL_SPACEDIM == 2)
  dstart = new int[5 * nsets];
  sstart = dstart + nsets;
  dstrid = dstart + 2 * nsets;
  sstrid = dstart + 3 * nsets;
  nvals  = dstart + 4 * nsets;
  for (int i = 0; i < nsets; i++)
    nvals[i] = 0;
#else
  dstart = new int[8 * nsets];
  sstart = dstart + nsets;
  dstrid1 = dstart + 2 * nsets;
  dstrid2 = dstart + 3 * nsets;
  sstrid1 = dstart + 4 * nsets;
  sstrid2 = dstart + 5 * nsets;
  nvals1  = dstart + 6 * nsets;
  nvals2  = dstart + 7 * nsets;
  for (int i = 0; i < nsets; i++) {
    nvals1[i] = 0;
    nvals2[i] = 0;
  }
#endif

  int iset = 0;
  for (iface = 0; iface < interface.nfaces(); iface++) {
    igrid = interface.fgrid(iface, 0);
    jgrid = interface.fgrid(iface, 1);
    if (igrid < 0 || jgrid < 0 || interface.fgeo(iface) != level_interface::ALL)
      break;
    const Box& b = interface.node_face(iface);
#if (BL_SPACEDIM == 2)
    int dstarti, dstartj, sstarti, sstartj, stridi, stridj, nvals;
    if (interface.fdim(iface) == 0) {
      nvals = b.length(1);
      stridi = r[igrid].box().length(0);
      stridj = r[jgrid].box().length(0);
      dstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
	        stridi * (b.smallEnd(1) - r[igrid].box().smallEnd(1));
      sstarti = dstarti - 2;
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
	        stridj * (b.smallEnd(1) - r[jgrid].box().smallEnd(1));
      sstartj = dstartj + 2;
    }
    else {
      nvals = b.length(0) + 2 * w;
      stridi = r[igrid].box().length(0);
      stridj = r[jgrid].box().length(0);
      dstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - w - r[igrid].box().smallEnd(0) +
	        stridi * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1));
      sstarti = dstarti - 2 * stridi;
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - w - r[jgrid].box().smallEnd(0) +
	        stridj * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1));
      sstartj = dstartj + 2 * stridj;
      stridi = 1;
      stridj = 1;
    }
    set(iset++, dstarti, sstartj, stridi, stridj, nvals);
    set(iset++, dstartj, sstarti, stridj, stridi, nvals);
#else
    int dstarti, dstartj, sstarti, sstartj;
    int stridi1, stridi2, stridj1, stridj2;
    if (interface.fdim(iface) == 0) {
      int nvals1, nvals2;
      nvals1 = b.length(1);
      nvals2 = b.length(2);
      stridi1 = r[igrid].box().length(0);
      stridi2 = stridi1 * r[igrid].box().length(1);
      stridj1 = r[jgrid].box().length(0);
      stridj2 = stridj1 * r[jgrid].box().length(1);
      dstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) + 1 - r[igrid].box().smallEnd(0) +
	        stridi1 * (b.smallEnd(1) - r[igrid].box().smallEnd(1)) +
	        stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
      sstarti = dstarti - 2;
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - 1 - r[jgrid].box().smallEnd(0) +
	        stridj1 * (b.smallEnd(1) - r[jgrid].box().smallEnd(1)) +
	        stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
      sstartj = dstartj + 2;
      set(iset++, dstarti, sstartj, stridi1, stridi2,
                 stridj1, stridj2, nvals1, nvals2);
      set(iset++, dstartj, sstarti, stridj1, stridj2,
                 stridi1, stridi2, nvals1, nvals2);
    }
    else if (interface.fdim(iface) == 1) {
      int nvals1i, nvals1j, nvals2;
      int il1 = 0, ih1 = 0, jl1 = 0, jh1 = 0;
      if (r.box(jgrid).smallEnd(0) == b.smallEnd(0))
        jl1 = w;
      if (r.box(jgrid).bigEnd(0) == b.bigEnd(0))
        jh1 = w;
      if (r.box(igrid).smallEnd(0) == b.smallEnd(0))
        il1 = w;
      if (r.box(igrid).bigEnd(0) == b.bigEnd(0))
        ih1 = w;
      nvals1i = b.length(0) + il1 + ih1;
      nvals1j = b.length(0) + jl1 + jh1;
      nvals2 = b.length(2);
      stridi1 = r[igrid].box().length(0);
      stridi2 = stridi1 * r[igrid].box().length(1);
      stridj1 = r[jgrid].box().length(0);
      stridj2 = stridj1 * r[jgrid].box().length(1);
      dstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - il1 - r[igrid].box().smallEnd(0) +
	        stridi1 * (b.smallEnd(1) + 1 - r[igrid].box().smallEnd(1)) +
	        stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
      sstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - jl1 - r[igrid].box().smallEnd(0) +
	        stridi1 * (b.smallEnd(1) - 1 - r[igrid].box().smallEnd(1)) +
	        stridi2 * (b.smallEnd(2) - r[igrid].box().smallEnd(2));
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - jl1 - r[jgrid].box().smallEnd(0) +
	        stridj1 * (b.smallEnd(1) - 1 - r[jgrid].box().smallEnd(1)) +
	        stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
      sstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - il1 - r[jgrid].box().smallEnd(0) +
	        stridj1 * (b.smallEnd(1) + 1 - r[jgrid].box().smallEnd(1)) +
	        stridj2 * (b.smallEnd(2) - r[jgrid].box().smallEnd(2));
      stridi1 = 1;
      stridj1 = 1;
      set(iset++, dstarti, sstartj, stridi1, stridi2,
	  stridj1, stridj2, nvals1i, nvals2);
      set(iset++, dstartj, sstarti, stridj1, stridj2,
	  stridi1, stridi2, nvals1j, nvals2);
    }
    else {
      int nvals1i, nvals1j, nvals2i, nvals2j;
      int il1 = 0, ih1 = 0, jl1 = 0, jh1 = 0;
      int il2 = 0, ih2 = 0, jl2 = 0, jh2 = 0;
      if (r.box(jgrid).smallEnd(0) == b.smallEnd(0))
        jl1 = w;
      if (r.box(jgrid).bigEnd(0) == b.bigEnd(0))
        jh1 = w;
      if (r.box(igrid).smallEnd(0) == b.smallEnd(0))
        il1 = w;
      if (r.box(igrid).bigEnd(0) == b.bigEnd(0))
        ih1 = w;
      if (r.box(jgrid).smallEnd(1) == b.smallEnd(1))
        jl2 = w;
      if (r.box(jgrid).bigEnd(1) == b.bigEnd(1))
        jh2 = w;
      if (r.box(igrid).smallEnd(1) == b.smallEnd(1))
        il2 = w;
      if (r.box(igrid).bigEnd(1) == b.bigEnd(1))
        ih2 = w;
      nvals1i = b.length(0) + il1 + ih1;
      nvals1j = b.length(0) + jl1 + jh1;
      nvals2i = b.length(1) + il2 + ih2;
      nvals2j = b.length(1) + jl2 + jh2;
      stridi1 = r[igrid].box().length(0);
      stridi2 = stridi1 * r[igrid].box().length(1);
      stridj1 = r[jgrid].box().length(0);
      stridj2 = stridj1 * r[jgrid].box().length(1);
      dstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - il1 - r[igrid].box().smallEnd(0) +
	        stridi1 * (b.smallEnd(1) - il2 - r[igrid].box().smallEnd(1)) +
	        stridi2 * (b.smallEnd(2) + 1 - r[igrid].box().smallEnd(2));
      sstarti = r[igrid].dataPtr() - baseptr +
	        b.smallEnd(0) - jl1 - r[igrid].box().smallEnd(0) +
	        stridi1 * (b.smallEnd(1) - jl2 - r[igrid].box().smallEnd(1)) +
	        stridi2 * (b.smallEnd(2) - 1 - r[igrid].box().smallEnd(2));
      dstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - jl1 - r[jgrid].box().smallEnd(0) +
	        stridj1 * (b.smallEnd(1) - jl2 - r[jgrid].box().smallEnd(1)) +
	        stridj2 * (b.smallEnd(2) - 1 - r[jgrid].box().smallEnd(2));
      sstartj = r[jgrid].dataPtr() - baseptr +
	        b.smallEnd(0) - il1 - r[jgrid].box().smallEnd(0) +
	        stridj1 * (b.smallEnd(1) - il2 - r[jgrid].box().smallEnd(1)) +
	        stridj2 * (b.smallEnd(2) + 1 - r[jgrid].box().smallEnd(2));
      stridi2 = stridi1;
      stridj2 = stridj1;
      stridi1 = 1;
      stridj1 = 1;
      set(iset++, dstarti, sstartj, stridi1, stridi2,
                 stridj1, stridj2, nvals1i, nvals2i);
      set(iset++, dstartj, sstarti, stridj1, stridj2,
                 stridi1, stridi2, nvals1j, nvals2j);
    }
#endif
  }

  bdy.set_border_cache(this, nsets, iset, r, interface, w);
  nsets = iset;
}

unroll_cache::unroll_cache(MultiFab& r)
{
  nsets = r.length();
  assert(nsets > 0);

  Real *baseptr = r[0].dataPtr();
  ptr = baseptr;

#if (BL_SPACEDIM == 2)
  start = new int[3 * nsets];
  strid = start + nsets;
  nvals = start + 2 * nsets;
#else
  start  = new int[4 * nsets];
  strid1 = start + nsets;
  strid2 = start + 2 * nsets;
  nvals  = start + 3 * nsets;
#endif

  for (int igrid = 0; igrid < r.length(); igrid++) {
    set(igrid, r[igrid].dataPtr() - baseptr,
	r[igrid].box().length(0),
#if (BL_SPACEDIM == 3)
	r[igrid].box().length(0) * r[igrid].box().length(1),
#endif
	r[igrid].box().numPts());
  }
}
