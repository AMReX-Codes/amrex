#ifndef _AMReX_EBDataVarMacros_H_
#define _AMReX_EBDataVarMacros_H_

namespace amrex {

#if (AMREX_SPACEDIM == 2)

static constexpr int V_VOLFRAC       = 0;
static constexpr int V_BNDAREA       = 1;
static constexpr int V_VOLCENTROIDX  = 2;
static constexpr int V_VOLCENTROIDY  = 3;
static constexpr int V_BNDCENTROIDX  = 4;
static constexpr int V_BNDCENTROIDY  = 5;
static constexpr int V_NORMALX       = 6;
static constexpr int V_NORMALY       = 7;
static constexpr int V_VOLNUMBER     = 8;

static constexpr int F_AREAFRAC      = 0;
static constexpr int F_FACECENTROIDX = 1;
static constexpr int F_FACECENTROIDY = 2;
static constexpr int F_FACENUMBER    = 3;

#elif (AMREX_SPACEDIM == 3)

static constexpr int V_VOLFRAC      = 0;
static constexpr int V_BNDAREA      = 1;
static constexpr int V_VOLCENTROIDX = 2;
static constexpr int V_VOLCENTROIDY = 3;
static constexpr int V_VOLCENTROIDZ = 4;
static constexpr int V_BNDCENTROIDX = 5;
static constexpr int V_BNDCENTROIDY = 6;
static constexpr int V_BNDCENTROIDZ = 7;
static constexpr int V_NORMALX      = 8;
static constexpr int V_NORMALY      = 9;
static constexpr int V_NORMALZ      = 10;
static constexpr int V_VOLNUMBER    = 11;

static constexpr int F_AREAFRAC      = 0;
static constexpr int F_FACECENTROIDX = 1;
static constexpr int F_FACECENTROIDY = 2;
static constexpr int F_FACECENTROIDZ = 3;
static constexpr int F_FACENUMBER    = 4;

#endif

}

#endif
