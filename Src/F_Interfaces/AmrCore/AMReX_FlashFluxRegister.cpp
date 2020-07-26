#include <AMReX_FlashFluxRegister.H>
#include <AMReX_OpenMP.H>

namespace amrex {

FlashFluxRegister::FlashFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                      const DistributionMapping& fdm, const DistributionMapping& cdm,
                                      const Geometry& fgeom, const Geometry& cgeom,
                                      IntVect const& ref_ratio, int nvar)
{
    define(fba,cba,fdm,cdm,fgeom,cgeom,ref_ratio,nvar);
}

void FlashFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
                                const DistributionMapping& fdm, const DistributionMapping& cdm,
                                const Geometry& fgeom, const Geometry& cgeom,
                                IntVect const& a_ref_ratio, int nvar)
{
    constexpr int ref_ratio = 2;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_ref_ratio == IntVect(ref_ratio),
                                     "FlashFluxRegister: refinement ratio != 2");

    m_fine_grids = fba;
    m_crse_grids = cba;
    m_fine_dmap = fdm;
    m_crse_dmap = cdm;
    m_fine_geom = fgeom;
    m_crse_geom = cgeom;
    m_ncomp = nvar;

    const int myproc = ParallelDescriptor::MyProc();

    // For a fine Box, there is at most one face per direction abutting coarse level.
    {
        BoxArray const& fndba = amrex::convert(fba,IntVect::TheNodeVector());
        Array<BoxList,AMREX_SPACEDIM> bl
        {{AMREX_D_DECL(BoxList(IndexType(IntVect::TheDimensionVector(0))),
                       BoxList(IndexType(IntVect::TheDimensionVector(1))),
                       BoxList(IndexType(IntVect::TheDimensionVector(2))))}};
        Array<Vector<int>,AMREX_SPACEDIM> procmap;
        Array<Vector<int>,AMREX_SPACEDIM> my_global_indices;
        std::vector<std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = m_fine_geom.periodicity().shiftIntVect();
        Vector<std::pair<Orientation,Box> > faces;
        for (int i = 0, nfines = fba.size(); i < nfines; ++i)
        {
            Box const& ccbx = fba[i];
            Box const& ndbx = amrex::surroundingNodes(ccbx);
            faces.clear();
            for (OrientationIter orit; orit.isValid(); ++orit) {
                const Orientation ori = orit();
                faces.emplace_back(std::make_pair(ori,amrex::bdryNode(ndbx,ori)));
            }
            for (auto const& shift : pshifts) {
                if (!faces.empty()) {
                    fndba.intersections(ndbx+shift, isects);
                    for (auto const& isect : isects) {
                        Box const& b = isect.second-shift;
                        faces.erase(std::remove_if(faces.begin(), faces.end(),
                                                   [&] (std::pair<Orientation,Box> const& x)
                                                   { return x.second == b; }),
                                    faces.end());
                    }
                }
            }
            for (auto const& face : faces) {  // coarse/fine boundary faces
                const int dir = face.first.coordDir();
                bl[dir].push_back(amrex::coarsen(amrex::bdryNode(ccbx,face.first),
                                                 ref_ratio));
                procmap[dir].push_back(fdm[i]);
                if (fdm[i] == myproc) my_global_indices[dir].push_back(i);
            }
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (bl[idim].isNotEmpty()) {
                BoxArray faces_ba(std::move(bl[idim]));
                DistributionMapping faces_dm(std::move(procmap[idim]));
                m_fine_fluxes[idim].define(faces_ba, faces_dm, m_ncomp, 0);
                for (MFIter mfi(m_fine_fluxes[idim]); mfi.isValid(); ++mfi) {
                    const int gi = my_global_indices[idim][mfi.LocalIndex()];
                    auto found = m_fine_map.find(gi);
                    if (found != m_fine_map.end()) {
                        found->second[idim] = &(m_fine_fluxes[idim][mfi]);
                    } else {
                        Array<FArrayBox*,AMREX_SPACEDIM> t{{AMREX_D_DECL(nullptr,nullptr,nullptr)}};
                        t[idim] = &(m_fine_fluxes[idim][mfi]);
                        m_fine_map.insert(std::make_pair(gi,t));
                    }
                }
            }
        }
    }

    // For a coarse Box, there are at most two faces per direction abutting fine level.
    {
        BoxArray const fba_coarsened = amrex::coarsen(fba,ref_ratio);
        Array<BoxList,AMREX_SPACEDIM> bl
        {{AMREX_D_DECL(BoxList(IndexType(IntVect::TheDimensionVector(0))),
                       BoxList(IndexType(IntVect::TheDimensionVector(1))),
                       BoxList(IndexType(IntVect::TheDimensionVector(2))))}};
        Array<Vector<int>,AMREX_SPACEDIM> procmap;
        Array<Vector<int>,AMREX_SPACEDIM> my_global_indices;
        std::vector<std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = m_crse_geom.periodicity().shiftIntVect();
        Vector<std::pair<Orientation,Box> > cell_faces;
        Vector<Orientation> crsefine_faces;
        for (int i = 0, ncrses = cba.size(); i < ncrses; ++i)
        {
            Box const& ccbx = cba[i];
            Box const& ccbxg1 = amrex::grow(ccbx,1);
            cell_faces.clear();
            crsefine_faces.clear();
            for (OrientationIter orit; orit.isValid(); ++orit) {
                const Orientation ori = orit();
                cell_faces.emplace_back(std::make_pair(ori,amrex::adjCell(ccbx,ori)));
            }
            for (auto const& shift : pshifts) {
                if (!cell_faces.empty()) {
                    fba_coarsened.intersections(ccbxg1+shift, isects);
                    for (auto const& isect : isects) {
                        Box const& b = isect.second-shift;
                        if (ccbx.intersects(b)) {
                            // This coarse box is covered by fine
                            cell_faces.clear();
                            crsefine_faces.clear();
                        } else {
                            auto found = std::find_if(cell_faces.begin(), cell_faces.end(),
                                                      [&] (std::pair<Orientation,Box> const& x)
                                                      { return x.second.contains(b); });
                            if (found != cell_faces.end()) {
                                crsefine_faces.push_back(found->first);
                                cell_faces.erase(found);
                            }
                        }
                    }
                }
            }
            for (auto const& face : crsefine_faces) {
                const int dir = face.coordDir();
                bl[dir].push_back(amrex::bdryNode(ccbx,face));
                procmap[dir].push_back(cdm[i]);
                if (cdm[i] == myproc) my_global_indices[dir].push_back(i);
            }
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (bl[idim].isNotEmpty()) {
                BoxArray faces_ba(std::move(bl[idim]));
                DistributionMapping faces_dm(std::move(procmap[idim]));
                m_crse_fluxes[idim].define(faces_ba, faces_dm, m_ncomp, 0);
                for (MFIter mfi(m_crse_fluxes[idim]); mfi.isValid(); ++mfi) {
                    const int gi = my_global_indices[idim][mfi.LocalIndex()];
                    const Orientation::Side side = (mfi.validbox() == amrex::bdryLo(cba[gi],idim))
                        ? Orientation::low : Orientation::high;
                    const int index = Orientation(idim,side);
                    auto found = m_crse_map.find(gi);
                    if (found != m_crse_map.end()) {
                        found->second[index] = &(m_crse_fluxes[idim][mfi]);
                    } else {
                        Array<FArrayBox*,2*AMREX_SPACEDIM> t{{AMREX_D_DECL(nullptr,nullptr,nullptr),
                                    AMREX_D_DECL(nullptr,nullptr,nullptr)}};
                        t[index] = &(m_crse_fluxes[idim][mfi]);
                        m_crse_map.insert(std::make_pair(gi,t));
                    }
                }
            }
        }
    }

    int nthreads = OpenMP::get_max_threads();
    m_h_ifd.resize(nthreads);
    m_d_ifd.resize(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        m_h_ifd[i].resize(m_ncomp, 0);
    }
}

void FlashFluxRegister::store (int fine_global_index, int dir, FArrayBox const& fine_flux, Real sf)
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            if (dir == 0) {
#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,0,0,n) = src(2*i,0,0,n)*sf;
                });
#endif
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,0,n) = (src(2*i,2*j  ,0,n) +
                                     src(2*i,2*j+1,0,n)) * (0.5*sf);
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = (src(2*i,2*j  ,2*k  ,n) +
                                     src(2*i,2*j+1,2*k  ,n) +
                                     src(2*i,2*j  ,2*k+1,n) +
                                     src(2*i,2*j+1,2*k+1,n)) * (0.25*sf);
                });
#endif
            } else if (dir == 1) {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,0,n) = (src(2*i  ,2*j,0,n) +
                                     src(2*i+1,2*j,0,n)) * (0.5*sf);
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = (src(2*i  ,2*j,2*k  ,n) +
                                     src(2*i+1,2*j,2*k  ,n) +
                                     src(2*i  ,2*j,2*k+1,n) +
                                     src(2*i+1,2*j,2*k+1,n)) * (0.25*sf);
                });
#endif
            } else {
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = (src(2*i  ,2*j  ,2*k,n) +
                                     src(2*i+1,2*j  ,2*k,n) +
                                     src(2*i  ,2*j+1,2*k,n) +
                                     src(2*i+1,2*j+1,2*k,n)) * (0.25*sf);
                });
#endif
            }
        }
    }
}

void FlashFluxRegister::store (int fine_global_index, int dir, FArrayBox const& fine_flux,
                               FArrayBox const& fine_area, Real sf)
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            Array4<Real const> const& area = fine_area.const_array();
            if (dir == 0) {
#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,0,0,n) = src(2*i,0,0,n)*area(2*i,0,0)*sf;
                });
#endif
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,0,n) = (src(2*i,2*j  ,0,n)*area(2*i,2*j  ,0) +
                                     src(2*i,2*j+1,0,n)*area(2*i,2*j+1,0)) * sf;
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = (src(2*i,2*j  ,2*k  ,n)*area(2*i,2*j  ,2*k  ) +
                                     src(2*i,2*j+1,2*k  ,n)*area(2*i,2*j+1,2*k  ) +
                                     src(2*i,2*j  ,2*k+1,n)*area(2*i,2*j  ,2*k+1) +
                                     src(2*i,2*j+1,2*k+1,n)*area(2*i,2*j+1,2*k+1)) * sf;
                });
#endif
            } else if (dir == 1) {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,0,n) = (src(2*i  ,2*j,0,n)*area(2*i  ,2*j,0) +
                                     src(2*i+1,2*j,0,n)*area(2*i+1,2*j,0)) * sf;
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = (src(2*i  ,2*j,2*k  ,n)*area(2*i  ,2*j,2*k  ) +
                                     src(2*i+1,2*j,2*k  ,n)*area(2*i+1,2*j,2*k  ) +
                                     src(2*i  ,2*j,2*k+1,n)*area(2*i  ,2*j,2*k+1) +
                                     src(2*i+1,2*j,2*k+1,n)*area(2*i+1,2*j,2*k+1)) * sf;
                });
#endif
            } else {
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = (src(2*i  ,2*j  ,2*k,n)*area(2*i  ,2*j  ,2*k) +
                                     src(2*i+1,2*j  ,2*k,n)*area(2*i+1,2*j  ,2*k) +
                                     src(2*i  ,2*j+1,2*k,n)*area(2*i  ,2*j+1,2*k) +
                                     src(2*i+1,2*j+1,2*k,n)*area(2*i+1,2*j+1,2*k)) * sf;
                });
#endif
            }
        }
    }
}

void FlashFluxRegister::store (int fine_global_index, int dir, FArrayBox const& fine_flux,
                               FArrayBox const& fine_area, const int* isFluxDensity, Real sf)
{
    auto& h_ifd = m_h_ifd[OpenMP::get_thread_num()];
    auto& d_ifd = m_d_ifd[OpenMP::get_thread_num()];

    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            bool allsame = true;
            for (int n = 0; n < m_ncomp; ++n) {
                if (h_ifd[n] != isFluxDensity[n]) {
                    allsame = false;
                    h_ifd[n] = isFluxDensity[n];
                }
            }
            if (d_ifd.empty()) {
                allsame = false;
                d_ifd.resize(m_ncomp);
            }
            if (not allsame) {
                Gpu::copyAsync(Gpu::HostToDevice(), h_ifd.begin(), h_ifd.end(), d_ifd.begin());
            }

            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            Array4<Real const> const& area = fine_area.const_array();
            const int* ifd = d_ifd.data();
            if (dir == 0) {
#if (AMREX_SPACEDIM == 1)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,0,0,n) = src(2*i,0,0,n)*area(2*i,0,0)*sf;
                    } else {
                        dest(i,0,0,n) = src(2*i,0,0,n)*sf;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,j,0,n) = (src(2*i,2*j  ,0,n)*area(2*i,2*j  ,0) +
                                         src(2*i,2*j+1,0,n)*area(2*i,2*j+1,0)) * sf;
                    } else {
                        dest(i,j,0,n) = (src(2*i,2*j  ,0,n) +
                                         src(2*i,2*j+1,0,n)) * sf;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,j,k,n) = (src(2*i,2*j  ,2*k  ,n)*area(2*i,2*j  ,2*k  ) +
                                         src(2*i,2*j+1,2*k  ,n)*area(2*i,2*j+1,2*k  ) +
                                         src(2*i,2*j  ,2*k+1,n)*area(2*i,2*j  ,2*k+1) +
                                         src(2*i,2*j+1,2*k+1,n)*area(2*i,2*j+1,2*k+1)) * sf;
                    } else {
                        dest(i,j,k,n) = (src(2*i,2*j  ,2*k  ,n) +
                                         src(2*i,2*j+1,2*k  ,n) +
                                         src(2*i,2*j  ,2*k+1,n) +
                                         src(2*i,2*j+1,2*k+1,n)) * sf;
                    }
                });
#endif
            } else if (dir == 1) {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,j,0,n) = (src(2*i  ,2*j,0,n)*area(2*i  ,2*j,0) +
                                         src(2*i+1,2*j,0,n)*area(2*i+1,2*j,0)) * sf;
                    } else {
                        dest(i,j,0,n) = (src(2*i  ,2*j,0,n) +
                                         src(2*i+1,2*j,0,n)) * sf;
                    }
                });
#endif
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,j,k,n) = (src(2*i  ,2*j,2*k  ,n)*area(2*i  ,2*j,2*k  ) +
                                         src(2*i+1,2*j,2*k  ,n)*area(2*i+1,2*j,2*k  ) +
                                         src(2*i  ,2*j,2*k+1,n)*area(2*i  ,2*j,2*k+1) +
                                         src(2*i+1,2*j,2*k+1,n)*area(2*i+1,2*j,2*k+1)) * sf;
                    } else {
                        dest(i,j,k,n) = (src(2*i  ,2*j,2*k  ,n) +
                                         src(2*i+1,2*j,2*k  ,n) +
                                         src(2*i  ,2*j,2*k+1,n) +
                                         src(2*i+1,2*j,2*k+1,n)) * sf;
                    }
                });
#endif
            } else {
#if (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,j,k,n) = (src(2*i  ,2*j  ,2*k,n)*area(2*i  ,2*j  ,2*k) +
                                         src(2*i+1,2*j  ,2*k,n)*area(2*i+1,2*j  ,2*k) +
                                         src(2*i  ,2*j+1,2*k,n)*area(2*i  ,2*j+1,2*k) +
                                         src(2*i+1,2*j+1,2*k,n)*area(2*i+1,2*j+1,2*k)) * sf;
                    } else {
                        dest(i,j,k,n) = (src(2*i  ,2*j  ,2*k,n) +
                                         src(2*i+1,2*j  ,2*k,n) +
                                         src(2*i  ,2*j+1,2*k,n) +
                                         src(2*i+1,2*j+1,2*k,n)) * sf;
                    }
                });
#endif
            }
        }
    }
}

void FlashFluxRegister::communicate ()
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_crse_fluxes[idim].ParallelCopy(m_fine_fluxes[idim], m_crse_geom.periodicity());
    }
}

void FlashFluxRegister::load (int crse_global_index, int dir, FArrayBox& crse_flux, Real sf) const
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_crse_map.find(crse_global_index);
    if (found != m_crse_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,2*AMREX_SPACEDIM> const& fab_a = found->second;
        Array4<Real> const& dest = crse_flux.array();
        for (int index = dir; index < 2*AMREX_SPACEDIM; index += AMREX_SPACEDIM) {
            if (fab_a[index]) {
                Box const& b = fab_a[index]->box();
                Array4<Real const> const& src = fab_a[index]->const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = src(i,j,k,n) * sf;
                });
            }
        }
    }
}

void FlashFluxRegister::load (int crse_global_index, int dir, FArrayBox& crse_flux,
                              FArrayBox const& cflux, Real sf_f, Real sf_c) const
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_crse_map.find(crse_global_index);
    if (found != m_crse_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,2*AMREX_SPACEDIM> const& fab_a = found->second;
        Array4<Real> const& dest = crse_flux.array();
        Array4<Real const> const& csrc = cflux.const_array();
        for (int index = dir; index < 2*AMREX_SPACEDIM; index += AMREX_SPACEDIM) {
            if (fab_a[index]) {
                Box const& b = fab_a[index]->box();
                Array4<Real const> const& src = fab_a[index]->const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = src(i,j,k,n) * sf_f + csrc(i,j,k,n) * sf_c;
                });
            }
        }
    }
}

void FlashFluxRegister::load (int crse_global_index, int dir, FArrayBox& crse_flux,
                              FArrayBox const& area_fab) const
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_crse_map.find(crse_global_index);
    if (found != m_crse_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,2*AMREX_SPACEDIM> const& fab_a = found->second;
        Array4<Real> const& dest = crse_flux.array();
        Array4<Real const> const& area = area_fab.const_array();
        for (int index = dir; index < 2*AMREX_SPACEDIM; index += AMREX_SPACEDIM) {
            if (fab_a[index]) {
                Box const& b = fab_a[index]->box();
                Array4<Real const> const& src = fab_a[index]->const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = src(i,j,k,n)/area(i,j,k);
                });
            }
        }
    }
}

void FlashFluxRegister::load (int crse_global_index, int dir, FArrayBox& crse_flux,
                              FArrayBox const& cflux, FArrayBox const& area_fab,
                              Real sf_f, Real sf_c) const
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_crse_map.find(crse_global_index);
    if (found != m_crse_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,2*AMREX_SPACEDIM> const& fab_a = found->second;
        Array4<Real> const& dest = crse_flux.array();
        Array4<Real const> const& csrc = cflux.const_array();
        Array4<Real const> const& area = area_fab.const_array();
        for (int index = dir; index < 2*AMREX_SPACEDIM; index += AMREX_SPACEDIM) {
            if (fab_a[index]) {
                Box const& b = fab_a[index]->box();
                Array4<Real const> const& src = fab_a[index]->const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = src(i,j,k,n)/area(i,j,k) * sf_f + csrc(i,j,k,n) * sf_c;
                });
            }
        }
    }
}

void FlashFluxRegister::load (int crse_global_index, int dir, FArrayBox& crse_flux,
                              FArrayBox const& cflux, FArrayBox const& area_fab,
                              const int* isFluxDensity, Real sf_f, Real sf_c) const
{
    auto& h_ifd = m_h_ifd[OpenMP::get_thread_num()];
    auto& d_ifd = m_d_ifd[OpenMP::get_thread_num()];

    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_crse_map.find(crse_global_index);
    if (found != m_crse_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,2*AMREX_SPACEDIM> const& fab_a = found->second;
        Array4<Real> const& dest = crse_flux.array();
        Array4<Real const> const& csrc = cflux.const_array();
        Array4<Real const> const& area = area_fab.const_array();
        int d_ifd_ready = false;
        for (int index = dir; index < 2*AMREX_SPACEDIM; index += AMREX_SPACEDIM) {
            if (fab_a[index]) {
                if (not d_ifd_ready) {
                    d_ifd_ready = true;
                    bool allsame = true;
                    for (int n = 0; n < m_ncomp; ++n) {
                        if (h_ifd[n] != isFluxDensity[n]) {
                            allsame = false;
                            h_ifd[n] = isFluxDensity[n];
                        }
                    }
                    if (d_ifd.empty()) {
                        allsame = false;
                        d_ifd.resize(m_ncomp);
                    }
                    if (not allsame) {
                        Gpu::copyAsync(Gpu::HostToDevice(), h_ifd.begin(), h_ifd.end(),
                                       d_ifd.begin());
                    }
                }
                Box const& b = fab_a[index]->box();
                Array4<Real const> const& src = fab_a[index]->const_array();
                const int* ifd = d_ifd.data();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    if (ifd[n]) {
                        dest(i,j,k,n) = src(i,j,k,n)/area(i,j,k) * sf_f + csrc(i,j,k,n) * sf_c;
                    } else {
                        dest(i,j,k,n) = src(i,j,k,n)             * sf_f + csrc(i,j,k,n) * sf_c;
                    }
                });
            }
        }
    }
}

}
