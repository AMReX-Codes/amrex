#ifndef AMREX_PARTICLE_LOCATOR_H_
#define AMREX_PARTICLE_LOCATOR_H_

#include <AMReX_ParGDB.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Tuple.H>

namespace amrex
{

template <class BinIteratorFactory>
struct AssignGrid
{
    BinIteratorFactory m_bif;

    Dim3 m_lo;
    Dim3 m_hi;
    Dim3 m_bin_size;
    Dim3 m_num_bins;

    Box m_domain;
    GpuArray<Real, AMREX_SPACEDIM> m_plo;
    GpuArray<Real, AMREX_SPACEDIM> m_dxi;

    AMREX_GPU_HOST_DEVICE
    AssignGrid () {}

    AssignGrid (BinIteratorFactory a_bif,
                const IntVect& a_bins_lo, const IntVect& a_bins_hi, const IntVect& a_bin_size,
                const IntVect& a_num_bins, const Geometry& a_geom)
        : m_bif(a_bif),
          m_lo(a_bins_lo.dim3()), m_hi(a_bins_hi.dim3()), m_bin_size(a_bin_size.dim3()),
          m_num_bins(a_num_bins.dim3()), m_domain(a_geom.Domain()),
          m_plo(a_geom.ProbLoArray()), m_dxi(a_geom.InvCellSizeArray())
        {
            // clamp bin size and num_bins to 1 for AMREX_SPACEDIM < 3
            m_bin_size.x = amrex::max(m_bin_size.x, 1);
            m_bin_size.y = amrex::max(m_bin_size.y, 1);
            m_bin_size.z = amrex::max(m_bin_size.z, 1);

            m_num_bins.x = amrex::max(m_num_bins.x, 1);
            m_num_bins.y = amrex::max(m_num_bins.y, 1);
            m_num_bins.z = amrex::max(m_num_bins.z, 1);
        }

    template <typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const P& p, int nGrow=0) const noexcept
    {
        const auto iv = getParticleCell(p, m_plo, m_dxi, m_domain);
        return this->operator()(iv, nGrow);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const IntVect& iv, int nGrow=0) const noexcept
    {
        const auto lo = iv.dim3();
        int ix_lo = amrex::max((lo.x - nGrow - m_lo.x) / m_bin_size.x - 1, 0);
        int iy_lo = amrex::max((lo.y - nGrow - m_lo.y) / m_bin_size.y - 1, 0);
        int iz_lo = amrex::max((lo.z - nGrow - m_lo.z) / m_bin_size.z - 1, 0);

        int ix_hi = amrex::min((lo.x + nGrow - m_lo.x) / m_bin_size.x, m_num_bins.x-1);
        int iy_hi = amrex::min((lo.y + nGrow - m_lo.y) / m_bin_size.y, m_num_bins.y-1);
        int iz_hi = amrex::min((lo.z + nGrow - m_lo.z) / m_bin_size.z, m_num_bins.z-1);

        for (int ii = ix_lo; ii <= ix_hi; ++ii) {
            for (int jj = iy_lo; jj <= iy_hi; ++jj) {
                for (int kk = iz_lo; kk <= iz_hi; ++kk) {
                    int index = (ii * m_num_bins.y + jj) * m_num_bins.z + kk;
                    for (const auto& nbor : m_bif.getBinIterator(index)) {
                        Box bx = nbor.second;
                        bx.grow(nGrow);
                        if (bx.contains(iv)) return nbor.first;
                    }
                }
            }
        }

        return -1;
    }
};

template <class Bins>
class ParticleLocator
{
public:

    using BinIteratorFactory = typename Bins::BinIteratorFactory;

    ParticleLocator () : m_defined(false) {}

    void build (const BoxArray& ba, const Geometry& geom)
    {
        m_defined = true;
        m_ba = ba;
        m_geom = geom;
        int num_boxes = ba.size();
        m_host_boxes.resize(0);
        for (int i = 0; i < num_boxes; ++i) m_host_boxes.push_back(ba[i]);

        m_device_boxes.resize(num_boxes);
        Gpu::copy(Gpu::hostToDevice, m_host_boxes.begin(), m_host_boxes.end(), m_device_boxes.begin());

        // compute the lo, hi and the max box size in each direction
        ReduceOps<AMREX_D_DECL(ReduceOpMin, ReduceOpMin, ReduceOpMin),
                  AMREX_D_DECL(ReduceOpMax, ReduceOpMax, ReduceOpMax),
                  AMREX_D_DECL(ReduceOpMax, ReduceOpMax, ReduceOpMax)> reduce_op;
        ReduceData<AMREX_D_DECL(int, int, int),
                   AMREX_D_DECL(int, int, int),
                   AMREX_D_DECL(int, int, int)> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const auto boxes_ptr = m_device_boxes.dataPtr();
        reduce_op.eval(num_boxes, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            const Box& box = boxes_ptr[i];
            IntVect lo = box.smallEnd();
            IntVect hi = box.bigEnd();
            IntVect si = box.length();
            return {AMREX_D_DECL(lo[0], lo[1], lo[2]),
                    AMREX_D_DECL(hi[0], hi[1], hi[2]),
                    AMREX_D_DECL(si[0], si[1], si[2])};
        });

        ReduceTuple hv = reduce_data.value();

        m_bins_lo  = IntVect(AMREX_D_DECL(amrex::get<0>(hv),
                                          amrex::get<1>(hv),
                                          amrex::get<2>(hv)));
        m_bins_hi  = IntVect(AMREX_D_DECL(amrex::get< AMREX_SPACEDIM  >(hv),
                                          amrex::get< AMREX_SPACEDIM+1>(hv),
                                          amrex::get< AMREX_SPACEDIM+2>(hv)));
        m_bin_size = IntVect(AMREX_D_DECL(amrex::get<2*AMREX_SPACEDIM>(hv),
                                          amrex::get<2*AMREX_SPACEDIM+1>(hv),
                                          amrex::get<2*AMREX_SPACEDIM+2>(hv)));

        m_num_bins = (m_bins_hi - m_bins_lo + m_bin_size) / m_bin_size;

        Box bins_box(IntVect::TheZeroVector(), m_num_bins-IntVect::TheUnitVector());
        IntVect bin_size = m_bin_size;
        IntVect bins_lo = m_bins_lo;
        m_bins.build(num_boxes, boxes_ptr, bins_box,
                     [=] AMREX_GPU_HOST_DEVICE (const Box& box) noexcept -> IntVect
                     {
                         return (box.smallEnd() - bins_lo) / bin_size;
                     });
    }

    void setGeometry (const Geometry& a_geom) noexcept
    {
        AMREX_ASSERT(m_defined);
        m_geom = a_geom;
    }

    AssignGrid<BinIteratorFactory> getGridAssignor () const noexcept
    {
        AMREX_ASSERT(m_defined);
        return AssignGrid<BinIteratorFactory>(m_bins.getBinIteratorFactory(),
                                              m_bins_lo, m_bins_hi, m_bin_size, m_num_bins, m_geom);
    }

    bool isValid (const BoxArray& ba) const noexcept
    {
        if (m_defined) return BoxArray::SameRefs(m_ba, ba);
        return false;
    }

protected:

    bool m_defined;

    BoxArray m_ba;
    Geometry m_geom;

    IntVect m_bins_lo;
    IntVect m_bins_hi;
    IntVect m_bin_size;
    IntVect m_num_bins;

    Bins m_bins;

    Gpu::HostVector<Box> m_host_boxes;
    Gpu::DeviceVector<Box> m_device_boxes;
};

template <class BinIteratorFactory>
struct AmrAssignGrid
{
    const AssignGrid<BinIteratorFactory>* m_funcs;
    std::size_t m_size;

    AmrAssignGrid(const AssignGrid<BinIteratorFactory>* a_funcs, std::size_t a_size)
        : m_funcs(a_funcs), m_size(a_size)
        {}

    template <typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    GpuTuple<int, int> operator() (const P& p, int lev_min=-1, int lev_max=-1, int nGrow=0) const noexcept
    {
        lev_min = (lev_min == -1) ? 0 : lev_min;
        lev_max = (lev_max == -1) ? m_size - 1 : lev_max;

        for (int lev = lev_max; lev >= lev_min; --lev)
        {
            int grid = m_funcs[lev](p);
            if (grid >= 0) return makeTuple(grid, lev);
        }

        for (int lev = lev_min; lev >= lev_min; --lev)
        {
            int grid = m_funcs[lev](p, nGrow);
            if (grid >= 0) return makeTuple(grid, lev);
        }

        return makeTuple(-1, -1);
    }
};

template <class Bins>
class AmrParticleLocator
{
public:
    using BinIteratorFactory = typename Bins::BinIteratorFactory;

private:
    Vector<ParticleLocator<Bins> > m_locators;
    Gpu::DeviceVector<AssignGrid<BinIteratorFactory> > m_grid_assignors;
    bool m_defined = false;

public:

    AmrParticleLocator() {}

    AmrParticleLocator(const Vector<BoxArray>& a_ba,
                       const Vector<Geometry>& a_geom)
    {
        build(a_ba, a_geom);
    }

    AmrParticleLocator(const ParGDBBase* a_gdb)
    {
        build(a_gdb);
    }

    void build (const Vector<BoxArray>& a_ba,
                const Vector<Geometry>& a_geom)
    {
        m_defined = true;
        int num_levels = a_ba.size();
        m_locators.resize(num_levels);
        m_grid_assignors.resize(num_levels);
        for (int lev = 0; lev < num_levels; ++lev)
        {
            m_locators[lev].build(a_ba[lev], a_geom[lev]);
            m_grid_assignors[lev] = m_locators[lev].getGridAssignor();
        }
    }

    void build (const ParGDBBase* a_gdb)
    {
        Vector<BoxArray> ba;
        Vector<Geometry> geom;
        int num_levels = a_gdb->finestLevel()+1;
        for (int lev = 0; lev < num_levels; ++lev)
        {
            ba.push_back(a_gdb->ParticleBoxArray(lev));
            geom.push_back(a_gdb->Geom(lev));
        }
        build(ba, geom);
    }

    bool isValid (const Vector<BoxArray>& a_ba) const
    {
        if ( !m_defined or (m_locators.size() == 0) or
             (m_locators.size() != a_ba.size()) ) return false;
        bool all_valid = true;
        int num_levels = m_locators.size();
        for (int lev = 0; lev < num_levels; ++lev)
            all_valid = all_valid && m_locators[lev].isValid(a_ba[lev]);
        return all_valid;
    }

    bool isValid (const ParGDBBase* a_gdb) const
    {
        Vector<BoxArray> ba;
        int num_levels = a_gdb->finestLevel()+1;
        for (int lev = 0; lev < num_levels; ++lev)
            ba.push_back(a_gdb->ParticleBoxArray(lev));
        return this->isValid(ba);
    }

    void setGeometry (const ParGDBBase* a_gdb)
    {
        int num_levels = a_gdb->finestLevel()+1;
        for (int lev = 0; lev < num_levels; ++lev)
        {
            m_locators[lev].setGeometry(a_gdb->Geom(lev));
            m_grid_assignors[lev] = m_locators[lev].getGridAssignor();
        }
    }

    AmrAssignGrid<BinIteratorFactory> getGridAssignor () const noexcept
    {
        AMREX_ASSERT(m_defined);
        return AmrAssignGrid<BinIteratorFactory>(m_grid_assignors.dataPtr(), m_locators.size());
    }
};

}

#endif
