#ifndef AMREX_ARRAYOFSTRUCTS_H_
#define AMREX_ARRAYOFSTRUCTS_H_
#include <AMReX_Config.H>

#include <AMReX_Particle.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>

namespace amrex {

template <int NReal, int NInt,
          template<class> class Allocator=DefaultAllocator>
class ArrayOfStructs {
public:
    using ParticleType  = Particle<NReal, NInt>;
    using RealType      = typename ParticleType::RealType;

    using ParticleVector = amrex::PODVector<ParticleType, Allocator<ParticleType> >;

    using Iterator      = typename ParticleVector::iterator;
    using ConstIterator = typename ParticleVector::const_iterator;

    static constexpr int SizeInReal = sizeof(ParticleType) / sizeof(RealType);

    [[nodiscard]] const ParticleVector& operator() () const { return m_data; }
    [[nodiscard]]       ParticleVector& operator() ()       { return m_data; }

    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    [[nodiscard]] std::size_t size () const { return m_data.size(); }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    [[nodiscard]] int numParticles () const { return numRealParticles(); }

    /**
    * \brief Returns the number of real particles (excluding neighbors)
    *
    */
    [[nodiscard]] int numRealParticles () const { return numTotalParticles()-m_num_neighbor_particles; }

    /**
    * \brief Returns the number of neighbor particles (excluding reals)
    *
    */
    [[nodiscard]] int numNeighborParticles () const { return m_num_neighbor_particles; }

    /**
    * \brief Returns the total number of particles (real and neighbor)
    *
    */
    [[nodiscard]] int numTotalParticles () const { return m_data.size(); }

    void setNumNeighbors (int num_neighbors)
    {
        auto nrp = numRealParticles();
        m_num_neighbor_particles = num_neighbors;
        resize(nrp + num_neighbors);
    }

    [[nodiscard]] int getNumNeighbors () { return m_num_neighbor_particles; }

    [[nodiscard]] bool empty () const { return m_data.empty(); }

    [[nodiscard]] const RealType* data () const { return &(m_data[0].m_pos[0]); }
    [[nodiscard]]       RealType* data ()       { return &(m_data[0].m_pos[0]); }

    [[nodiscard]] const RealType* dataPtr () const { return data(); }
    [[nodiscard]]       RealType* dataPtr ()       { return data(); }

    [[nodiscard]] std::pair<int,int> dataShape () const {
        return std::make_pair(SizeInReal, static_cast<int>(m_data.size()));
    }

    void push_back (const ParticleType& p) { return m_data.push_back(p); }
    void pop_back() {m_data.pop_back(); }
    [[nodiscard]] bool empty() {return m_data.empty(); }

    [[nodiscard]] const ParticleType& back() const {return m_data.back(); }
    [[nodiscard]]       ParticleType& back()       {return m_data.back(); }

    [[nodiscard]] const ParticleType& operator[] (int i) const { return m_data[i]; }
    [[nodiscard]]       ParticleType& operator[] (int i)       { return m_data[i]; }

    void swap (ArrayOfStructs<NReal, NInt>& other)
    {
        m_data.swap(other.m_data);
    }

    void resize (size_t count) { m_data.resize(count); }

    Iterator erase ( ConstIterator first, ConstIterator second) { return m_data.erase(first, second); }

    template< class InputIt >
    void insert ( Iterator pos, InputIt first, InputIt last ) { m_data.insert(pos, first, last); }

    [[nodiscard]] typename ParticleVector::iterator begin () { return m_data.begin(); }
    [[nodiscard]] typename ParticleVector::const_iterator begin () const { return m_data.begin(); }
    [[nodiscard]] typename ParticleVector::const_iterator cbegin () const { return m_data.cbegin(); }

    [[nodiscard]] typename ParticleVector::iterator end () { return m_data.end(); }
    [[nodiscard]] typename ParticleVector::const_iterator end () const { return m_data.end(); }
    [[nodiscard]] typename ParticleVector::const_iterator cend () const { return m_data.cend(); }

    int m_num_neighbor_particles{0};

private:
    ParticleVector m_data;
};

#if __cplusplus < 201703L
template <int NReal, int NInt,
          template<class> class Allocator>
constexpr int ArrayOfStructs<NReal, NInt, Allocator>::SizeInReal;
#endif

} // namespace amrex

#endif // AMREX_ARRAYOFSTRUCTS_H_
