/*

  This program compares the particle data stored in plt files using the 
  Checkpoint() method of the ParticleContainer. 

  To compile, navigate to AMREX_HOME/Tools/Postprocessing/C_Src and type
  "make".

  Usage:

      mpirun -np 4 ./particle_compare.exe old00000 new00000 Tracer

  This compares the particle type "Tracer" between the old00000 and 
  new00000 plt files. For this to work, the plt files must have been
  run with same grids / number of processes.

 */

#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdlib>

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

///
/// This stores the metadata associated with an AMReX particle header file.
///
struct ParticleHeader {
    
    ///
    /// Initializes the struct to store the header information associated with 
    /// the given particle header file. For example, if plt_file = "plt00000"
    /// and particle_type = "Tracer", this will read the Header file at
    /// "plt00000/Tracer/Header"
    ///
    ParticleHeader(const std::string& plt_file, 
                   const std::string& particle_type) {
        
        plt_file_name = plt_file;
        par_file_name = plt_file_name + "/" + particle_type;
        hdr_file_name = par_file_name + "/Header";

        std::ifstream file(hdr_file_name.c_str());
        if ( not file.is_open() ) {
#ifdef BL_USE_MPI
            int myproc;
            MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#else
            int myproc = 0;
#endif
            if ( myproc == 0) {
                std::cout << "Error! Could not open file " << hdr_file_name << "." << std::endl;
            }
#ifdef USE_MPI
            MPI_Finalize();
#endif 
            exit(1);
        }

        file >> *this;
        if (not file.is_open() ) {
#ifdef BL_USE_MPI
            int myproc;
            MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#else
            int myproc = 0;
#endif
            if ( myproc == 0) {
                std::cout << "File header is corrupt." << std::endl;
            }
#ifdef USE_MPI
            MPI_Finalize();
#endif 
            exit(1);
        }

        file.close();
        
        num_real = ndim + num_real_extra;
        num_int  = 2*is_checkpoint + num_int_extra;
        num_comp = num_real + num_int;
        
        std::string directions = "xyz";    
        for (int i = 0; i < ndim; ++i) {
            std::stringstream ss;
            ss << particle_type << "_position_" << directions[i];
            comp_names.push_back(ss.str());
        }
        
        for (int i = 0; i < num_real_extra; ++i) {
            std::stringstream ss;
            ss << particle_type << "_" << real_comps[i];
            comp_names.push_back(ss.str());
        }
        
        if (is_checkpoint) {
            {
                std::stringstream ss;
                ss << particle_type << "_id";
                comp_names.push_back(ss.str());
            }
            {
                std::stringstream ss;
                ss << particle_type << "_cpu";
                comp_names.push_back(ss.str());
            }
        }
        
        for (int i = 0; i < num_int_extra; ++i) {
            std::stringstream ss;
            ss << particle_type << "_" << int_comps[i];
            comp_names.push_back(ss.str());
        }
    }
    
    // This is the metadata actually stored in the header files
    std::string version;
    int ndim;
    int num_real_extra;
    std::vector<std::string> real_comps;
    int num_int_extra;
    std::vector<std::string> int_comps;
    bool is_checkpoint;
    int nparticles;
    int next_id;
    int finest_level;
    std::vector<int> num_grids;
    std::vector<std::vector<int> > file_nums;
    std::vector<std::vector<int> > particle_counts;
    std::vector<std::vector<int> > offsets;

    // This is additional metadata derived from the above
    int num_real;
    int num_int;
    int num_comp;
    std::string plt_file_name;
    std::string par_file_name;
    std::string hdr_file_name;
    std::vector<std::string> comp_names;
    
    // These operators only use the data actually in the headers, not the derived stuff.
    friend std::ostream& operator<< (std::ostream& stream, const ParticleHeader& header);
    friend std::istream& operator>> (std::istream& stream, ParticleHeader& header);
    friend bool operator==(const ParticleHeader& lhs, const ParticleHeader& rhs);
    friend bool operator!=(const ParticleHeader& lhs, const ParticleHeader& rhs);
};

bool operator==(const ParticleHeader& lhs, const ParticleHeader& rhs) {    
    return lhs.version         == rhs.version &&
           lhs.ndim            == rhs.ndim &&
           lhs.num_real_extra  == rhs.num_real_extra &&
           lhs.real_comps      == rhs.real_comps &&
           lhs.num_int_extra   == rhs.num_int_extra &&
           lhs.int_comps       == rhs.int_comps &&
           lhs.is_checkpoint   == rhs.is_checkpoint &&
           lhs.nparticles      == rhs.nparticles &&
           lhs.next_id         == rhs.next_id && 
           lhs.finest_level    == rhs.finest_level &&
           lhs.num_grids       == rhs.num_grids &&
           lhs.particle_counts == rhs.particle_counts &&
           lhs.offsets         == rhs.offsets;
}

bool operator!=(const ParticleHeader& lhs, const ParticleHeader& rhs) {
    return not (lhs == rhs);
}

std::ostream& operator<< (std::ostream& stream, const ParticleHeader& header) {

    stream << header.version << std::endl;
    stream << header.ndim << std::endl;

    stream << header.num_real_extra << std::endl;
    for (int i = 0; i < header.num_real_extra; ++i) {
        stream << header.real_comps[i] << std::endl;
    }

    stream << header.num_int_extra << std::endl;
    for (int i = 0; i < header.num_int_extra; ++i) {
        stream << header.int_comps[i] << std::endl;
    }

    stream << header.is_checkpoint << std::endl;
    stream << header.nparticles << std::endl;
    stream << header.next_id << std::endl;
    stream << header.finest_level << std::endl;

    for (int i = 0; i <= header.finest_level; ++i) {
        stream << header.num_grids[i] << std::endl;
    }
  
    for (int j = 0; j <= header.finest_level; ++j) {
        for (unsigned i = 0; i < header.file_nums.size(); ++i) {
            stream << header.file_nums[j][i] << " ";
            stream << header.particle_counts[j][i] << " ";
            stream << header.offsets[j][i] << std::endl;
        }
    }

    return stream;
}

std::istream& operator>> (std::istream& stream, ParticleHeader& header) {

    stream >> header.version;
    stream >> header.ndim;

    stream >> header.num_real_extra;
    for (int i = 0; i < header.num_real_extra; ++i) {
        std::string comp;
        stream >> comp;
        header.real_comps.push_back(comp);
    }

    stream >> header.num_int_extra;
    for (int i = 0; i < header.num_int_extra; ++i) {
        std::string comp;
        stream >> comp;
        header.int_comps.push_back(comp);
    }

    stream >> header.is_checkpoint;
    stream >> header.nparticles;
    stream >> header.next_id;
    stream >> header.finest_level;

    for (int i = 0; i <= header.finest_level; ++i) {
        int grid_count;
        stream >> grid_count;
        header.num_grids.push_back(grid_count);
    }

    header.file_nums.resize(header.finest_level+1);
    header.particle_counts.resize(header.finest_level+1);
    header.offsets.resize(header.finest_level+1);

    for (int i = 0; i <= header.finest_level; ++i) {
        for (int j = 0; j < header.num_grids[i]; ++j) {
            int num;
            stream >> num;
            header.file_nums[i].push_back(num);
            stream >> num;
            header.particle_counts[i].push_back(num);
            stream >> num;
            header.offsets[i].push_back(num);
        }
    }

    return stream;
}

std::string getDataFileName(const std::string& prefix, int level, int file_num) {
    std::stringstream ss;
    ss << prefix << "/Level_" << level << "/DATA_";
    ss << std::setfill('0');
    ss << std::setw(5);
    ss << file_num;
    return ss.str();
}

void getDataBuffer(std::vector<char>& buffer, const std::string& file,
                   size_t buffer_size, int offset) {
    std::ifstream is(file.c_str(), std::ifstream::binary);
    assert(is.is_open());
    is.seekg(offset);
    is.read(buffer.data(), buffer_size);
    is.close();
}

void compare_particle_chunk(const ParticleHeader& header1,
                            const ParticleHeader& header2,
                            std::vector<double>&  norms,
                            int level, int file_num, int np, int offset) {
    
    if (np == 0) return;

    std::string read_file1 = getDataFileName(header1.par_file_name, level, file_num);
    std::string read_file2 = getDataFileName(header2.par_file_name, level, file_num);

    int idata_size = header1.num_int*sizeof(int);
    int rdata_size = header1.num_real*sizeof(double);
    int pdata_size = rdata_size + idata_size;    
    size_t buffer_size = pdata_size * np;

    std::vector<char> data1(buffer_size);
    getDataBuffer(data1, read_file1, buffer_size, offset);

    std::vector<char> data2(buffer_size);
    getDataBuffer(data2, read_file2, buffer_size, offset);

    char* tmp1 = data1.data();
    char* tmp2 = data2.data();
    for (int i = 0; i < np; ++i) {
        for (int j = 0; j < header1.num_int; ++j) {
            int val1, val2;
            std::memcpy(&val1, tmp1, sizeof(int));
            std::memcpy(&val2, tmp2, sizeof(int));
            norms[j+header1.num_real] = std::max((double) std::abs(val2 - val1), 
                                                 norms[j+header1.num_real]);
            if (val1 == 0) {
                norms[header1.num_comp+j+header1.num_real] = 
                    norms[j+header1.num_real];
            } else {
                norms[header1.num_comp+j+header1.num_real] = 
                    norms[j+header1.num_real] / std::abs(val1);
            }
            tmp1 += sizeof(int);
            tmp2 += sizeof(int);
        }
    }
    
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < header1.num_real; ++j) {
            double val1, val2;
            std::memcpy(&val1, tmp1, sizeof(double));
            std::memcpy(&val2, tmp2, sizeof(double));
            norms[j] = std::max(std::abs(val2 - val1), norms[j]);
            if (val1 == 0) {
                norms[header1.num_comp+j] = norms[j];
            } else {
                norms[header1.num_comp+j] = norms[j] / std::abs(val1);
            }
            tmp1 += sizeof(double);
            tmp2 += sizeof(double);
        }
    }
}

int main(int argc, char* argv[])
{

#ifdef BL_USE_MPI
    MPI_Init(NULL, NULL);

    int nprocs, myproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
#else
    int nprocs = 1;
    int myproc = 0;
#endif

    if (argc < 4) {
        if (myproc == 0)
            std::cerr << "Usage: " << argv[0] << " file1 file2 ptype" << std::endl;
#ifdef BL_USE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    std::string fn1 = argv[1];
    std::string fn2 = argv[2];
    std::string pt  = argv[3];
    
    ParticleHeader header1(fn1, pt);
    ParticleHeader header2(fn2, pt);
    
    if (header1 != header2) {
        if (myproc == 0) 
            std::cout << "FAIL - Particle data headers do not agree." << std::endl;
    }

    // for each grid, store the corresponding information about where to look up the 
    // particle data
    std::vector<int> levels;
    std::vector<int> file_nums;
    std::vector<int> particle_counts;
    std::vector<int> offsets;
    for (int lev = 0; lev <= header1.finest_level; ++lev) {
        for (unsigned gid = 0; gid < header1.file_nums[lev].size(); ++gid) {
            levels.push_back(lev);
            file_nums.push_back(header1.file_nums[lev][gid]);
            particle_counts.push_back(header1.particle_counts[lev][gid]);
            offsets.push_back(header1.offsets[lev][gid]);
        }
    }

    // cut the list of grids based on the number of MPI tasks
    int n = file_nums.size();
    int ibegin, iend;
    int navg = n/nprocs;
    int nleft = n - navg * nprocs;
    if (myproc < nleft) {
        ibegin = myproc*(navg+1);
        iend = ibegin + navg+1;
    } else {
        ibegin = myproc*navg + nleft;
        iend = ibegin + navg;
    }

#ifdef BL_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Each proc computes the max norm of the particle data diff over its grids    
    // The first num_comp values are the abs norms, the second are the rel norms
    std::vector<double> norms(2*header1.num_comp, 0.0);
    for (int i = ibegin; i < iend; ++i) {
        compare_particle_chunk(header1, header2, norms, 
                               levels[i], file_nums[i], 
                               particle_counts[i], offsets[i]);
    }

#ifdef BL_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    std::vector<double> global_norms(2*header1.num_comp, 0.0);
#ifdef BL_USE_MPI
    // parallel reduction
    MPI_Reduce(norms.data(), global_norms.data(), 2*header1.num_comp, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
#else
    for (unsigned i = 0; i < global_norms.size(); ++i) {
        global_norms[i] = norms[i];
    }
#endif

    // write out results
    if (myproc == 0) {
        std::cout << std::endl << std::string(71, '-') << std::endl;
        std::cout << pt << std::endl;
        for (int i = 0; i < header1.num_comp; ++i) {
            std::cout << std::scientific << std::left << std::setw(36) << std::setfill(' ') << std::setprecision(8) << header1.comp_names[i];
            std::cout << std::scientific << std::left << std::setw(22) << std::setfill(' ') << std::setprecision(8) << global_norms[i];
            std::cout << std::scientific << std::left << std::setw(22) << std::setfill(' ') << std::setprecision(8) << global_norms[i+header1.num_comp];
            std::cout << std::endl;
        }
       std::cout << std::endl;
    }

#ifdef BL_USE_MPI
    MPI_Finalize();
#endif

    int exit_code = 0;
    for (unsigned i = 0; i < global_norms.size(); ++i) {
        if (global_norms[i] > 0.0) exit_code = 1;
    }

    return exit_code;
}
