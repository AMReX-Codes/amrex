#ifndef PLOTFILEUTIL_H_
#define PLOTFILEUTIL_H_

#include <vector>
#include <string>

#include <Geometry.H>
#include <MultiFab.H>

namespace BoxLib
{
    void WriteSingleLevelPlotfile (const std::string& plotfilename,
				   const MultiFab& mf, const std::vector<std::string>& varnames,
				   const Geometry& geom, Real t, int level_step);

    void WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
				  const Array<const MultiFab*>& mf,
				  const std::vector<std::string>& varnames,
				  const Array<Geometry>& geom, Real time, const Array<int>& level_steps,
				  const Array<IntVect>& ref_ratio);
}

#endif
