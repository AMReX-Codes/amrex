#ifndef WARPX_FULLDIAGNOSTICS_H_
#define WARPX_FULLDIAGNOSTICS_H_

#include "Diagnostics.H"

class FullDiagnostics : public Diagnostics
{
public:
    FullDiagnostics (int i, std::string name);
};

#endif // WARPX_FULLDIAGNOSTICS_H_
