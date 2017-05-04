#ifndef AMREX_BUILDINFO_H_
#define AMREX_BUILDINFO_H_

namespace amrex {
const char* buildInfoGetBuildDate();
const char* buildInfoGetBuildDir();
const char* buildInfoGetBuildMachine();
const char* buildInfoGetAMReXDir();
const char* buildInfoGetComp();
const char* buildInfoGetCompVersion();
const char* buildInfoGetCompilerLoc(int i);
const char* buildInfoGetCompilerVer(int i);
const char* buildInfoGetCompilerFlags(int i);
const char* buildInfoGetFcomp();
const char* buildInfoGetFcompVersion();
const char* buildInfoGetCXXName();
const char* buildInfoGetFName();
const char* buildInfoGetCXXFlags();
const char* buildInfoGetFFlags();
const char* buildInfoGetLinkFlags();
const char* buildInfoGetLibraries();
const char* buildInfoGetAux(int i);
int buildInfoGetNumModules();
const char* buildInfoGetModuleName(int i);
const char* buildInfoGetModuleVal(int i);
const char* buildInfoGetDif(int i);
const char* buildInfoGetGitHash(int i);
const char* buildInfoGetGitBranch(int i);
const char* buildInfoGetLinkerInfo(int i);
const char* buildInfoGetBuildGitName();
const char* buildInfoGetBuildGitHash();
}

#endif
