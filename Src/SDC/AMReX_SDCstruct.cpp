#include "AMReX_SDCstruct.H"

/*
SDCstruct::SDCstruct()
      {
       int qtype_in=1;
       pf_quadrature(&qtype_in, &nnodes, &nnodes,qnodes,nflags, &qmats[0][0][0]);
     }
*/

SDCstruct::SDCstruct(int Nnodes_in,int Npieces_in, MultiFab& phi_in)
      {

       Nnodes=Nnodes_in;
       Npieces=Npieces_in;       


       int qtype_in=1;
       pf_quadrature(&qtype_in, &Nnodes, &Nnodes,qnodes,Nflags, &qmats[0][0][0]);
       
       phi_sdc.resize(Nnodes);
       res_sdc.resize(Nnodes);
       //       f_sdc.resize(Npieces);
       f_sdc.resize(Nnodes);       

       //       for (auto& v : f_sdc) v.resize(Nnodes);
       const BoxArray &ba=phi_in.boxArray();
       const DistributionMapping &dm=phi_in.DistributionMap();
       const int Nghost=phi_in.nGrow();
       const int Ncomp=phi_in.nComp();       
       for (int sdc_m = 0; sdc_m < Nnodes; sdc_m++)
	 {
	   phi_sdc[sdc_m].define(ba, dm, Ncomp, Nghost);
	   res_sdc[sdc_m].define(ba, dm, Ncomp, Nghost);
	   f_sdc[sdc_m].define(ba, dm, Ncomp, Nghost);	   
	   //  for (int i = 0; i < Npieces; i++)
	   //	     {
	   //	       f_sdc[i][sdc_m].define(ba, dm, Ncomp, Nghost);
	   //	     }
	 }
       
     }
