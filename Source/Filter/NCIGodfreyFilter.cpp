#include <WarpX.H>
#include <NCIGodfreyFilter.H>
#include <NCIGodfreyTables.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

NCIGodfreyFilter::NCIGodfreyFilter(godfrey_coeff_set coeff_set_, amrex::Real cdtodz_, amrex::Real l_lower_order_in_v_){
    coeff_set = coeff_set_;
    cdtodz = cdtodz_;
    l_lower_order_in_v = l_lower_order_in_v_;
#if (AMREX_SPACEDIM == 3)
    stencil_length_each_dir = {1,1,5};
    slen = {1,1,5};
#else
    stencil_length_each_dir = {1,5};
    slen = {1,5,1};
#endif    
}

void NCIGodfreyFilter::ComputeStencils(){
    Print()<<"slen "<<slen<<std::endl;
    Print()<<"slen.x "<<slen.x<<std::endl;
    Print()<<"slen.y "<<slen.y<<std::endl;
    Print()<<"slen.z "<<slen.z<<std::endl;
    Print()<<"cdtodz "<<cdtodz<<std::endl;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
#if ( AMREX_SPACEDIM == 3 )
        slen.z==5,
#else
        slen.y==5, 
#endif
        "ERROR: NCI filter requires 5 points in z");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l_lower_order_in_v==1, 
        "ERROR: NCI corrector requires l_lower_order_in_v=1, i.e., Galerkin scheme");

    Real table_nci[tab_length][tab_width];
    if        (coeff_set == godfrey_coeff_set::Ex_Ey_Bz) {
        std::copy(&table_nci_godfrey_Ex_Ey_Bz[0][0], 
                  &table_nci_godfrey_Ex_Ey_Bz[0][0]+tab_width*tab_length, 
                  &table_nci[0][0]);
    } else if (coeff_set == godfrey_coeff_set::Bx_By_Ez) {
        std::copy(&table_nci_godfrey_Bx_By_Ez[0][0], 
                  &table_nci_godfrey_Bx_By_Ez[0][0]+tab_width*tab_length, 
                  &table_nci[0][0]);
    }
    
    int index = tab_length*cdtodz;
    index = min(index, tab_length-2);
    index = max(index, 0);
    Real weight_right = cdtodz - index/tab_length;
    Real prestencil[4];
    for(int i=0; i<tab_width; i++){
        Print()<<"i "<<i<<std::endl;
        prestencil[i] = (1-weight_right)*table_nci[index][i] + weight_right*table_nci[index+1][i];
    }
    stencil_z.resize( 5 );
    stencil_z[0] =  (256 + 128*prestencil[0] + 96*prestencil[1] + 80*prestencil[2] + 70*prestencil[3]) / 256;
    stencil_z[1] = -(       64*prestencil[0] + 64*prestencil[1] + 60*prestencil[2] + 56*prestencil[3]) / 256;
    stencil_z[2] =  (                          16*prestencil[1] + 24*prestencil[2] + 28*prestencil[3]) / 256;
    stencil_z[3] = -(                                              4*prestencil[2] +  8*prestencil[3]) / 256;
    stencil_z[4] =  (                                                                 1*prestencil[3]) / 256;
    stencil_z[0] /= 2.;

    stencil_x.resize(1);
    stencil_x[0] = 1.;
    stencil_x[0] /= 2.;
#if (AMREX_SPACEDIM == 3)
    stencil_y.resize(1);
    stencil_y[0] = 1.;
    stencil_y[0] /= 2.;
#endif
    
}
