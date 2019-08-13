
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define CKRP CKRP
#define CKPY CKPY
#define CKRHOY CKRHOY
#define CKWT CKWT
#define CKMMWY CKMMWY
#define CKYTX CKYTX
#define CKXTY CKXTY
#define CKCPBS CKCPBS
#define CKHMS CKHMS
#define CKCVBL CKCVBL
#define CKCVBS CKCVBS
#define CKUBMS CKUBMS
#define VCKWYR VCKWYR
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
#elif defined(BL_FORT_USE_LOWERCASE)
#define CKRP ckrp
#define CKPY ckpy
#define CKRHOY ckrhoy
#define CKWT ckwt
#define CKMMWY ckmmwy
#define CKYTX ckytx
#define CKXTY ckxty
#define CKCPBS ckcpbs
#define CKHMS ckhms
#define CKCVBL ckcvbl
#define CKCVBS ckcvbs
#define CKUBMS ckubms
#define VCKWYR vckwyr
#define GET_T_GIVEN_EY get_t_given_ey
#elif defined(BL_FORT_USE_UNDERSCORE)
#define CKRP ckrp_
#define CKPY ckpy_
#define CKRHOY ckrhoy_
#define CKWT ckwt_
#define CKMMWY ckmmwy_
#define CKYTX ckytx_
#define CKXTY ckxty_
#define CKCPBS ckcpbs_
#define CKHMS ckhms_
#define CKCVBL ckcvbl_
#define CKCVBS ckcvbs_
#define CKUBMS ckubms_
#define VCKWYR vckwyr_
#define GET_T_GIVEN_EY get_t_given_ey_
#endif

/*function declarations */
void molecularWeight(double * restrict wt);
void gibbs(double * restrict species, double * restrict tc);
void speciesInternalEnergy(double * restrict species, double * restrict tc);
void speciesEnthalpy(double * restrict species, double * restrict tc);
void cp_R(double * restrict species, double * restrict tc);
void cv_R(double * restrict species, double * restrict tc);
void vproductionRate(int npt, double * restrict wdot, double * restrict c, double * restrict T);
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa);
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P);
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho);
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt);
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm);
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x);
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y);
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict ums);
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cpbl);
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs);
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms);
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
            double * restrict y, int * restrict iwrk, double * restrict rwrk,
            double * restrict wdot);
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int *ierr);

/*
** Inverse of molecular weights.
** So we can multiply instead of divide in places.
*/
static const double imw[9] = { 1.0 / 2.01594,   /*H2 */
                               1.0 / 31.9988,   /*O2 */
                               1.0 / 18.01534,  /*H2O */
                               1.0 / 1.00797,   /*H */
                               1.0 / 15.9994,   /*O */
                               1.0 / 17.00737,  /*OH */
                               1.0 / 33.00677,  /*HO2 */
                               1.0 / 34.01474,  /*H2O2 */
                               1.0 / 28.0134 }; /*N2 */


/* Returns R, Rc, Patm */
void CKRP(int * ickwrk, double * restrict rckwrk, double * restrict ru, double * restrict ruc, double * restrict pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(y) */
void CKPY(double * restrict rho, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double * restrict P, double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict rho)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (83145100 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*get molecular weight for all species */
void CKWT(int * iwrk, double * restrict rwrk, double * restrict wt)
{
    molecularWeight(wt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict wtm)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double * restrict y, int * iwrk, double * restrict rwrk, double * restrict x)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 9; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double * restrict x, int * iwrk, double * restrict rwrk, double * restrict y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*31.998800/XW; 
    y[2] = x[2]*18.015340/XW; 
    y[3] = x[3]*1.007970/XW; 
    y[4] = x[4]*15.999400/XW; 
    y[5] = x[5]*17.007370/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*28.013400/XW; 

    return;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9], tresult[9]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 9; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 9; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double * restrict T, int * iwrk, double * restrict rwrk, double * restrict hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 9; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*H2O */
    result += cvor[3]*y[3]*imw[3]; /*H */
    result += cvor[4]*y[4]*imw[4]; /*O */
    result += cvor[5]*y[5]*imw[5]; /*OH */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double * restrict T, double * restrict x, int * iwrk, double * restrict rwrk, double * restrict cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*get mean internal energy in mass units */
void CKUBMS(double * restrict T, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[9]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*H2O */
    result += y[3]*ums[3]*imw[3]; /*H */
    result += y[4]*ums[4]*imw[4]; /*O */
    result += y[5]*ums[5]*imw[5]; /*OH */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

    *ubms = result * RT;
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int * restrict np, double * restrict rho, double * restrict T,
	    double * restrict y, int * restrict iwrk, double * restrict rwrk,
	    double * restrict wdot)
{
    double c[9*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<9; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<9*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*save molecular weights into array */
void molecularWeight(double * restrict wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 31.998800; /*O2 */
    wt[2] = 18.015340; /*H2O */
    wt[3] = 1.007970; /*H */
    wt[4] = 15.999400; /*O */
    wt[5] = 17.007370; /*OH */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 34.014740; /*H2O2 */
    wt[8] = 28.013400; /*N2 */

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double * restrict wdot, double * restrict sc, double * restrict T)
{
    double k_f_s[21][npt], Kc_s[21][npt], mixture[npt], g_RT[9*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0][i] = 1e-06 * 3.547e+15*exp(-0.406*tc[i]-8352.8934356925419706*invT[i]);
        k_f_s[1][i] = 1e-06 * 50800*exp(2.67*tc[i]-3165.2328279116868543*invT[i]);
        k_f_s[2][i] = 1e-06 * 2.16e+08*exp(1.51*tc[i]-1726.0331637101885462*invT[i]);
        k_f_s[3][i] = 1e-06 * 2.97e+06*exp(2.02*tc[i]-6743.1033217832446098*invT[i]);
        k_f_s[4][i] = 1e-06 * 4.577e+19*exp(-1.4*tc[i]-52525.75557669664704*invT[i]);
        k_f_s[5][i] = 1e-12 * 6.165e+15*exp(-0.5*tc[i]);
        k_f_s[6][i] = 1e-12 * 4.714e+18*exp(-1*tc[i]);
        k_f_s[7][i] = 1e-12 * 3.8e+22*exp(-2*tc[i]);
        k_f_s[8][i] = 1e-06 * 1.475e+12*exp(0.6*tc[i]);
        k_f_s[9][i] = 1e-06 * 1.66e+13*exp(-414.14731595728432012*invT[i]);
        k_f_s[10][i] = 1e-06 * 7.079e+13*exp(-148.44891641239232172*invT[i]);
        k_f_s[11][i] = 1e-06 * 3.25e+13;
        k_f_s[12][i] = 1e-06 * 2.89e+13*exp(+250.09868290494571852*invT[i]);
        k_f_s[13][i] = 1e-06 * 4.2e+14*exp(-6029.5420896721516328*invT[i]);
        k_f_s[14][i] = 1e-06 * 1.3e+11*exp(+819.89091359562974048*invT[i]);
        k_f_s[15][i] = 1 * 2.951e+14*exp(-24370.783124922574643*invT[i]);
        k_f_s[16][i] = 1e-06 * 2.41e+13*exp(-1997.7701632447369775*invT[i]);
        k_f_s[17][i] = 1e-06 * 4.82e+13*exp(-4000.5724931475210724*invT[i]);
        k_f_s[18][i] = 1e-06 * 9.55e+06*exp(2*tc[i]-1997.7701632447369775*invT[i]);
        k_f_s[19][i] = 1e-06 * 1e+12;
        k_f_s[20][i] = 1e-06 * 5.8e+14*exp(-4809.2416750957054319*invT[i]);
    }

    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[9];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = 101325. / 8.31451 / T[i];

        Kc_s[0][i] = exp((g_RT[3*npt+i] + g_RT[1*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[1][i] = exp((g_RT[4*npt+i] + g_RT[0*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[2][i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[3][i] = exp((g_RT[4*npt+i] + g_RT[2*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[4][i] = refC * exp((g_RT[0*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[5][i] = 1.0 / (refC) * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i]));
        Kc_s[6][i] = 1.0 / (refC) * exp((g_RT[4*npt+i] + g_RT[3*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[7][i] = 1.0 / (refC) * exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[8][i] = 1.0 / (refC) * exp((g_RT[3*npt+i] + g_RT[1*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[9][i] = exp((g_RT[6*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[1*npt+i]));
        Kc_s[10][i] = exp((g_RT[6*npt+i] + g_RT[3*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[11][i] = exp((g_RT[6*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[12][i] = exp((g_RT[6*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[1*npt+i]));
        Kc_s[13][i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[7*npt+i] + g_RT[1*npt+i]));
        Kc_s[14][i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[7*npt+i] + g_RT[1*npt+i]));
        Kc_s[15][i] = refC * exp((g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[5*npt+i]));
        Kc_s[16][i] = exp((g_RT[7*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[17][i] = exp((g_RT[7*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i] + g_RT[0*npt+i]));
        Kc_s[18][i] = exp((g_RT[7*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[19][i] = exp((g_RT[7*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i] + g_RT[2*npt+i]));
        Kc_s[20][i] = exp((g_RT[7*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i] + g_RT[2*npt+i]));
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<9; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha, redP, F, logPred, logFcent;
        double troe_c, troe_n, troe, F_troe;

        /*reaction 1: H + O2 <=> O + OH */
        phi_f = sc[3*npt+i]*sc[1*npt+i];
        k_f = k_f_s[0][i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[0][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[1*npt+i] -= 1 * qdot;
        wdot[4*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 2: O + H2 <=> H + OH */
        phi_f = sc[4*npt+i]*sc[0*npt+i];
        k_f = k_f_s[1][i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[1][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 3: H2 + OH <=> H2O + H */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[2][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[2][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;

        /*reaction 4: O + H2O <=> OH + OH */
        phi_f = sc[4*npt+i]*sc[2*npt+i];
        k_f = k_f_s[3][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[3][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[2*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 5: H2 + M <=> H + H + M */
        phi_f = sc[0*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[4][i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[4][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;
        wdot[3*npt+i] += 1 * qdot;

        /*reaction 6: O + O + M <=> O2 + M */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[5][i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i];
        Kc = Kc_s[5][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 7: O + H + M <=> OH + M */
        phi_f = sc[4*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[6][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[6][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 8: H + OH + M <=> H2O + M */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = alpha * k_f_s[7][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[7][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;

        /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[3*npt+i]*sc[1*npt+i];
        alpha = mixture[i] + sc[0*npt+i] + 10*sc[2*npt+i] + -0.22*sc[1*npt+i];
        k_f = k_f_s[8][i];
        redP = 1e-12 * alpha / k_f * 6.366e+20*exp(-1.72*tc[i]-264.08810621431689469*invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10((0.2*exp(T[i]/-1e-30))+ (0.8*exp(T[i]/-1e+30)));
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[8][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[1*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;

        /*reaction 10: HO2 + H <=> H2 + O2 */
        phi_f = sc[6*npt+i]*sc[3*npt+i];
        k_f = k_f_s[9][i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[1*npt+i];
        Kc = Kc_s[9][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[0*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 11: HO2 + H <=> OH + OH */
        phi_f = sc[6*npt+i]*sc[3*npt+i];
        k_f = k_f_s[10][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[10][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 12: HO2 + O <=> O2 + OH */
        phi_f = sc[6*npt+i]*sc[4*npt+i];
        k_f = k_f_s[11][i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[11][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 13: HO2 + OH <=> H2O + O2 */
        phi_f = sc[6*npt+i]*sc[5*npt+i];
        k_f = k_f_s[12][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[1*npt+i];
        Kc = Kc_s[12][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[13][i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[1*npt+i];
        Kc = Kc_s[13][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[7*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[14][i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[1*npt+i];
        Kc = Kc_s[14][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[6*npt+i] -= 1 * qdot;
        wdot[7*npt+i] += 1 * qdot;
        wdot[1*npt+i] += 1 * qdot;

        /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + 1.5*sc[0*npt+i] + 11*sc[2*npt+i];
        k_f = k_f_s[15][i];
        redP = 1e-6 * alpha / k_f * 1.202e+17*exp(-22896.358294114746968*invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10((0.5*exp(T[i]/-1e-30))+ (0.5*exp(T[i]/-1e+30)));
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[5*npt+i];
        Kc = Kc_s[15][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 17: H2O2 + H <=> H2O + OH */
        phi_f = sc[7*npt+i]*sc[3*npt+i];
        k_f = k_f_s[16][i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[16][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;

        /*reaction 18: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[7*npt+i]*sc[3*npt+i];
        k_f = k_f_s[17][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[0*npt+i];
        Kc = Kc_s[17][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[3*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[0*npt+i] += 1 * qdot;

        /*reaction 19: H2O2 + O <=> OH + HO2 */
        phi_f = sc[7*npt+i]*sc[4*npt+i];
        k_f = k_f_s[18][i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[18][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[4*npt+i] -= 1 * qdot;
        wdot[5*npt+i] += 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;

        /*reaction 20: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[7*npt+i]*sc[5*npt+i];
        k_f = k_f_s[19][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[2*npt+i];
        Kc = Kc_s[19][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;

        /*reaction 21: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[7*npt+i]*sc[5*npt+i];
        k_f = k_f_s[20][i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[2*npt+i];
        Kc = Kc_s[20][i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= 1 * qdot;
        wdot[5*npt+i] -= 1 * qdot;
        wdot[6*npt+i] += 1 * qdot;
        wdot[2*npt+i] += 1 * qdot;
    }
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.012520870000000e+03 * invT
            +6.592218400000000e+00
            -3.298124310000000e+00 * tc[0]
            -4.124720870000000e-04 * tc[1]
            +1.357169215000000e-07 * tc[2]
            +7.896195275000000e-12 * tc[3]
            -2.067436120000000e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.005249020000000e+03 * invT
            -2.821801190000000e+00
            -3.212936400000000e+00 * tc[0]
            -5.637431750000000e-04 * tc[1]
            +9.593584116666666e-08 * tc[2]
            -1.094897691666667e-10 * tc[3]
            +4.384276960000000e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.020811330000000e+04 * invT
            +7.966096399999998e-01
            -3.386842490000000e+00 * tc[0]
            -1.737491230000000e-03 * tc[1]
            +1.059116055000000e-06 * tc[2]
            -5.807151058333333e-10 * tc[3]
            +1.253294235000000e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117608000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.914764450000000e+04 * invT
            -1.756619999999964e-02
            -2.946428780000000e+00 * tc[0]
            +8.190832450000000e-04 * tc[1]
            -4.035052833333333e-07 * tc[2]
            +1.335702658333333e-10 * tc[3]
            -1.945348180000000e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766314650000000e+04 * invT
            -3.396609550000000e+00
            -3.388753650000000e+00 * tc[0]
            -3.284612905000000e-03 * tc[1]
            +2.475020966666667e-08 * tc[2]
            +3.854837933333333e-10 * tc[3]
            -1.235757375000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.350339970000000e+02 * invT
            +4.346533540000000e+00
            -2.991423370000000e+00 * tc[0]
            -3.500322055000000e-04 * tc[1]
            +9.389714483333333e-09 * tc[2]
            +7.692981816666667e-13 * tc[3]
            -7.913758950000000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.233930180000000e+03 * invT
            +5.084126000000002e-01
            -3.697578190000000e+00 * tc[0]
            -3.067598445000000e-04 * tc[1]
            +2.098069983333333e-08 * tc[2]
            -1.479401233333333e-12 * tc[3]
            +5.682176550000000e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.989920900000000e+04 * invT
            -4.190671200000001e+00
            -2.672145610000000e+00 * tc[0]
            -1.528146445000000e-03 * tc[1]
            +1.455043351666667e-07 * tc[2]
            -1.000830325000000e-11 * tc[3]
            +3.195808935000000e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117638000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.923080270000000e+04 * invT
            -2.378248450000000e+00
            -2.542059660000000e+00 * tc[0]
            +1.377530955000000e-05 * tc[1]
            +5.171338916666667e-10 * tc[2]
            -3.792556183333333e-13 * tc[3]
            +2.184025750000000e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696090000000e+04 * invT
            +4.072029891000000e+00
            -4.573166850000000e+00 * tc[0]
            -2.168068195000000e-03 * tc[1]
            +2.457814700000000e-07 * tc[2]
            -1.957419641666667e-11 * tc[3]
            +7.158267800000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double * restrict species, double * restrict tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double * restrict e, double * restrict y, int * iwrk, double * restrict rwrk, double * restrict t, int * ierr)
{
    const int maxiter = 200;
    const double tol  = 1.e-6;
    double ein  = *e;
    double tmin = 250;/*max lower bound for thermo def */
    double tmax = 3500;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, iwrk, rwrk, &emin);
    CKUBMS(&tmax, y, iwrk, rwrk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, iwrk, rwrk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, iwrk, rwrk, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,iwrk,rwrk,&e1);
        CKCVBS(&t1,y,iwrk,rwrk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

