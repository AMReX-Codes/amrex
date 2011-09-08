      real function backtr(vrgs,nt,vr,vrg,zmin,zmax,ltail,ltpar,
     +                     utail,utpar)
c-----------------------------------------------------------------------
c
c           Back Transform Univariate Data from Normal Scores
c           *************************************************
c
c This subroutine backtransforms a standard normal deviate from a
c specified back transform table and option for the tails of the
c distribution.  Call once with "first" set to true then set to false
c unless one of the options for the tail changes.
c
c
c
c INPUT VARIABLES:
c
c   vrgs             normal score value to be back transformed
c   nt               number of values in the back transform tbale
c   vr(nt)           original data values that were transformed
c   vrg(nt)          the corresponding transformed values
c   zmin,zmax        limits possibly used for linear or power model
c   ltail            option to handle values less than vrg(1):
c   ltpar            parameter required for option ltail
c   utail            option to handle values greater than vrg(nt):
c   utpar            parameter required for option utail
c
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      parameter(EPSLON=1.0e-20)
      dimension vr(nt),vrg(nt)
      real*8      ltpar,utpar,lambda
      integer   ltail,utail
c
c Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
c
      if(vrgs.le.vrg(1)) then
            backtr = vr(1)
            cdflo  = gcum(vrg(1))
            cdfbt  = gcum(vrgs)
            if(ltail.eq.1) then
                  backtr = powint(0.d0,cdflo,zmin,vr(1),cdfbt,1.0)
            else if(ltail.eq.2) then
                  cpow   = 1.0 / ltpar
                  backtr = powint(0.d0,cdflo,zmin,vr(1),cdfbt,cpow)
            endif
c
c Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
c
      else if(vrgs.ge.vrg(nt)) then
            backtr = vr(nt)
            cdfhi  = gcum(vrg(nt))
            cdfbt  = gcum(vrgs)
            if(utail.eq.1) then
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,1.0)
            else if(utail.eq.2) then
                  cpow   = 1.0 / utpar
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,cpow)
            else if(utail.eq.4) then
                  lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
                  backtr = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
            endif
      else
c
c Value within the transformation table:
c
            call locate(vrg,nt,1,nt,vrgs,j)
            j = max(min((nt-1),j),1)
            backtr = powint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,1.0)
      endif
      return
      end
