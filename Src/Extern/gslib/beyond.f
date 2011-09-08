      subroutine beyond(ivtype,nccut,ccut,ccdf,ncut,cut,cdf,zmin,zmax,
     +           ltail,ltpar,middle,mpar,utail,utpar,zval,cdfval,ierr)
c-----------------------------------------------------------------------
c
c                     Go Beyond a Discrete CDF
c                     ************************
c
c This subroutine is a general purpose subroutine to interpolate within
c and extrapolate beyond discrete points on a conditional CDF.  If the
c Z value "zval" is specified then the corresponding CDF value "cdfval"
c will be computed, if the CDF value "cdfval" is specified the
c corresponding Z value "zval" will be computed.
c
c
c
c INPUT/OUTPUT VARIABLES:
c
c   ivtype           variable type (1=continuous, 0=categorical)
c   nccut            number of cutoffs defining the conditional CDF
c   ccut()           real array of the nccut cutoffs
c   ccdf()           real array of the conditional cdf values
c   ncut             number of cutoffs defining the global CDF
c   cut()            real array of the ncut cutoffs
c   cdf()            real array of the global cdf values
c
c   zmin,zmax        minimum and maximum allowable data values
c   ltail            option to handle values in lower tail
c   ltpar            parameter required for option ltail
c   middle           option to handle values in the middle
c   mpar             parameter required for option middle
c   utail            option to handle values in upper tail
c   utpar            parameter required for option utail
c
c   zval             interesting cutoff (if -1 then it is calculated)
c   cdfval           interesting CDF (if -1 then it is calculated)
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      parameter(EPSLON=1.0e-20,UNEST=-1.0)
      dimension ccut(nccut),ccdf(nccut),cut(1),cdf(1)
      real*8      utpar,mpar,ltpar,lambda
      integer   ltail,utail,middle,cclow,cchigh
c
c Check for both "zval" and "cdfval" defined or undefined:
c
      ierr  = 1
      if(zval.gt.UNEST.and.cdfval.gt.UNEST) return
      if(zval.le.UNEST.and.cdfval.le.UNEST) return
c
c Handle the case of a categorical variable:
c
      if(ivtype.eq.0) then
            cum = 0
            do i=1,nccut
                  cum = cum + ccdf(i)
                  if(cdfval.le.cum) then
                         zval = ccut(i)
                         return
                  endif
            end do
            return
      end if
c
c Figure out what part of distribution: ipart = 0 - lower tail
c                                       ipart = 1 - middle
c                                       ipart = 2 - upper tail
      ierr  = 0
      ipart = 1
      if(zval.gt.UNEST) then
            if(zval.le.ccut(1))       ipart = 0
            if(zval.ge.ccut(nccut))   ipart = 2
      else
            if(cdfval.le.ccdf(1))     ipart = 0
            if(cdfval.ge.ccdf(nccut)) ipart = 2
      endif
c
c ARE WE IN THE LOWER TAIL?
c
      if(ipart.eq.0) then
            if(ltail.eq.1) then
c
c Straight Linear Interpolation:
c
                  powr = 1.0
                  if(zval.gt.UNEST) then
                        cdfval = powint(zmin,ccut(1),0.d0,ccdf(1),
     +                                  zval,powr)
                  else
                        zval = powint(0.d0,ccdf(1),zmin,ccut(1),
     +                                cdfval,powr)
                  endif
            else if(ltail.eq.2) then
c
c Power Model interpolation to lower limit "zmin"?
c
                  if(zval.gt.UNEST) then
                        cdfval = powint(zmin,ccut(1),0.d0,ccdf(1),
     +                                  zval,ltpar)
                  else
                        powr = 1.0 / ltpar
                        zval = powint(0.d0,ccdf(1),zmin,ccut(1),
     +                                cdfval,powr)
                  endif
c
c Linear interpolation between the rescaled global cdf?
c
            else if(ltail.eq.3) then
                  if(zval.gt.UNEST) then
c
c Computing the cdf value. Locate the point and the class bound:
c
                        call locate(cut,ncut,1,ncut,zval,idat)
                        call locate(cut,ncut,1,ncut,ccut(1),iupp)
c
c Straight linear interpolation if no data; otherwise, linear:
c
                        if(idat.le.0.or.idat.ge.ncut.or.
     +                     iupp.le.0.or.iupp.ge.ncut) then
                               cdfval = powint(zmin,cut(1),0.d0,cdf(1),
     +                                         zval,1.)
                         else
                               temp   = powint(cut(idat),cut(idat+1),
     +                                    cdf(idat),cdf(idat+1),zval,1.)
                               cdfval = temp*ccdf(1)/cdf(iupp)
                         endif
                   else
c
c Computing Z value: Are there any data out in the tail?
c
                         call locate(cut,ncut,1,ncut,ccut(1),iupp)
c
c Straight linear interpolation if no data; otherwise, local linear
c interpolation:
c
                        if(iupp.le.0.or.iupp.ge.ncut) then
                              zval = powint(0.d0,cdf(1),zmin,cut(1),
     +                                      cdfval,1.)
                        else
                              temp = cdfval*cdf(iupp)/ccdf(1)
                              call locate(cdf,ncut,1,ncut,temp,idat)
                              if(idat.le.0.or.idat.ge.ncut) then
                                    zval = powint(0.d0,cdf(1),zmin,
     +                                            cut(1),cdfval,1.)
                              else
                                    zval = powint(cdf(idat),cdf(idat+1),
     +                                   cut(idat),cut(idat+1),temp,1.)
                              end if
                        endif
                  endif
            else
c
c Error situation - unacceptable option:
c
                  ierr = 2
                  return
            endif
      endif
c
c FINISHED THE LOWER TAIL,  ARE WE IN THE MIDDLE?
c
      if(ipart.eq.1) then
c
c Establish the lower and upper limits:
c
            if(zval.gt.UNEST) then
                  call locate(ccut,nccut,1,nccut,zval,cclow)
            else
                  call locate(ccdf,nccut,1,nccut,cdfval,cclow)
            endif
            cchigh = cclow + 1
            if(middle.eq.1) then
c
c Straight Linear Interpolation:
c
                  powr = 1.0
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(cclow),ccut(cchigh),
     +                           ccdf(cclow),ccdf(cchigh),zval,powr)
                  else
                        zval = powint(ccdf(cclow),ccdf(cchigh),
     +                         ccut(cclow),ccut(cchigh),cdfval,powr)
                  endif
c
c Power interpolation between class bounds?
c
            else if(middle.eq.2) then
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(cclow),ccut(cchigh),
     +                           ccdf(cclow),ccdf(cchigh),zval,mpar)
                  else
                        powr = 1.0 / mpar
                        zval = powint(ccdf(cclow),ccdf(cchigh),
     +                         ccut(cclow),ccut(cchigh),cdfval,powr)
                  endif
c
c Linear interpolation between the rescaled global cdf?
c
            else if(middle.eq.3) then
                  call locate(cut,ncut,1,ncut,ccut(cclow),ilow)
                  call locate(cut,ncut,1,ncut,ccut(cchigh),iupp)
                  if(cut(ilow).lt.ccut(cclow))  ilow = ilow + 1
                  if(cut(iupp).gt.ccut(cchigh)) iupp = iupp - 1
                  if(zval.gt.UNEST) then
                        call locate(cut,ncut,1,ncut,zval,idat)
c
c Straight linear interpolation if no data; otherwise, local linear
c interpolation:
c
                        if(idat.le.0.or.idat.ge.ncut.or.
     +                     ilow.le.0.or.ilow.ge.ncut.or.
     +                     iupp.le.0.or.iupp.ge.ncut.or.
     +                     iupp.le.ilow) then
                              cdfval=powint(ccut(cclow),ccut(cchigh),
     +                                ccdf(cclow),ccdf(cchigh),zval,1.)
                        else
                              temp = powint(cut(idat),cut(idat+1),
     +                                cdf(idat),cdf(idat+1),zval,1.)
                              cdfval=powint(cdf(ilow),cdf(iupp),
     +                                ccdf(cclow),ccdf(cchigh),temp,1.)
                        endif
                  else
c
c Straight linear interpolation if no data; otherwise, local linear
c interpolation:
c
                        if(ilow.le.0.or.ilow.ge.ncut.or.
     +                     iupp.le.0.or.iupp.ge.ncut.or.
     +                     iupp.le.ilow) then
                              zval=powint(ccdf(cclow),ccdf(cchigh),
     +                              ccut(cclow),ccut(cchigh),cdfval,1.)
                        else
                              temp=powint(ccdf(cclow),ccdf(cchigh),
     +                              cdf(ilow),cdf(iupp),cdfval,1.)
                              call locate(cdf,ncut,1,ncut,temp,idat)
                              if(cut(idat).lt.ccut(cclow)) idat=idat+1
                              if(idat.le.0.or.idat.ge.ncut.or.
     +                           cut(idat+1).gt.ccut(cchigh)) then
                                    zval = powint(ccdf(cclow),
     +                                     ccdf(cchigh),ccut(cclow),
     +                                     ccut(cchigh),cdfval,1.)
                              else
                                    zval = powint(cdf(idat),cdf(idat+1),
     +                                   cut(idat),cut(idat+1),temp,1.)
                              end if
                              zval = powint(cdf(idat),cdf(idat+1),
     +                                cut(idat),cut(idat+1),temp,1.)
                        endif
                  endif
            else
c
c Error situation - unacceptable option:
c
                  ierr = 2
                  return
            endif
      endif
c
c FINISHED THE MIDDLE,  ARE WE IN THE UPPER TAIL?
c
      if(ipart.eq.2) then
            if(utail.eq.1) then
                  powr = 1.0
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(nccut),zmax,ccdf(nccut),
     +                                  1.0,zval,powr)
                  else
                        zval   = powint(ccdf(nccut),1.0,ccut(nccut),
     +                                  zmax,cdfval,powr)
                  endif

            else if(utail.eq.2) then
c
c Power interpolation to upper limit "utpar"?
c
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(nccut),zmax,ccdf(nccut),
     +                                  1.0,zval,utpar)
                  else
                        powr = 1.0 / utpar
                        zval   = powint(ccdf(nccut),1.0,ccut(nccut),
     +                                  zmax,cdfval,powr)
                  endif
c
c Linear interpolation between the rescaled global cdf?
c
            else if(utail.eq.3) then
                  if(zval.gt.UNEST) then
c
c Approximately Locate the point and the class bound:
c
                        call locate(cut,ncut,1,ncut,zval,idat)
                        call locate(cut,ncut,1,ncut,ccut(nccut),ilow)
                        if(cut(idat).lt.zval)        idat = idat + 1
                        if(cut(ilow).lt.ccut(nccut)) ilow = ilow + 1
c
c Straight linear interpolation if no data; otherwise, local linear
c interpolation:
c
                        if(idat.le.0.or.idat.ge.ncut.or.
     +                     ilow.le.0.or.ilow.ge.ncut) then
                              cdfval = powint(ccut(nccut),zmax,
     +                                  ccdf(nccut),1.0,zval,1.)
                        else
                              temp   = powint(cut(idat),cut(idat+1),
     +                                  cdf(idat),cdf(idat+1),zval,1.)
                              cdfval = powint(cdf(ilow),1.0,
     +                                  ccdf(nccut),1.0,temp,1.)
                        endif
                  else
c
c Computing Z value: Are there any data out in the tail?
c
                        call locate(cut,ncut,1,ncut,ccut(nccut),ilow)
                        if(cut(ilow).lt.ccut(nccut)) ilow = ilow + 1
c
c Straight linear interpolation if no data; otherwise, local linear
c interpolation:
c
                        if(ilow.le.0.or.ilow.ge.ncut) then
                              zval   = powint(ccdf(nccut),1.0,
     +                                  ccut(nccut),zmax,cdfval,1.)
                        else
                              temp = powint(ccdf(nccut),1.0,
     +                                cdf(ilow),1.0,cdfval,1.)
                              call locate(cdf,ncut,1,ncut,temp,idat)
                              if(cut(idat).lt.ccut(nccut)) idat=idat+1
                              if(idat.ge.ncut) then
                                    zval   = powint(ccdf(nccut),1.0,
     +                                       ccut(nccut),zmax,cdfval,1.)
                              else
                                    zval = powint(cdf(idat),cdf(idat+1),
     +                                   cut(idat),cut(idat+1),temp,1.)
                              endif
                        endif
                  endif
c
c Fit a Hyperbolic Distribution?
c
            else if(utail.eq.4) then
c
c Figure out "lambda" and required info:
c
                  lambda = (ccut(nccut)**utpar)*(1.0-ccdf(nccut))
                  if(zval.gt.UNEST) then
                        cdfval = 1.0 - (lambda/(zval**utpar))
                  else
                        zval = (lambda/(1.0-cdfval))**(1.0/utpar)
                  endif
            else
c
c Error situation - unacceptable option:
c
                  ierr = 2
                  return
            endif
      endif
      if(zval.lt.zmin) zval = zmin
      if(zval.gt.zmax) zval = zmax
c
c All finished - return:
c
      return
      end
