#include "fintrf.h"

c ------------------------------------------------------------------------
c     hwm93.f - MATLAB gateway to hwm93
c     version 1.0 (06 August 2010)
c     gateway written by Tim Duly <duly2@illinois.edu>
c ------------------------------------------------------------------------
c     this is for a 64 bit Linux machine (xbuntu 10.04)
c ------------------------------------------------------------------------
      subroutine mexfunction(nlhs,plhs,nrhs,prhs)
c      implicit real(4) (a-h,q-z) ------ get rid of this, not needed!
      real(4) w(2),ap(2)
      real*8 w_8bit(2),ap_8bit(2)

c Required declarations for MATLAB
	mwPointer plhs(*), prhs(*)
	integer nlhs,nrhs
c Prototype the return types of these functions
        mwPointer mxGetPr, mxCreateDoubleMatrix
        real*8 mxGetScalar
	
	mwsize size1,size2

c Get variables passed from MATLAB, properly typed
        iyd = int(mxgetscalar(prhs(1)))
        sec = sngl(mxgetscalar(prhs(2)))
        alt = sngl(mxgetscalar(prhs(3)))
        glat = sngl(mxgetscalar(prhs(4)))
        glong = sngl(mxgetscalar(prhs(5)))
        stl = mxgetscalar(prhs(6))
	f107a = mxgetscalar(prhs(7))
	f107 = mxgetscalar(prhs(8))
	ap(1) = mxgetscalar(prhs(9))
	ap(2) = mxgetscalar(prhs(10))

c ap(1) is daily, ap(2) is 3hr

c Do the deed: call hwm93
        call gws5(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,w)

c Cast into 8-bit
	w_8bit = w;
	print *, w_8bit	
c Create MATLAB return matrices
	size1 = 1
	size2 = 2
        plhs(1) = mxCreateDoubleMatrix(size1,size2,0)

c Put the returned HWM93 data into the MATLAB return matrices we created        
        call mxCopyReal8ToPtr(w_8bit,mxGetPr(plhs(1)),size2)

c Done
        return
        end

        
        
