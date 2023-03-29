      subroutine qpdifmatl(ly,ldeg,rr,mat)
      implicit none
      integer*4 ly,ldeg
      real*8 rr
      complex*16 mat(6,6)
c
      include 'qpglobal.h'
c
c     4x4 coefficient matrix for spheroidal mode l > 0 in liquid media
c
      real*8 f0,f,mass,rorr,beta,dro,ro1,grrr,vprr,n2rr
      complex*16 cldeg,clp1,cll1
      complex*16 cup,clw,crr,crorr,cgrrr,cvprr
      complex*16 cgarr,gamma,cn2rr,cqpfac,cqnfac
c
      complex*16 c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      cll1=cldeg*clp1
c
      crr=dcmplx(rr,0.d0)
      cup=(crr-crrlw(ly))/(crrup(ly)-crrlw(ly))
      clw=c1-cup
      cvprr=cup*cvpup(ly)+clw*cvplw(ly)
c
      f=dreal(comi)/PI2
      vprr=cdabs(cvprr)
      cqpfac=cvprr/dcmplx(vprr,0.d0)
      beta=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
c
      if(rr.le.rrlw(ly))then
        mass=0.d0
        crorr=crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
      else if(ly.ge.lyos)then
        crorr=cup*croup(ly)+clw*crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
        dro=(rorr-rolw(ly))/(rr-rrlw(ly))
        ro1=rorr-dro*rrlw(ly) 
        mass=PI*(rr-rrlw(ly))*((4.d0/3.d0)*ro1
     &        *(rr**2+rr*rrlw(ly)+rrlw(ly)**2)
     &        +dro*(rr**3+rr**2*rrlw(ly)
     &        +rr*rrlw(ly)**2+rrlw(ly)**3))
      else
        rorr=rolw(ly)*dexp(dlog(roup(ly)/rolw(ly))
     &                  *(rr-rrlw(ly))/(rrup(ly)-rrlw(ly)))
        crorr=dcmplx(rorr,0.d0)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
        mass=2.d0*PI2/beta**3
     &      *(rorr*(beta*rr*(beta*rr-2.d0)+2.d0)
     &       -rolw(ly)*(beta*rrlw(ly)*(beta*rrlw(ly)-2.d0)+2.d0))
      endif
      cgrrr=cgrlw(ly)*(crrlw(ly)/crr)**2
     &     +dcmplx(BIGG*mass/rr**2,0.d0)
      grrr=cdabs(cgrrr)
      n2rr=dmax1(0.d0,-grrr*(beta+grrr/vprr**2))
      cn2rr=-cgrrr*(dcmplx(beta,0.d0)+cgrrr/cvprr**2)
      cqnfac=cn2rr/dcmplx(cdabs(cn2rr),0.d0)
c
c     stabilization at the static limit
c
      if(ly.lt.lyos)then
        f0=dmax1(fbvatm,0.5d0*dsqrt(n2rr)/PI2)
        if(fbvatm.le.0.d0)then
          n2rr=0.d0
        else if(f.le.f0)then
c          n2rr=n2rr*(f/f0)**2.5d0
c          vprr=dsqrt(-grrr/(n2rr/grrr+beta))
        endif
      else if(ly.ge.lyos.and.ly.lt.lyob)then
        f0=dmax1(fbvocean,0.5d0*dsqrt(n2rr)/PI2)
        if(fbvocean.le.0.d0)then
          n2rr=0.d0
        else if(f.le.f0)then
          n2rr=n2rr*(f/f0)**2.5d0
          vprr=dsqrt(-grrr/(n2rr/grrr+beta))
        endif
      else if(ly.ge.lycm.and.ly.lt.lycc)then
        f0=dmax1(fbvcore,0.5d0*dsqrt(n2rr)/PI2)
        if(fbvcore.le.0.d0)then
          n2rr=0.d0
        else if(f.le.f0)then
          n2rr=n2rr*(f/f0)**2.5d0
          vprr=dsqrt(-grrr/(n2rr/grrr+beta))
        endif
      else
        stop 'Error in qpdifmatl: unknown fiquid layer!'
      endif
c
c     Brunt-Väisälä frequency
c
      cn2rr=dcmplx(n2rr,0.d0)*cqnfac
      cvprr=dcmplx(vprr,0.d0)*cqpfac
c
      mat(1,1)=cgrrr/cvprr**2-c1/crr
      mat(1,2)=cll1/crr-comi2*crr/cvprr**2
      mat(1,3)=-crr/cvprr**2
      mat(1,4)=(0.d0,0.d0)
c
      if(n2rr.le.0.d0)then
        mat(2,1)=c1/crr
        mat(2,2)=(0.d0,0.d0)
        mat(2,3)=(0.d0,0.d0)
        mat(2,4)=(0.d0,0.d0)
      else
        mat(2,1)=(c1-cn2rr/comi2)/crr
        mat(2,2)=cn2rr/cgrrr
        mat(2,3)=cn2rr/(comi2*cgrrr)
        mat(2,4)=(0.d0,0.d0)
      endif
c
      mat(3,1)=cgarr/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=-clp1/crr
      mat(3,4)=c1/crr
c
      mat(4,1)=cgarr*clp1/crr
      mat(4,2)=-cgarr*cll1/crr
      mat(4,3)=(0.d0,0.d0)
      mat(4,4)=cldeg/crr
      return
      end