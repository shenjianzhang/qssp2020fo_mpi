      subroutine qpshkern(f,ldeg,ysh)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ysh(2,2): sh solution vector (complex)
c
      integer*4 ldeg
      real*8 f
      complex*16 ysh(2,2)
c
      include 'qpglobal.h'
c
      integer*4 i,j,ly,istp,lylw,lwup
      real*8 xmin
      complex*16 cldeg,cll1,cruplw,delta
      complex*16 spbphj,spbphy
c
      complex*16 c1,c2
      data c1,c2/(1.d0,0.d0),(2.d0,0.d0)/
c
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      lylw=lylwsh(ldeg)
c
      if(ldeg.eq.0.or.lyr.lt.lyob.or.lylw.lt.max0(lys,lyr))return
c
      if(ldeg.eq.1.and.cdabs(comi).le.0.d0)then
        if(lyr.le.lys)then
          ysh(1,1)= cmomlw*(crrup(lyr)/crrup(lys))/cmomsum
        else
          ysh(1,1)=-cmomup*(crrup(lyr)/crrup(lys))/cmomsum
        endif
        return
      endif
c
      xmin=dsqrt(2.d0*dble(2*ldeg+1))
      cldeg=dcmplx(dble(ldeg),0.d0)
      cll1=cldeg*(cldeg+c1)
c
      do ly=lyob,min0(lylwsh(ldeg),lycm-1)
        if(ldeg.le.ldegsh(ly))then
          if(cdabs(kt(ly))*rrup(ly).le.xmin)then
            ksmallsh(ly)=.true.
            if(rrlw(ly).gt.0.d0)then
              cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
            else
              cruplw=(0.d0,0.d0)
            endif
            wj(ldeg,ly,3)=cldeg*cruplw
            wh(ldeg,ly,3)=-(cldeg+c1)*cruplw
          else if(ksmallsh(ly))then
            ksmallsh(ly)=.false.
            call spbjh(ldegsh(ly),kt(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,3),zjlw(0,ly,3),wj(0,ly,3),
     &                 zhup(0,ly,3),zhlw(0,ly,3),wh(0,ly,3))
          endif
        endif
      enddo
c
      do ly=lyob,min0(lylwsh(ldeg),lycm-1)
        if(ly.lt.lys)then
          lwup=0
          cpt(1,ly)=-wj(ldeg,ly,3)
          cpt(2,ly)=-wh(ldeg,ly,3)
        else
          lwup=1
          cpt(1,ly)=wj(ldeg,ly,3)
          cpt(2,ly)=wh(ldeg,ly,3)
        endif
        call qptmat(ldeg,ly,lylw,lwup)
      enddo
      call qptprop(ysh,lylw)
      return
      end

