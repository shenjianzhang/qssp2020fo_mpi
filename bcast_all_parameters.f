      subroutine bcast_all_parameters()

      use my_mpi
      implicit none
c
      include 'qpglobal.h'
      integer*4 myrank,ierr
      integer*4 nparam_i
      parameter(nparam_i=7)
      integer*4 nparam_r
      parameter(nparam_r=15)
      integer*4 nparam_l
      parameter(nparam_l=4)
      integer*4 bcast_integer(nparam_i)
      real*8 bcast_real(nparam_r)
      logical*2 bcast_logical(nparam_l)
c
c     initialize
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
      bcast_integer(:) = 0
      bcast_real(:) = 0.d0
      bcast_logical(:) = .false.
c
c     master process
      if(myrank .eq. 0) then
          bcast_integer = (/ngrn, nf, ldeg1,ldeg2,dldeg,
     & ldeggr, l0/)
          bcast_real = (/grndep, dpr, flw1,fup1,flw2,fup2,
     & fcut, fi, qmmax,pkratio, fgr, rplanet, depatmos,rratmos,qsmin/)
            
          bcast_logical = (/nogravity, selpsv, selsh,dispersion/)
      endif

      call bcast_all_i(bcast_integer, nparam_i)
      call bcast_all_r(bcast_real, nparam_r)
      call bcast_all_l(bcast_logical, nparam_l)
c
      call bcast_all_r(dp0, 2*lymax)
      call bcast_all_r(vp0, 2*lymax)
      call bcast_all_r(vs0, 2*lymax)
      call bcast_all_r(ro0, 2*lymax)
      call bcast_all_r(qp0, 2*lymax)
      call bcast_all_r(qs0, 2*lymax)
c
      call bcast_all_r(dp0up, 2*lymax)
      call bcast_all_r(vp0up, 2*lymax)
      call bcast_all_r(vs0up, 2*lymax)
      call bcast_all_r(ro0up, 2*lymax)
      call bcast_all_r(qp0up, 2*lymax)
      call bcast_all_r(qs0up, 2*lymax)
      call bcast_all_r(dp0lw, 2*lymax)
      call bcast_all_r(vp0lw, 2*lymax)
      call bcast_all_r(vs0lw, 2*lymax)
      call bcast_all_r(ro0lw, 2*lymax)
      call bcast_all_r(qp0lw, 2*lymax)
      call bcast_all_r(qs0lw, 2*lymax)
c 
      if(myrank .ne. 0) then
          ngrn = bcast_integer(1)
          nf = bcast_integer(2)
          ldeg1 = bcast_integer(3)
          ldeg2 = bcast_integer(4)
          dldeg = bcast_integer(5)
          ldeggr = bcast_integer(6)
          l0 = bcast_integer(7)
c
          grndep = bcast_real(1)
          dpr = bcast_real(2)
          flw1 = bcast_real(3)
          fup1 = bcast_real(4)
          flw2 = bcast_real(5)
          fup2 = bcast_real(6)
          fcut = bcast_real(7)
          fi = bcast_real(8)
          qmmax = bcast_real(9)
          pkratio = bcast_real(10)
          fgr = bcast_real(11)
          rplanet = bcast_real(12)
          depatmos = bcast_real(13)
          rratmos = bcast_real(14)
          qsmin = bcast_real(15)
c
          nogravity = bcast_logical(1)
          selpsv = bcast_logical(2)
          selsh = bcast_logical(3)
          dispersion = bcast_logical(4)
      endif
      end subroutine bcast_all_parameters