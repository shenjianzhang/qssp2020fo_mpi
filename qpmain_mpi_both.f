      program qssp
      use mpi
      use my_mpi
      implicit none
c
      include 'qpglobal.h'
      integer*4 myrank
      integer*4 ig,ierr,runtime
      integer*4 time
      character*80 inputfile, sphefile
      character*10 ldeg1_str,ldeg2_str
      integer*4 i,j,ldeg,nm,ldeg_w
      real*8 f, freq_w, amp_w, phase_w, Q_w
      integer ldeg_start, ldeg_end, ldeg_step, ldeg_per_proc
      integer ldeg_rem, ldeg_local_start, ldeg_local_end, ldeg_total
      integer*4 ldegl,ldegr
      integer*4 nl_local, indx_local
      integer*4,allocatable :: l_local(:),l_master(:),nl_master(:)
      integer*4,allocatable :: counts(:), displs(:)
      integer*4,allocatable :: nm_loc_psv(:),nm_mas_psv(:)
      integer*4,allocatable :: nm_loc_sh(:),nm_mas_sh(:)
      real*8,allocatable :: freq_loc_psv(:,:),amp_loc_psv(:,:)
      real*8,allocatable :: phase_loc_psv(:,:),Q_loc_psv(:,:)
      real*8,allocatable :: freq_loc_sh(:,:),amp_loc_sh(:,:)
      real*8,allocatable :: phase_loc_sh(:,:),Q_loc_sh(:,:)
      real*8,allocatable :: freq_mas_psv(:,:),amp_mas_psv(:,:)
      real*8,allocatable :: phase_mas_psv(:,:),Q_mas_psv(:,:)
      real*8,allocatable :: freq_mas_sh(:,:),amp_mas_sh(:,:)
      real*8,allocatable :: phase_mas_sh(:,:),Q_mas_sh(:,:)
c
      real*8 fm(nmmax),qm(nmmax)
      complex*16 cam(nmmax)
      complex*16 cpsvres(nfmax),cshres(nfmax)
      complex*16 cpsvspc(nfmax),cshspc(nfmax)
c
c     MPI variables
      integer*4 numprocs, ierr_mpi
      integer*4 tag, status(MPI_STATUS_SIZE)
c      inputfile='test.inp'
c     Initialize MPI
      call MPI_Init(ierr_mpi)
      
c     Get the rank and number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr_mpi)
      call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr_mpi)
      my_local_mpi_comm_world = MPI_COMM_WORLD
c
c      write(*,'(i10,i10)')numprocs,myrank
c     read input file file
c
      if (myrank == 0) then      
        write(*,'(a,$)')' the input data file is '
        read(*,'(a)')inputfile
        open(10,file=inputfile,status='old')
        call qpgetinp(10)
        close(10)
        runtime=time()
      endif
c     broadcast parameters read from master to all processes
      call synchronize_all()
      call bcast_all_parameters()
c
c     calculate at each process
      call qpsublayer(ierr)
c
      do i=1,nf
          cpsvspc(i)=(0.d0,0.d0)
          cpsvres(i)=(0.d0,0.d0)
          cshspc(i)=(0.d0,0.d0)
          cshres(i)=(0.d0,0.d0)
      enddo
c
      lys=lygrn(1)
      selsh=ldeg2.ge.1.and.lys.ge.lyob.and.lyr.ge.lyob
c     
c     Distribute work among processes
      ldeg_start = ldeg1
      ldeg_end = ldeg2
      ldeg_step = 1
      ldeg_per_proc = (ldeg_end - ldeg_start + 1) / numprocs
      ldeg_rem = (ldeg_end - ldeg_start + 1) - numprocs* ldeg_per_proc
      if (myrank .lt. ldeg_rem) then
        ldeg_local_start = myrank * (ldeg_per_proc + 1) + ldeg_start
      else
        ldeg_local_start=myrank*ldeg_per_proc+ldeg_start+ldeg_rem
      endif     
      if (myrank+1 .lt. ldeg_rem) then
        ldeg_local_end=(myrank+1)*ldeg_per_proc+ldeg_start+myrank
      else
        ldeg_local_end=(myrank+1)*ldeg_per_proc+ldeg_start+ldeg_rem-1
      endif
      if (myrank == numprocs - 1) ldeg_local_end = ldeg_end
c
c     prepare arrays to store l,n,f,amplitude, phase, Q
      nl_local = ldeg_local_end - ldeg_local_start +1 
      allocate(l_local(nl_local))
      allocate(nm_loc_psv(nl_local),nm_loc_sh(nl_local))
      allocate(freq_loc_psv(nmmax,nl_local), amp_loc_psv(nmmax,nl_local))
      allocate(phase_loc_psv(nmmax,nl_local), Q_loc_psv(nmmax,nl_local))
      allocate(freq_loc_sh(nmmax,nl_local), amp_loc_sh(nmmax,nl_local))
      allocate(phase_loc_sh(nmmax,nl_local), Q_loc_sh(nmmax,nl_local))
      l_local = -1
      nm_loc_psv = 0
      nm_loc_sh = 0
      freq_loc_psv = 0.d0
      amp_loc_psv = 0.d0
      phase_loc_psv = 0.d0
      Q_loc_psv = 0.d0
      freq_loc_sh = 0.d0
      amp_loc_sh = 0.d0
      phase_loc_sh = 0.d0
      Q_loc_sh = 0.d0
c
c     most important part
      do ldeg=ldeg_local_start,ldeg_local_end,ldeg_step
        indx_local = ldeg - ldeg_local_start +1
        if(ldeg2 .eq.ldeg1)then
          flw=flw1
          fup=fup1
        else
          flw=flw1+(flw2-flw1)*dble(ldeg-ldeg1)/dble(ldeg2-ldeg1)
          fup=fup1+(fup2-fup1)*dble(ldeg-ldeg1)/dble(ldeg2-ldeg1)
        endif
        df=(fup-flw)/dble(nf-1)
c       calculate spectrum of both P-SV and SH waves
        call qpgrnspec(1, ldeg)
c       find peaks in the spectrum of ldeg (P-SV)
        call fmodes(ldeg,flw,fup,fi,nf,psvspecr,psvspeci,
     &              cam,nm,nmmax,qmmax,pkratio,freq_loc_psv(:,indx_local), 
     &             amp_loc_psv(:,indx_local), phase_loc_psv(:,indx_local),
     &              Q_loc_psv(:,indx_local))
        l_local(indx_local) = ldeg
        nm_loc_psv(indx_local) = nm
c
c       calculate spectrum of SH waves
        if(ldeg.lt.1.or..not.selsh)goto 50
c       find peaks in the spectrum of ldeg (SH)
        call fmodes(ldeg,flw,fup,fi,nf,shspecr,shspeci,
     &              cam,nm,nmmax,qmmax,pkratio,freq_loc_sh(:,indx_local), 
     &              amp_loc_sh(:,indx_local), phase_loc_sh(:,indx_local),
     &              Q_loc_sh(:,indx_local))
        nm_loc_sh(indx_local) = nm

50      continue
      enddo
c
c     gather result arrays from processes to root
      allocate(nl_master(numprocs))
      nl_master(:) = 0
c     define N-array for the columns of each process
      call MPI_ALLGATHER(nl_local, 1, MPI_INTEGER, nl_master, 1,
     & MPI_INTEGER,  MPI_COMM_WORLD, ierr_mpi)
c     define counts and displs
      allocate(counts(numprocs), displs(numprocs))
      do i = 1, numprocs
        counts(i) = nl_master(i) * nmmax
        displs(i) = sum(counts(1:i-1))
      enddo
c     NOTE: nl_master, counts, displs MUST be defined at all processes.
c     define arrays in master-process to store results
      if(myrank .eq. 0) then
        ldeg_total = ldeg_end - ldeg_start + 1
        allocate( l_master(ldeg_total) )
        allocate( nm_mas_psv(ldeg_total),nm_mas_sh(ldeg_total) )
        allocate( freq_mas_psv(nmmax, ldeg_total) )
        allocate( amp_mas_psv(nmmax, ldeg_total) )
        allocate( phase_mas_psv(nmmax, ldeg_total) )
        allocate( Q_mas_psv(nmmax, ldeg_total) )
        allocate( freq_mas_sh(nmmax, ldeg_total) )
        allocate( amp_mas_sh(nmmax, ldeg_total) )
        allocate( phase_mas_sh(nmmax, ldeg_total) )
        allocate( Q_mas_sh(nmmax, ldeg_total) )
      endif
c
c     gather
c     l-degree
      call MPI_Gatherv(l_local,counts(myrank+1)/nmmax,MPI_INTEGER,
     & l_master,counts/nmmax,displs/nmmax,MPI_INTEGER, 
     & 0, MPI_COMM_WORLD, ierr_mpi)
c     modes found for P-SV waves
      call MPI_Gatherv(nm_loc_psv,counts(myrank+1)/nmmax,MPI_INTEGER,
     & nm_mas_psv,counts/nmmax,displs/nmmax,MPI_INTEGER, 
     & 0, MPI_COMM_WORLD, ierr_mpi)
c     freq., amplitude, phase, Q-value for P-SV waves
      call MPI_Gatherv(freq_loc_psv,counts(myrank+1),
     & MPI_DOUBLE_PRECISION,freq_mas_psv,counts,displs,
     & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Gatherv(amp_loc_psv,counts(myrank+1),
     & MPI_DOUBLE_PRECISION,amp_mas_psv,counts,displs,
     & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Gatherv(phase_loc_psv,counts(myrank+1),
     & MPI_DOUBLE_PRECISION,phase_mas_psv,counts,displs,
     & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Gatherv(Q_loc_psv,counts(myrank+1),MPI_DOUBLE_PRECISION,
     & Q_mas_psv,counts,displs,MPI_DOUBLE_PRECISION, 
     & 0, MPI_COMM_WORLD, ierr_mpi)
c     modes found for SH waves
      call MPI_Gatherv(nm_loc_sh,counts(myrank+1)/nmmax,MPI_INTEGER,
     & nm_mas_sh,counts/nmmax,displs/nmmax,MPI_INTEGER, 
     & 0, MPI_COMM_WORLD, ierr_mpi)
c     freq., amplitude, phase, Q-value for SH waves
      call MPI_Gatherv(freq_loc_sh,counts(myrank+1),
     & MPI_DOUBLE_PRECISION,freq_mas_sh,counts,displs,
     & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Gatherv(amp_loc_sh,counts(myrank+1),
     & MPI_DOUBLE_PRECISION,amp_mas_sh,counts,displs,
     & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Gatherv(phase_loc_sh,counts(myrank+1),
     & MPI_DOUBLE_PRECISION,phase_mas_sh,counts,displs,
     & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Gatherv(Q_loc_sh,counts(myrank+1),MPI_DOUBLE_PRECISION,
     & Q_mas_sh,counts,displs,MPI_DOUBLE_PRECISION, 
     & 0, MPI_COMM_WORLD, ierr_mpi)
c
c     write the final results
      if(myrank.eq.0)then
c     prepare file to write
        write(ldeg1_str, '(I5)') ldeg1
        write(ldeg2_str, '(I5)') ldeg2
        sphefile=trim(adjustl(smodesfile))//trim(adjustl(ldeg1_str))//
     & '.'//trim(adjustl(ldeg2_str))//'.dat'
        open(21,file=sphefile,status='unknown')
        write(21,'(a)')'   l   n   f[Hz]       Amplitude'
     &             //'  Phase[deg]             Q'
        if(ldeg2 .gt. 0 .and. selsh)
          torofile=trim(adjustl(smodesfile))//trim(adjustl(ldeg1_str))//
     &    '.'//trim(adjustl(ldeg2_str))//'.dat'
          open(22,file=torofile,status='unknown')
          write(22,'(a)')'   l   n   f[Hz]       Amplitude'
     &             //'  Phase[deg]             Q'
        endif
c     write to file
        do i=1,ldeg_total
c         result of P-SV waves
          ldeg_w = l_master(i)
          do j = 1,nm_mas_psv(i)
            freq_w = freq_mas_psv(j,i)
            amp_w = amp_mas_psv(j,i)
            phase_w = phase_mas_psv(j,i)
            Q_w = Q_mas_psv(j,i)
            write(21,1000)ldeg_w, j-1, freq_w, amp_w, phase_w, Q_w
          enddo
c         result of SH waves
          if(ldeg.lt.1.or..not.selsh)goto 50
          do j = 1,nm_mas_sh(i)
            freq_w = freq_mas_sh(j,i)
            amp_w = amp_mas_sh(j,i)
            phase_w = phase_mas_sh(j,i)
            Q_w = Q_mas_sh(j,i)
            write(22,1000)ldeg_w, j-1, freq_w, amp_w, phase_w, Q_w
          enddo
        enddo
        close(21)
        if(ldeg2 .gt. 0 .and. selsh)
          close(22)
        endif
c       release memory at master process
        deallocate( l_master,nm_master)
        deallocate( freq_master,amp_master )
        deallocate( phase_master,Q_master )
        runtime=time()-runtime
        write(*,'(a)')' #############################################'
        write(*,'(a)')' #                                           #'
        write(*,'(a)')' #      End of computations with qsspfo      #'
        write(*,'(a)')' #                                           #'
        write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
        write(*,'(a)')' #############################################'
      endif
c       release memory at all processes
      deallocate(nl_master,counts, displs)
      deallocate(l_local,nm_local)
      deallocate(freq_local, amp_local)
      deallocate(phase_local, Q_local)
c
1000  format(2i4,2E16.8,f12.4,E16.8)
c     Finalize MPI
      call MPI_Finalize(ierr_mpi)     
c
c
      stop
      end
