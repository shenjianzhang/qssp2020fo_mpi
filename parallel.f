      module my_mpi
      use mpi
      implicit none
      integer*4 my_local_mpi_comm_world
c
      end module
c
      subroutine bcast_all_i(buffer, countval)
      use my_mpi
      implicit none
      integer*4 countval
      integer*4 buffer(countval)
      integer*4 ierr
c
      call MPI_BCAST(buffer, countval, MPI_INTEGER, 0,
     &     my_local_mpi_comm_world, ierr)
c
      end subroutine bcast_all_i

      subroutine bcast_all_r(buffer, countval)
      use my_mpi
      implicit none
      integer*4 countval
      real*8 buffer(countval)
      integer*4 ierr
c
      call MPI_BCAST(buffer, countval, MPI_DOUBLE_PRECISION, 0,
     &     my_local_mpi_comm_world, ierr)
c
      end subroutine bcast_all_r

      subroutine bcast_all_l(buffer, countval)
      use my_mpi
      implicit none
      integer*4 countval
      logical*2 buffer(countval)
      integer*4 ierr
c
      call MPI_BCAST(buffer, countval, MPI_LOGICAL, 0,
     &     my_local_mpi_comm_world, ierr)
c
      end subroutine bcast_all_l

      subroutine synchronize_all()
      use my_mpi

      integer*4 ierr
      call MPI_BARRIER(my_local_mpi_comm_world, ierr)
      if(ierr .ne. 0) then
          stop 'Error synchronize'
      endif
      end subroutine synchronize_all