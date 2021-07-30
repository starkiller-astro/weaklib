C23456789012345678901234567890123456789012345678901234567890123456789012
C
c      Program to compute spline fits to fermi integrals
cc  Must provide data file 94
      subroutine initferm(FilePath)

!     use parallel_module
!     use mpi

      implicit  double precision (a-h,o-z)
      parameter (n=201)
      character*128 FilePath
      dimension f32(n),f12(n),fm12(n),eta(n),fr(n)
      dimension f32a(n),f12a(n),fra(n),fia(n)
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia
      dimension send_buf(4,n)
C
!     if (myid.eq.0) then
      open(94,file=TRIM(FilePath)//'/fermi.atb',
     & status='old')
C
      do 10 i=1,n
       read(94,*)eta(i),f32(i),f12(i),fm12(i)
       send_buf(1,i) = eta(i)
       send_buf(2,i) = f32(i)
       send_buf(3,i) = f12(i)
       send_buf(4,i) = fm12(i)
 10   continue
C
      close(94,status='keep')
      
!     endif

!     i_extent=n*4
!     CALL MPI_BCAST( send_buf , i_extent, MPI_DOUBLE_PRECISION, 0,
!    &                     MPI_COMM_WORLD, ierr)

      do 20 i=1,n
       eta(i) = send_buf(1,i)
       f32(i) = send_buf(2,i)
       f12(i) = send_buf(3,i)
       fm12(i) = send_buf(4,i)
 20    fr(i)=f12(i)/fm12(i)

C
       call spline(eta,f12,n,f12a)
       call spline(eta,f32,n,f32a)
       call spline(eta,fr,n,fra)
       call spline(f12,eta,n,fia)
       return
       end subroutine initferm
