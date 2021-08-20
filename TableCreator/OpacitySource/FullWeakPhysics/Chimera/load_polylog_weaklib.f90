SUBROUTINE load_polylog_weaklib

   USE polylog_module_weaklib

   implicit none

   integer :: i !loop counter

   CHARACTER (len=128)                  :: c_diagnostic  ! description of problem terminating run

 1001 FORMAT (' nnu=',i3,' exceeds nnud=',i3,' in subroutine sctekrnl')

   xa(0)=0.d0
   pl2a(0)=0.d0
   pl3a(0)=0.d0

   open(unit=20,file='../Chimera/polylog.tab',status='old',position="rewind")

   do i=1,nt
     read(20,*) xa(i),pl2a(i),pl3a(i)
   enddo

   close(20)

end subroutine load_polylog_weaklib
