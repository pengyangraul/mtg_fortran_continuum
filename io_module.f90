module io_module
!  USE netcdf
  USE mat_module
  implicit none
  contains
!    subroutine nfcheck(status)
!      implicit none
!      integer, intent (in) :: status
!      if(status /= nf90_noerr) then
!        print *, trim(nf90_strerror(status))
!        stop 2
!      end if
!    end subroutine nfcheck

  function read_param_string(nth)
  implicit none
      integer, intent(in) :: nth
      character(len=40) :: read_param_string
      call GETARG(nth,read_param_string)
  end function read_param_string

    function read_param_real(nth)
      implicit none
      real(kind = np) :: read_param_real
      integer, intent(in) :: nth
      character(len=40) :: param_str
      call GETARG(nth,param_str)
      read(param_str,*)  read_param_real
    end function

    function read_param_int(nth)
      implicit none
      integer :: read_param_int
      integer, intent(in) :: nth
      character(len=40) :: param_str
      call GETARG(nth,param_str)
      read(param_str,*)  read_param_int
    end function


    function read_param_complex(nth)
      implicit none
      complex(kind = np) :: read_param_complex
      integer, intent(in) :: nth
      character(len=40) :: param_str
      call GETARG(nth,param_str)
      read(param_str,*)  read_param_complex
    end function
  end module
