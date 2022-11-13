! PURPOSE:
!  - Write arrays to a binary files (array_to_binary)
!  - Read binary files into arrays (binary_to_array)

! Remember to include BINARY_IO before using the module :)

! NOTE: form='binary' is not portable and should **never be used** in any
!       new codes. The outdated 'binary' functionality is emulated in
!       F2003 standard combination of 'unformatted' and access='stream';
!       status='replace' makes sure the file is opened strictly for writing
!       and not for appending (default for stream is system dependent).

module binary_IO

  ! We get only the basic precision value RP from the module
  ! parameter_values
  ! This precision parameter is defined as (can also be set directly here):
  ! integer, parameter, public :: RP  = selected_real_kind(15,  307)
  use parameter_values, only : RP

  ! We get the standard error unit, which may be system-specific. It is
  ! necessary for error output. Error messages must go to the standard
  ! error device rather than standard output. This adheres to the standards
  ! and allows error redirection. In most cases, standard error goes to the
  ! terminal along with the standard output.
  use, intrinsic :: ISO_FORTRAN_ENV, only : ERROR_UNIT

  ! Getting path to result folder and renaming it to FIX_FOLDER
  ! NB! Spesific for the Hormone Model
  use folderpath, only : FIX_FOLDER => foldername


  implicit none

  ! Global constants.
  ! used only within the module
  integer, parameter, private :: MAX_UNIT=255      ! Maximum unit number
  integer, parameter, private :: FAKE_ERR = -9999  ! Fake error array value

  ! Error code returned in case of rank mismatch between the array in the file
  ! and main the output array passed to the procedure
  integer, parameter, private :: ERROR_RANK_MISMATCH = FAKE_ERR

  ! File form for saving binary files
  character(len=*), parameter, private :: file_form = 'unformatted'

  ! Fixed folder path, example: FIX_FOLDER = "./data/" (note slashes on both
  ! sides, otherwise FIX_FOLDER is just the prefix for the file name).
  ! Warning: do not forget to use trim(FIX_FOLDER) function because FIX_FOLDER
  !          is fixed length string!
  !character(255), public :: FIX_FOLDER = ""
  ! Note: FIX_FOLDER can also be declared allocatable, but then may need
  !       assigning explicit value. Without assigning an uninitialized variable
  !       can be dangerous to use
  !character(len=:), allocatable, public :: FIX_FOLDER

  ! Generic interface for printing  can use array_to_binary in the calling
  ! program.
  interface array_to_binary
    module procedure array_to_binary_real_1d
    module procedure array_to_binary_real_2d
    module procedure array_to_binary_real_3d
    module procedure array_to_binary_real_4d
    module procedure array_to_binary_real_5d
    module procedure array_to_binary_int_1d
    module procedure array_to_binary_int_2d
    module procedure array_to_binary_int_3d
    module procedure array_to_binary_int_4d
    module procedure array_to_binary_int_5d
  end interface array_to_binary

  ! Generic interface for reading, can use binary_to_array in the calling
  ! program.
  interface binary_to_array
    module procedure binary_to_array_real_1d
    module procedure binary_to_array_real_2d
    module procedure binary_to_array_real_3d
    module procedure binary_to_array_real_4d
    module procedure binary_to_array_real_5d
    module procedure binary_to_array_int_1d
    module procedure binary_to_array_int_2d
    module procedure binary_to_array_int_3d
    module procedure binary_to_array_int_4d
    module procedure binary_to_array_int_5d
  end interface binary_to_array

  contains

  !----------------------------------------------------------------------------------
  ! TOOLS FOR WRITING TO FILES
  !----------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! # GET_FREE_FUNIT #
  ! Purpose: returns the first free Fortran unit number (search from 1 to
  ! MAX_UNIT).
  ! ## RETURNS: ##
  !  - Integer unit number
  ! ## CALL PARAMETERS ##
  !  - optional logical execution error status (.TRUE.)
  !  - optional integer max_funit to search (default MAX_UNIT defined in
  !    module)
  !
  ! Author: John Burkardt : This code is distributed under the GNU LGPL license.
  ! Modified by Sergey Budaev.
  !
  ! This subroutine is a part of the *HEDTOOLS* suite. The stable SVN trunk
  ! address is here: https://tegsvn.uib.no/svn/tegsvn/tags/HEDTOOLS/1.1
  !-----------------------------------------------------------------------------
  function GET_FREE_FUNIT (file_status, max_funit) result (file_unit)

    use, intrinsic :: ISO_FORTRAN_ENV     !Provides system-wide scalar constants
                                          !INPUT_UNIT OUTPUT_UNIT ERROR_UNIT
    implicit none

    !Function value
    integer :: file_unit

    !Calling parameters
    logical, optional, intent(out) :: file_status
    integer, optional, intent(in)  :: max_funit

    !Local variables: copies of optional parameters we need co copy optional
    logical  :: file_status_here    !Variables in case they are absent, so always
    integer  :: max_funit_here      !Work with copies of optionals inside
    !Other local variables
    integer :: i
    integer :: ios
    logical :: lopen

    !Subroutine name for DEBUG LOGGER
    character (len=*), parameter :: PROCNAME = "GET_FREE_FUNIT"

    file_unit = 0
    file_status_here = .FALSE.

    if (present(max_funit)) then
        max_funit_here=max_funit
      else
        max_funit_here=MAX_UNIT !Max from globals
    end if

    do i=1, max_funit_here
      if (i /= INPUT_UNIT .and. i /= OUTPUT_UNIT .and. &
            i /= ERROR_UNIT) then             !Exclude standard console units
        inquire (unit=i, opened=lopen, iostat=ios)
        if (ios == 0) then
          if (.not. lopen) then
            file_unit = i                     !First free unit found
            file_status_here=.TRUE.
            if (present(file_status)) file_status=file_status_here
            return
          end if
        end if
      end if
    end do

    if (.not. file_status_here) file_unit=-1    !If no free unit found return -1
    if (present(file_status)) file_status=file_status_here ! and error flag

  end function GET_FREE_FUNIT

  !----------------------------------------------------------------------------------
  ! DIAGNOSTIC FUNCTIONS FOR BINARY FILES
  !----------------------------------------------------------------------------------

  ! # Diagnostic functions for binary files #
  ! ## Notes ##
  ! First, there are two diagnostic functions:
  ! - `binary_get_rank` determines what is the **rank** of the multidimensional
  !    array saved in the binary file. If the file cannot be opened for any
  !    reason, returns -1.
  ! - `binary_get_dims` returns the sizes of all the **dimensions** of the
  !    multidimensional array saved into the binary file. If the file cannot
  !    be opened for any reason, return an empty zero-size array.
  ! These functions should work with matrices of any dimensionality and any
  ! type (real and integer) because they don't read the actual data, only the
  ! header part of the binary file.
  ! These functions can be used as any other functions, e.g. in the
  ! `if`-statements: `if (binary_get_rank("data_file.bin")==4) call something`.

  !Determine what is the rank of the array from the binary file.
  function binary_get_rank(file_name, funit) result (array_rank_get)
    !Return value
    integer :: array_rank_get

    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    !Local variables
    integer :: funit_loc, error_status

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the file.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !Check file open error. If okay, read header.
    if (error_status == 0) then
      !Get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
    else
      !If there was an error opening the file return -1 as error code.
      array_rank_get = -1
    end if
    !Finally, close the file.
    close(funit_loc)

  end function binary_get_rank


  !Determine what are the dimensions of the data array saved into the
  ! binary file. Dimensions of the array are returned as an allocatable
  ! array of the size equal to the rank of the array.
  function binary_get_dims(file_name, funit) result (dim_length)
    !Return value is actually a 4d array
    integer, allocatable, dimension(:) :: dim_length
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    !Local variables
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !Check file open error. If okay, read header.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Knowing the rank, we can now allocate the output array to this size
      allocate(dim_length(array_rank_get))
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
    else
      !If there was an error opening the file return null-size output array.
      allocate(dim_length(0))
    end if
    !Finally, close the file.
    close(funit_loc)

  end function binary_get_dims

  !----------------------------------------------------------------------------------
  ! WRITE ARRAY TO BINARY FILE
  !----------------------------------------------------------------------------------

  !Write to binary for 1D array with real values
  subroutine array_to_binary_real_1d(array, file_name, funit, error)
    integer, parameter :: array_rank = 1             !For this version it is 1
    !Calling parameters of the subroutine
    real(kind=RP), dimension(:), intent(in) :: array !1D real array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_real_1d


  !Write to binary for 2D array with real values
  subroutine array_to_binary_real_2d(array, file_name, funit, error)
    integer, parameter :: array_rank = 2               !For this version it is 2
    !Calling parameters of the subroutine
    real(kind=RP), dimension(:,:), intent(in) :: array !2D real array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_real_2d


  !Write to binary for 3D array with real values
  subroutine array_to_binary_real_3d(array, file_name, funit, error)
    integer, parameter :: array_rank = 3                 !For this version it is 3
    !Calling parameters of the subroutine
    real(kind=RP), dimension(:,:,:), intent(in) :: array !3D real array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_real_3d


  !Write to binary for 4D array with real values
  subroutine array_to_binary_real_4d(array, file_name, funit, error)
    integer, parameter :: array_rank = 4                   !For this version it is 4
    !Calling parameters of the subroutine
    real(kind=RP), dimension(:,:,:,:), intent(in) :: array !4D real array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_real_4d


  !Write to binary for 5D array with real values
  subroutine array_to_binary_real_5d(array, file_name, funit, error)
    integer, parameter :: array_rank = 5                    !For this version it is 5
    !Calling parameters of the subroutine
    real(kind=RP), dimension(:,:,:,:,:), intent(in) :: array !5D real array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_real_5d


  !Write to binary for 1D array with integer values
  subroutine array_to_binary_int_1d(array, file_name, funit, error)
    integer, parameter :: array_rank = 1       !For this version it is 1
    !Calling parameters of the subroutine
    integer, dimension(:), intent(in) :: array !1D integer array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_int_1d


  !Write to binary for 2D array with integer values
  subroutine array_to_binary_int_2d(array, file_name, funit, error)
    integer, parameter :: array_rank = 2         !For this version it is 2
    !Calling parameters of the subroutine
    integer, dimension(:,:), intent(in) :: array !2D integer array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_int_2d


  !Write to binary for 3D array with integer values
  subroutine array_to_binary_int_3d(array, file_name, funit, error)
    integer, parameter :: array_rank = 3           !For this version it is 3
    !Calling parameters of the subroutine
    integer, dimension(:,:,:), intent(in) :: array !3D integer array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_int_3d


  !Write to binary for 4D array with integer values
  subroutine array_to_binary_int_4d(array, file_name, funit, error)
    integer, parameter :: array_rank = 4             !For this version it is 4
    !Calling parameters of the subroutine
    integer, dimension(:,:,:,:), intent(in) :: array !4D integer array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_int_4d


  !Write to binary for 5D array with integer values
  subroutine array_to_binary_int_5d(array, file_name, funit, error)
    integer, parameter :: array_rank = 5               !For this version it is 5
    !Calling parameters of the subroutine
    integer, dimension(:,:,:,:,:), intent(in) :: array !5D array
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    character(255) :: error_message

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if
    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='replace', iostat=error_status,             &
          iomsg=error_message )
    !If there was no error opening the file, proceed to write it.
    if (error_status == 0) then
      dim_length(1:array_rank) =  shape(array)
      write(funit_loc) array_rank, (dim_length(i), i = 1, array_rank), array
    else
      write(ERROR_UNIT,*) "Error writing the file <<",file_name, ">>: ",      &
                          error_message
    end if
    close(funit_loc)
    !Report error optionally.
    if (present(error)) then
      error = error_status
    end if

  end subroutine array_to_binary_int_5d

  !If you want to print another array type, just add a new submodule similar
  ! to array_to_binary_real_5d or array_to_binary_int_5d :)

  !-----------------------------------------------------------------------------
  ! READ BINARY FILE TO ARRAY
  !-----------------------------------------------------------------------------

  ! ``fortran
  ! ! Here is a short test program that illustrates how to
  ! ! use the binary read procedures:
  ! program test
  !
  !     use binary_IO
  !
  !     integer, dimension(5,15,25,35,10) :: ZZ = -9999
  !     integer, allocatable, dimension(:,:,:,:,:) :: II
  !     integer :: ier
  !
  !     ! Write data to a binary file.
  !     call array_to_binary(ZZ, "zzz.bin")
  !
  !     ! Use the diagnostic functions to check what the files are like.
  !     print *, "RANK = ", binary_get_rank("zzz.bin")
  !     print *, "DIMS = ", binary_get_dims("zzz.bin")
  !
  !     ! Print the whole data array halved using the read function:
  !     print *, binary_read_i5d("zzz.bin") / 2
  !
  !     ! Use the read function in another way
  !     II= binary_read_i5d("zzz.bin")
  !     print *, II
  !
  !     ! Use the read subroutine:
  !     call binary_to_array(II, "zzz.bin", error=ier)
  !     print *, II
  !
  ! end program test
  ! ```



  ! Read binary file into a 1D array with real values.
  subroutine binary_to_array_real_1d(array, file_name, funit, error)
    integer, parameter :: array_rank = 1               !For this version it is 1
    !Return value is actually a 1D array
    real(kind=RP), dimension(:), intent(inout) :: array !1D real array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                          FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_real_1d

 !Read binary file into a 2D array with real values.
  subroutine binary_to_array_real_2d(array, file_name, funit, error)
    integer, parameter :: array_rank = 2               !For this version it is 2
    !Return value is actually a 2D array
    real(kind=RP), dimension(:,:), intent(inout) :: array !2D real array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if

        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_real_2d

 !Read binary file into a 3D array with real values.
  subroutine binary_to_array_real_3d(array, file_name, funit, error)
    integer, parameter :: array_rank = 3               !For this version it is 3
    !Return value is actually a 3D array
    real(kind=RP), dimension(:,:,:), intent(inout) :: array !3D real array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_real_3d

  !Read binary file into a 4D array with real values.
  subroutine binary_to_array_real_4d(array, file_name, funit, error)
    integer, parameter :: array_rank = 4               !For this version it is 4
    !Return value is actually a 4D array
    real(kind=RP), dimension(:,:,:,:), intent(inout) :: array !4D real array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_real_4d

  !Read binary file into a 5D array with real values.
  subroutine binary_to_array_real_5d(array, file_name, funit, error)
    integer, parameter :: array_rank = 5               !For this version it is 5
    !Return value is actually a 5D array
    real(kind=RP), dimension(:,:,:,:,:), intent(inout) :: array !5D real array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_real_5d

  !Read binary file into a 1D array with integer values.
  subroutine binary_to_array_int_1d(array, file_name, funit, error)
    integer, parameter :: array_rank = 1              !For this version it is 1
    !Return value is actually a 1D array
    integer, dimension(:), intent(inout) :: array     !1D integer array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_int_1d

  !Read binary file into a 2D array with integer values.
  subroutine binary_to_array_int_2d(array, file_name, funit, error)
    integer, parameter :: array_rank = 2            !For this version it is 2
    !Return value is actually a 2D array
    integer, dimension(:,:), intent(inout) :: array !2D integer array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_int_2d

  !Read binary file into a 3D array with integer values.
  subroutine binary_to_array_int_3d(array, file_name, funit, error)
    integer, parameter :: array_rank = 3              !For this version it is 3
    !Return value is actually a 3D array
    integer, dimension(:,:,:), intent(inout) :: array !3D integer array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_int_3d

  !Read binary file into a 4D array with integer values.
  subroutine binary_to_array_int_4d(array, file_name, funit, error)
    integer, parameter :: array_rank = 4               !For this version it is 4
    !Return value is actually a 4D array
    integer, dimension(:,:,:,:), intent(inout) :: array !4D integer array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_int_4d

  !Read binary file into a 5D array with integer values.
  subroutine binary_to_array_int_5d(array, file_name, funit, error)
    integer, parameter :: array_rank = 5              !For this version it is 5
    !Return value is actually a 5D array
    integer, dimension(:,:,:,:,:), intent(inout) :: array !5D array
    !Calling parameters of the subroutine
    character(len=*), intent(in) :: file_name
    integer, optional, intent(in) :: funit
    integer, optional, intent(out) :: error
    !Local variables
    integer, dimension(array_rank) :: dim_length
    integer :: i, funit_loc, error_status
    integer :: array_rank_get

    !Check if optional unit name is provided, if yes, use it if not use
    ! the first free unit. This avoids possible units conflicts if many files
    ! are opened in parallel.
    if (present(funit)) then
      funit_loc=funit
    else
      funit_loc=GET_FREE_FUNIT()
    end if

    !Open the files and write. Note that local variable funit_loc keeps unit.
    open( unit=funit_loc, file=trim(FIX_FOLDER) // file_name, form=file_form, &
          access='stream', status='old', iostat=error_status )
    !If there was no error opening the file, proceed to read it.
    if (error_status == 0) then
      !First, get the rank of the array that is saved in the binary file.
      read(unit=funit_loc) array_rank_get
      !Here check if the rank from the file coincides with the input array
      if ( array_rank_get /= array_rank ) then
        write(ERROR_UNIT,*) "Wrong rank of the parameter array in file ",     &
                            array_rank_get, ", but must be ", array_rank,     &
                            ", fake array returned: ", FAKE_ERR
        array = FAKE_ERR
        if (present(error)) then
          error = ERROR_RANK_MISMATCH
        end if
        close(funit_loc)
        return
      end if
      !Second, get the shape for the dimensions rank of the array
      read(unit=funit_loc) (dim_length(i), i = 1, array_rank_get)
      !Third, read the actual values of the data array now.
      read(unit=funit_loc) array
    !If array read is not successful...
    else
      ! ... return a fake array with clearly unusual values, so the error
      ! is easy to notice in debugging.
      write(ERROR_UNIT,*) "Error reading the file <<",file_name, ">>, ",      &
                           FAKE_ERR, " fake array returned."
      array = FAKE_ERR
    end if
    !Report error optionally
    if (present(error)) then
      error = error_status
    end if
    !Finally, close the file.
    close(funit_loc)

  end subroutine binary_to_array_int_5d

end module binary_IO
