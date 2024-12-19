module ic_bc_helper

  use pmc_netcdf
  use pmc_mpi
  use pmc_aero_data

  integer, parameter :: num_mode_mam3 = 3
  ! Map MAM3 modes to PartMC modes
  integer, dimension(3), parameter :: pmc_to_mam3_map = (/ 2, 1, 3 /)
  integer, parameter :: num_spec_mam3 = 6
  character(len=3), dimension(num_spec_mam3), parameter :: aero_names_mam3 = &
       ["so4", "ncl", "bc ", "pom", "soa", "dst"]
  character(len=4), dimension(num_spec_mam3), parameter :: aero_names_pmc = &
       ["SO4 ","Na  ", "BC  ", "OC  ", "API1", "OIN "]
  real(kind=dp), parameter :: mw_na = 22.990d0
  real(kind=dp), parameter :: mw_cl = 35.450d0
  real(kind=dp), dimension(num_spec_mam3), parameter :: aero_factors = &
       [96.0d0/115.0d0, mw_na / (mw_na + mw_cl), 1.0d0, 1.0d0, 1.0d0, 1.0d0]
  ! Species contained in each MAM3 mode
  logical, dimension(num_mode_mam3, num_spec_mam3), parameter :: &
      mode_contains_species = &
      reshape([.true., .true., .true., .true., .true., .true., &
      .true., .false., .false., .true., .false., .false., &
      .true., .true., .false., .true., .false., .true.], &
      shape(mode_contains_species))

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Makes a string all lowercase.
  function make_lower(strIn)

    !> Input string
    character(len=*), intent(in) :: strIn

    character(len=len(strIn)) :: make_lower
    integer :: i,j

    do i = 1, len(strIn)
       j = iachar(strIn(i:i))
       if (j<= iachar("9")) then
          make_lower(i:i) = strIn(i:i)
       else if (j>= iachar("a") .and. j<=iachar("z") ) then
          make_lower(i:i) = strIn(i:i)
       else
          make_lower(i:i) = achar(iachar(strIn(i:i))+32)
       end if
    end do

  end function make_lower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple real array from a NetCDF file.
  subroutine pmc_nc_read_string_array(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    character(len=19), allocatable, dimension(:), intent(inout) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present
    integer :: dimids(2)
    integer :: dim_size(2)
    integer :: start(2)
    integer :: count(2)
    integer :: i, dim

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = ""
       return
    end if

    ! Find the varid
    status = nf90_inq_varid(ncid, name, varid)
    ! Inquire for the dimension IDs
    call pmc_nc_check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    ! Loop over all the dimensions to get their lengths
    do dim = 1, 2
       call pmc_nc_check(nf90_inquire_dimension(ncid, dimids(dim), &
            len=dim_size(dim)))
    end do

    allocate(var(dim_size(2)))
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))

  end subroutine pmc_nc_read_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_mode_parameters(aero_mode_name, mode_diams, mode_std, &
       mode_vol_fracs, mode_source, aero_data)

    character(len=100), allocatable, dimension(:) :: aero_mode_name
    real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std
    real(kind=dp), allocatable, dimension(:,:) :: mode_vol_fracs
    integer, allocatable, dimension(:) :: mode_source
    type(aero_data_t) :: aero_data

    integer :: n_modes
    integer :: n_aero_species

    n_modes = 3
    n_aero_species = aero_data_n_spec(aero_data)

    allocate(mode_diams(n_modes))
    allocate(mode_std(n_modes))
    allocate(mode_source(n_modes))
    allocate(aero_mode_name(n_modes))
    allocate(mode_vol_fracs(n_modes,n_aero_species))

    aero_mode_name(1) = "aitken"
    mode_diams(1) = 0.0260d-6
    mode_std(1) =  1.6d0
    mode_source(1) = 1

    aero_mode_name(2) = "accumulation"
    mode_diams(2) =  0.11d-6
    mode_std(2) =  1.8d0
    mode_source(2) = 2

    aero_mode_name(3) = "coarse"
    mode_diams(3) =  2.0d-6
    mode_std(3) =  1.8d0
    mode_source(3) = 3

  end subroutine get_mode_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ic_bc_helper
