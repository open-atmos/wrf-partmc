program make_ics

  use netcdf
  use pmc_netcdf
  use pmc_aero_dist
  use pmc_aero_mode
  use pmc_aero_data
  use pmc_spec_file
  use pmc_gas_data
  use pmc_gas_state
  use pmc_mpi
  use mpi

  implicit none
 
  !> Filename of NetCDF file to open.
  character(len=100) :: filename_ic
  !> NetCDF file ID, in data mode.
  integer :: ncid_ic

  character(len=100), allocatable, dimension(:) :: aero_spec_name
  character(len=100), allocatable, dimension(:) :: gas_spec_name

  real(kind=dp), allocatable, dimension(:,:,:,:) :: aero_values, gas_values
  real(kind=dp), allocatable, dimension(:,:,:) :: density
  character(len=100) :: name
  integer :: varid, dimid_times
  integer :: status
  integer :: dim
  integer :: n_modes
  integer :: nx, ny, nz, nt
  integer :: i,j,k
  integer :: i_mode
  type(aero_data_t) :: aero_data
  type(gas_data_t) :: gas_data
  type(gas_state_t) :: gas_state
  type(gas_state_t), allocatable, dimension(:) :: gas_state_init
  type(aero_dist_t) :: aero_dist
  type(aero_dist_t), allocatable, dimension(:) :: aero_dist_init
  integer :: ncid
  real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std, &
       mode_num_concs
  integer, allocatable, dimension(:) :: mode_source
  real(kind=dp), allocatable, dimension(:,:) :: mode_vol_fracs
  character(len=300) :: output_file_prefix, output_file_path, file_prefix, &
       prefix, sub_filename
  character(len=100) ::  suffix
  type(spec_file_t) :: sub_file, file
  character(len=300) :: file_path, filename, command, spec_name
  integer :: start(1), count(1)
  ! 
  integer :: rem, n_proc, is_local, ie_local, dimid_n_gas_specs
  integer :: n_gas_specs, n_aero_specs

  call pmc_mpi_init()

  if (pmc_mpi_rank() == 0) then
     call get_command_argument(1, spec_name)
  end if

  call spec_file_open(spec_name, file)

  ! read in aero_data
  call spec_file_read_string(file, 'aerosol_data', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_data(sub_file, aero_data)
  call spec_file_close(sub_file)

  ! read in gas_data
  call spec_file_read_string(file, 'gas_data', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_gas_data(sub_file, gas_data)
  call spec_file_close(sub_file)

  allocate(aero_dist_init(1))
  allocate(gas_state_init(1))

  call spec_file_read_string(file, 'aerosol_init', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_dist(sub_file, aero_data, aero_dist_init(1))
  call spec_file_close(sub_file)

  call spec_file_read_string(file, 'gas_init', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_gas_state(sub_file, gas_data, &
       gas_state_init(1))
  call spec_file_close(sub_file)

  suffix = '.nc'

  call spec_file_read_integer(file, 'nz', nz)

  n_aero_specs = aero_data_n_spec(aero_data)
  n_gas_specs = gas_data_n_spec(gas_data)

  prefix = "ideal/ics_"
  do i = 1,4
  do j = 1,4
  do k = 1,nz
    ! Output to the NetCDF file
    write(filename, '(a,i3.3,a,i3.3,a,i3.3,a)') trim(prefix),(i), &
       '_',j,'_',k,trim(suffix)
    call pmc_nc_open_write(filename, ncid)
    aero_dist = aero_dist_init(1)
    ! Output the boundary condition data
    call aero_dist_output_netcdf(aero_dist,ncid)
    gas_state = gas_state_init(1)

    ! Output gas state
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "n_gas_specs", &
         n_gas_specs, dimid_n_gas_specs))
    call pmc_nc_check(nf90_enddef(ncid))
    name='gas_ic'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, &
         (/dimid_n_gas_specs/), varid))
    call pmc_nc_write_atts(ncid, varid, unit='ppb', &
         long_name='gas initial condition', &
         standard_name='gas boundary condition', &
         description='value of gas at boundary')
    start = (/ 1/)
    count = (/n_gas_specs/)
    call pmc_nc_check(nf90_put_var(ncid, varid, gas_state%mix_rat, &
         start = start, count = count))

    call pmc_nc_check(nf90_close(ncid))
  end do
  end do
  end do

  call pmc_mpi_finalize()
 
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subroutine to output to a file the aerosol distributions.
  subroutine aero_dist_output_netcdf(aero_dists,ncid)

    !> Array of aerosol distributions.
    type(aero_dist_t) :: aero_dists
    !> NetCDF file ID, in data mode
    integer, intent(in) :: ncid

    character(len=100) :: unit, long_name, standard_name, description
    integer :: dimid_n_times, dimid_n_modes, dimid_n_specs
    integer :: n_modes, n_times, n_specs
    integer :: status, check_dim_size
    character(len=NF90_MAX_NAME) :: check_name
    integer, allocatable, dimension(:) :: source
    real(kind=dp), allocatable, dimension(:) :: radius, std
    real(kind=dp), allocatable, dimension(:) :: num_conc
    real(kind=dp), allocatable, dimension(:,:) :: vol_frac
    integer :: i_time, i_mode
    integer :: start(1), count(1)
    integer :: start2(2), count2(2)
    integer :: start3(3), count3(3)

    ! number of modes isn't defined
    n_modes = size(aero_dists%mode)
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "n_modes", &
         n_modes, dimid_n_modes))
    call pmc_nc_check(nf90_enddef(ncid))
    !
    n_specs = size(aero_dists%mode(1)%vol_frac)
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "n_aero_specs", &
         n_specs, dimid_n_specs))
    call pmc_nc_check(nf90_enddef(ncid))

    ! We then create some variables with descriptions
    ! Things that are 1D
    !   - diameter (n_modes)
    !   - names (n_modes)
    !       - names should be done similar to aero_data species names
    !   - standard deviations (n_modes)
    !
    name='char_radius'
    unit='m'
    long_name = 'characteristic_radius'
    standard_name = 'characteristic_radius'
    description = 'Characteristic radius, with meaning dependent on mode type'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid_n_modes, &
         varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))
    allocate(radius(n_modes))
    do i_mode = 1, n_modes
       radius(i_mode) = aero_dists%mode(i_mode)%char_radius
    end do

    start = (/ 1 /)
    count = (/ n_modes /)
    call pmc_nc_check(nf90_put_var(ncid, varid, radius, &
         start = start, count = count))

    name='log10_std_dev_radius'
    unit='m'
    long_name = 'log10_std_dev_radius' 
    standard_name = 'log10_std_dev_radius' 
    description = 'Log base 10 of geometric standard deviation of radius, (m).'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimid_n_modes, &
         varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))
    allocate(std(n_modes))
    do i_mode = 1, n_modes
       std(i_mode) = aero_dists%mode(i_mode)%log10_std_dev_radius
    end do

    start = (/ 1 /)
    count = (/ n_modes /)
    call pmc_nc_check(nf90_put_var(ncid, varid, std, &
         start = start, count = count))
  
    name='source'
    unit='(1)'
    long_name = 'Source number.'
    standard_name = 'Source number'
    description = 'Source number.'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimid_n_modes, &
         varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))
    allocate(source(n_modes))
    do i_mode = 1, n_modes
       source(i_mode) = aero_dists%mode(i_mode)%source
    end do

    start = (/ 1 /)
    count = (/ n_modes /)
    call pmc_nc_check(nf90_put_var(ncid, varid, source, &
         start = start, count = count))   
    ! Things that are 1D
    !   - number concentration (n_modes)
    name='num_conc'
    unit='#m/^3'
    long_name = 'total number concentration'
    standard_name = 'total number concentration'
    description = 'Total number concentration of mode (#/m^3).'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, &
         (/dimid_n_modes/), varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))
    allocate(num_conc(n_modes))
    do i_mode = 1,n_modes
       num_conc(i_mode) = aero_dists%mode(i_mode)%num_conc
    end do

    ! Error checking
    call assert_msg(754432751, sum(num_conc) > 0.0d0, &
         "number concentrations are zero. check input files")

    start = (/ 1 /)
    count = (/ n_modes/)
    call pmc_nc_check(nf90_put_var(ncid, varid, num_conc, &
         start = start, count = count))
    ! Things that are 2D
    !   - vol_frac (n_modes, n_spec)
    name='vol_frac'
    unit='(1)'
    long_name = 'species fractions'
    standard_name = 'species_fractions'
    description = 'Species fractions by volume [length \c aero_data%%n_spec].'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, &
         (/dimid_n_specs, dimid_n_modes/), varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))    

    allocate(vol_frac(n_specs, n_modes))
       do i_mode = 1,n_modes
          vol_frac(:,i_mode) = aero_dists%mode(i_mode)%vol_frac
       end do

    start2 = (/1,1/)
    count2 = (/n_specs, n_modes/)
    call pmc_nc_check(nf90_put_var(ncid, varid, vol_frac, &
         start = start2, count = count2))

  end subroutine aero_dist_output_netcdf

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

  !> Print the usage text to stderr.
  subroutine print_usage()

    write(*,*) 'Usage: ./create_ics <wrf file path> <outfile path> &
        <outfile prefix>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program make_ics
