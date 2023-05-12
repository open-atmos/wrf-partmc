program make_emissions_ideal

  use netcdf
  use pmc_netcdf
  use pmc_aero_dist
  use pmc_aero_mode
  use pmc_aero_data
  use pmc_spec_file
  use pmc_gas_state
  use pmc_mpi
  use pmc_scenario
  use mpi

  implicit none
  
  !> Filename of NetCDF file to open.
  character(len=200) :: filename_grid
  character(len=200) :: filename_emissions
  !> NetCDF file ID, in data mode.
  real(kind=dp), allocatable, dimension(:,:,:,:,:) :: gas_emissions, &
       gas_emissions_patch
  real(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: aero_emissions, &
       aero_emissions_patch
  character(len=100) :: name
  integer :: varid, dimid_times, send_size_a, send_size_g, ierr
  integer :: status, dimid_n_gas_specs
  integer :: var_size(4)
  integer :: temp_size(2)
  integer :: dim
  integer :: t
  integer :: n_source_classes
  real(kind=dp) :: total_mass
  real(kind=dp) :: diam_i, diam_j
  real(kind=dp) :: std_i, std_j
  real(kind=dp) :: num_conc_i, num_conc_j
  integer :: x, y, z, nx, ny, nz, nt, s
  integer :: i,j,k, i_source, n_species, i_mode, iclass, i_type
  integer :: n_gas_spec_partmc, n_aero_spec_partmc
  type(aero_data_t) :: aero_data
  type(gas_data_t) :: gas_data
  character(len=9), allocatable :: partmc_species(:)
  character(len=9), allocatable :: smoke_species(:)
  character(len=9), allocatable :: wrf_species(:)

  ! output structures
  type(aero_dist_t), allocatable :: aero_emission(:)
  real(kind=dp), allocatable, dimension(:) :: aero_emission_rate, &
       aero_emission_time
  real(kind=dp), allocatable, dimension(:) :: gas_emission_rate, &
       gas_emission_time
  integer :: n_mode_counter
  real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std, &
       mode_num_concs
  integer, allocatable, dimension(:) ::  mode_source
  real(kind=dp), allocatable, dimension(:,:) :: mode_vol_fracs, temp_array, &
       mode_vol_frac_std
  real(kind=dp), allocatable, dimension(:) :: vol_frac
  type(gas_state_t), allocatable :: gas_emission(:)
  character(len=100) :: prefix, filename, case_name, dir_path
  character(len=3), parameter  :: suffix = '.nc'
  integer :: ncid
  character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
  type(spec_file_t) :: sub_file, file
  character(len=300) :: spec_name
  character(len=300) :: grid_name
  integer :: start_date, n_days_per_file, n_days, i_day
  logical :: do_point, do_nonpoint, do_bio, do_nonroad, do_mobile
  logical :: only_sectional
  real(kind=dp) :: gas_scale_factor, aerosol_scale_factor
  real(kind=dp) :: dx

  integer :: i_send_s, i_send_e
  integer :: rem, is_local, ie_local, n_proc, i_proc, root
  integer :: recv_size_a, recv_size_g
  integer, allocatable, dimension(:) :: displs_a, send_counts_a
  integer, allocatable, dimension(:) :: displs_g, send_counts_g

!  character(len=200),parameter,dimension(4) :: mobile_sectors = &
!      ["RPH", "RPD", "RPV", "RPP"]
  character(len=200),parameter,dimension(3) :: mobile_sectors = &
      ["RPH", "RPV", "RPD"]

  integer :: stride, block_length
  integer :: grid_emission
  integer :: n_smoke_species

  ! Use this one for emissions
  ! FIXME: Limited to the surface
  type(scenario_t) :: scenario
 
  type(aero_dist_t), allocatable, dimension(:) :: aero_dist_init, aero_dist
  type(gas_state_t), allocatable, dimension(:) :: gas_state_init
  type(gas_state_t) :: gas_state
  real(kind=dp), allocatable, dimension(:) :: altitudes
  integer :: i_level, n_levels
  integer :: start(1), count(1)
  integer :: start2(2), count2(2)
  integer :: start3(3), count3(3)

  call pmc_mpi_init()

  if (command_argument_count() /= 1) then
        call print_usage()
        call die_msg(739173198, "invalid commandline arguments")
  end if

  call get_command_argument(1, spec_name)

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

  ! Read in emissions
  call spec_file_read_string(file, "aero_emissions", sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_dists_times_rates(sub_file, aero_data, .false., &
       scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
       scenario%aero_emission)
  call spec_file_close(sub_file)

  ! gas emissions profile
  call spec_file_read_string(file, "gas_emissions", sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
       scenario%gas_emission_time, scenario%gas_emission_rate_scale, &
       scenario%gas_emission)
  call spec_file_close(sub_file)

  call spec_file_read_integer(file, 'nz', nz)

  n_aero_spec_partmc = aero_data_n_spec(aero_data)
  n_gas_spec_partmc = gas_data_n_spec(gas_data)

  ! FIXME: We should set this somehow and then put it in namelist.
  ! Output to the NetCDF file
  i = 2
  j = 2
  k = 1
  nt = size(scenario%gas_emission_time)
  prefix = "ideal/aero_emit_dist_"
  write(filename, '(a,i3.3,a,i3.3,a,i3.3,a)') trim(prefix),(i), &
       '_',j,'_',k,trim(suffix)
  call pmc_nc_open_write((filename), ncid)
  call pmc_nc_check(nf90_redef(ncid))
  call pmc_nc_check(nf90_def_dim(ncid, "n_times", &
       nt, dimid_times))
  call pmc_nc_check(nf90_enddef(ncid))
  call pmc_nc_write_real_1d(ncid, scenario%aero_emission_rate_scale, &
       "aero_emission_rate_scale", (/ dimid_times /), &
       unit="(1)", &
       long_name="Aerosol emission rate", &
       description="Aerosol emission rate scales at set-points")
  call pmc_nc_write_real_1d(ncid, scenario%aero_emission_time, &
       "aero_emission_time", (/ dimid_times /), &
       unit="s", &
       long_name="Aerosol emission time", &
       description="Aerosol emission set-points times (s).")
  call aero_dist_output_netcdf(scenario%aero_emission, ncid)

  ! Output the gas emissions
  call pmc_nc_write_real_1d(ncid, scenario%gas_emission_rate_scale, &
        "gas_emission_rate_scale", (/ dimid_times /), &
        unit="(1)", &
        long_name="Gas emission rate", &
        description="Gas emission rate scales at set-points")
  call pmc_nc_write_real_1d(ncid, scenario%gas_emission_time, &
        "gas_emission_time", (/ dimid_times /), &
        unit="s", &
        long_name="Gas emission time", &
        description="Gas emission set-points times (s).")
  call gas_emissions_output_netcdf(scenario%gas_emission, ncid)

  call pmc_nc_check(nf90_close(ncid))

  ! FIXME: Gases are specified at some levels. This is a short cut til we decide
  ! what we want to do.
  n_levels = 3
  allocate(gas_state_init(n_levels))
  allocate(altitudes(n_levels))
  do i_level = 1, n_levels
     call spec_file_read_string(file, 'gas_init', sub_filename)
     call spec_file_open(sub_filename, sub_file)
     call spec_file_read_gas_state(sub_file, gas_data, gas_state_init(i_level))
  end do

  ! FIXME: We can make this variable with height somehow
  allocate(aero_dist_init(1)) 
  call spec_file_read_string(file, 'aerosol_init', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_dist(sub_file, aero_data, .false., aero_dist_init(1))
  call spec_file_close(sub_file)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subroutine to output to a file the aerosol emission
  subroutine aero_dist_output_netcdf(aero_emissions, ncid)

    !> Aerosols emissions for all time.
    type(aero_dist_t), allocatable, intent(in) :: aero_emissions(:)
    !> NetCDF file ID.
    integer, intent(in) :: ncid

    character(len=100) :: unit, long_name, standard_name, description
    integer :: dimid_n_times, dimid_n_modes, dimid_n_specs
    integer :: n_modes, n_times, n_specs
    integer :: status, check_dim_size
    character(len=NF90_MAX_NAME) :: check_name
    integer, allocatable, dimension(:) :: source
    real(kind=dp), allocatable, dimension(:) :: radius, std
    real(kind=dp), allocatable, dimension(:,:) :: num_conc
    real(kind=dp), allocatable, dimension(:,:,:) :: vol_frac
    integer :: i_time, i_mode

    integer :: start(1), count(1)
    integer :: start2(2), count2(2)
    integer :: start3(3), count3(3)

    ! First we have to define some dimensions
    ! Time was defined already, so let's find it again
    status = nf90_inq_dimid(ncid, "n_times", dimid_n_times)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_times, check_name, &
            check_dim_size))
    n_times = check_dim_size
    ! number of modes isn't defined
    n_modes = size(aero_emissions(1)%mode) 
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "n_modes", &
         n_modes, dimid_n_modes))
    call pmc_nc_check(nf90_enddef(ncid))
    !
    n_specs = size(aero_emissions(1)%mode(1)%vol_frac)
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
       radius(i_mode) = aero_emissions(1)%mode(i_mode)%char_radius
    end do

    start = (/ 1 /)
    count = (/ n_modes/)
    call pmc_nc_check(nf90_put_var(ncid, varid, radius, &
         start = start, count = count))

    ! 
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
       std(i_mode) = aero_emissions(1)%mode(i_mode)%log10_std_dev_radius
    end do

    start = (/ 1 /)
    count = (/ n_modes/)
    call pmc_nc_check(nf90_put_var(ncid, varid, std, &
         start = start, count = count))
  
    ! 
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
       source(i_mode) = aero_emissions(1)%mode(i_mode)%source
    end do

    start = (/ 1 /)
    count = (/ n_modes/)
    call pmc_nc_check(nf90_put_var(ncid, varid, source, &
         start = start, count = count))   
    ! Things that are 2D
    !   - number concentration (times, n_modes)
    name='num_conc'
    unit='# m^{-3}'
    long_name = 'total number concentration'
    standard_name = 'total number concentration'
    description = 'Total number concentration of mode (#/m^3).'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, &
         (/dimid_n_modes, dimid_n_times/), varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))
    allocate(num_conc(n_modes, n_times))
    do i_time = 1, n_times
    do i_mode = 1, n_modes
       num_conc(i_mode,i_time) = aero_emissions(i_time)%mode(i_mode)%num_conc
    end do
    end do
    start2 = (/ 1, 1 /)
    count2 = (/ n_modes, n_times/)
    call pmc_nc_check(nf90_put_var(ncid, varid, num_conc, &
         start = start2, count = count2))
    ! Things that are 3D
    !   - vol_frac (times, n_modes, n_spec)
    name='vol_frac'
    unit='(1)'
    long_name = 'species fractions'
    standard_name = 'species_fractions'
    description = 'Species fractions by volume [length \c aero_data%%n_spec].'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, &
         (/dimid_n_specs, dimid_n_modes, dimid_n_times/), varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))    

    allocate(vol_frac(n_specs, n_modes, n_times))
    do i_time = 1, n_times
    do i_mode = 1, n_modes
       vol_frac(:,i_mode,i_time) = aero_emissions(i_time)%mode(i_mode)%vol_frac
    end do
    end do

    start3 = (/1,1,1/)
    count3 = (/n_specs, n_modes, n_times/)
    call pmc_nc_check(nf90_put_var(ncid, varid, vol_frac, &
         start = start3, count = count3))

  end subroutine aero_dist_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a number concentration based off a log-normal distribution with
  !! a given diameter and standard deviation
  real(kind=dp) function get_num_conc(mass, diam, std)

    !> Mass concentrations.
    real(kind=dp) :: mass
    !> Geometric mean diameter.
    real(kind=dp) :: diam
    !> Geometric standard deviation.
    real(kind=dp) :: std

    integer :: k
    real(kind=dp) :: density, tmp

    k = 3
    ! FIXME: This needs the be handled differently as some species have
    ! different densities.
    density = 1800.0d0 ! kg / m^3
    density = density * 1000.0d0 ! g / m^3
    tmp = density*(const%pi/ 6.0d0)*diam**3.0 * &
         exp(k**2.0d0 / 2.0d0 * log(std)**2.0d0)
    get_num_conc = mass / tmp

  end function get_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine get_vol_frac(species, vol_frac, n_specs, aero_data)

    !> Masses for all species.
    real(kind=dp) :: species(:)
    !> Volume fractions.
    real(kind=dp), allocatable :: vol_frac(:)
    !> Number of aerosol species.
    integer :: n_specs
    !> Aerosol data.
    type(aero_data_t) :: aero_data

    integer :: i

    integer :: spec_index(5)
    real(kind=dp) :: tot_vol_frac

    if (allocated(vol_frac)) deallocate(vol_frac)
    !if (allocated(vol_frac_std)) deallocate(vol_frac_std)
    allocate(vol_frac(n_specs))
    !allocate(vol_frac_std(n_specs))
    vol_frac = 0.0d0

    if (sum(species) > 0.0d0) then
       ! FIXME: Check ordering 
       ! Map the input species to PartMC species
       ! SO4 = 1
       ! NO3 = 2
       ! OIN = 17
       ! OC = 18
       ! BC = 19
       !spec_index = [1,19,2,18,17]
       spec_index = [17,1,2,18,19]
       do i = 1, size(species)
          vol_frac(spec_index(i)) = species(i)
       end do
       vol_frac = vol_frac / aero_data%density
       tot_vol_frac = sum(vol_frac)

       ! convert mass fractions to volume fractions
       vol_frac = vol_frac / tot_vol_frac
       !vol_frac_std = vol_frac_std / aero_data%density
   else ! There isn't anything for this mode so set a dummy value
       vol_frac = 1.0d0
   end if

  end subroutine get_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs gas emissions for all time.
  subroutine gas_emissions_output_netcdf(gas_emissions, ncid)

    !> Gas states containing emissions for all time.
    type(gas_state_t) :: gas_emissions(:)
    !> NetCDF file ID.
    integer :: ncid

    character(len=100) :: unit, long_name, standard_name, description
    integer :: dimid_n_times, dimid_n_gas_specs
    integer :: n_times, n_gas_specs
    integer :: status, check_dim_size
    character(len=NF90_MAX_NAME) :: check_name
    real(kind=dp), allocatable, dimension(:,:) :: mix_rat 
    integer :: i_time
    integer :: varid
    integer :: start2(2), count2(2)

    ! First we have to define some dimensions
    ! Time was defined already, so let's find it again
    status = nf90_inq_dimid(ncid, "n_times", dimid_n_times)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_times, check_name, &
            check_dim_size))
    n_times = check_dim_size

    n_gas_specs = size(gas_emissions(1)%mix_rat)
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "n_gas_specs", &
         n_gas_specs, dimid_n_gas_specs))
    call pmc_nc_check(nf90_enddef(ncid))

    ! Things that are 2D
    name='gas_emit'
    unit='mol m^{-2} s^{-1}'
    long_name = 'gas emissions'
    standard_name = 'gas emissions'
    description = 'gas phase emissions.'
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, &
         (/dimid_n_gas_specs, dimid_n_times/), varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))
    allocate(mix_rat(n_gas_specs, n_times))
    do i_time = 1,n_times
       mix_rat(:,i_time) = gas_emissions(i_time)%mix_rat
    end do
    start2 = (/ 1, 1 /)
    count2 = (/ n_gas_specs, n_times/)
    call pmc_nc_check(nf90_put_var(ncid, varid, mix_rat, &
         start = start2, count = count2))


  end subroutine gas_emissions_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the usage text to stderr.
  subroutine print_usage()

    write(*,*) 'Usage: ./create_emissions <spec-file>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program make_emissions_ideal
