program make_ics

  use netcdf
  use pmc_netcdf
  use pmc_aero_dist
  use pmc_aero_mode
  use pmc_aero_data
  use pmc_spec_file
  use pmc_mpi
  use mpi

  implicit none
 
  !> Filename of NetCDF file to open.
  character(len=100) :: filename_ic
  !> NetCDF file ID, in data mode.
  integer :: ncid_ic

  character(len=100), allocatable, dimension(:) :: aero_spec_name

  real(kind=dp), allocatable, dimension(:,:,:,:,:) :: aero_values
  real(kind=dp), allocatable, dimension(:,:,:) :: density
  character(len=100) :: name
  integer :: varid, dimid_times
  integer :: status
  integer :: dim
  integer :: n_modes
  integer :: nx, ny, nz, nt
  integer :: i,j,k
  integer :: i_mode
  integer :: n_aero_species
  type(aero_data_t) :: aero_data
  integer :: is, ie, js, je
  real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std, &
       mode_num_concs
  integer, allocatable, dimension(:) :: mode_source
  real(kind=dp), allocatable, dimension(:,:) :: mode_vol_fracs
  character(len=300) :: output_file_prefix, output_file_path, file_prefix
  character(len=100) ::  suffix
  type(spec_file_t) :: sub_file
  character(len=300) :: file_path, filename, command

  ! 
  integer :: rem, n_proc, is_local, ie_local

  logical, parameter :: mam3 = .true.

  call pmc_mpi_init()

  if (pmc_mpi_rank() == 0) then
     ! only the root process accesses the commandline

     if (command_argument_count() /= 3) then
        call print_usage()
        call die_msg(739173192, "invalid commandline arguments")
     end if

     call get_command_argument(1, file_path)
     call get_command_argument(2, output_file_path)
     call get_command_argument(3, output_file_prefix)
  end if

  call pmc_mpi_bcast_string(file_path)
  if (pmc_mpi_rank() == 0) then
     write(command,'(A,A)') "mkdir ", trim(output_file_path)
     call system(command)
     write(file_prefix,'(A,A,A)') trim(output_file_path), '/', &
          trim(output_file_prefix)
  end if

  call pmc_mpi_bcast_string(file_prefix)

  ! read in aero_data
  write(filename,'(A,A)') trim(file_path),'aero_data.dat'
  call spec_file_open(filename, sub_file)
  call spec_file_read_aero_data(sub_file, aero_data)
  call spec_file_close(sub_file)

  suffix = '.nc'

  ! What is our grid - we should get this from the netcdf file
  write(filename_ic,'(A,A)') trim(file_path), "wrfinput_d01"
  call pmc_nc_check(nf90_open(filename_ic, NF90_NOWRITE, ncid_ic))
  status = nf90_get_att(ncid_ic,NF90_GLOBAL,"WEST-EAST_GRID_DIMENSION",ie)
  call pmc_nc_check(nf90_open(filename_ic, NF90_NOWRITE, ncid_ic))
  status = nf90_get_att(ncid_ic,NF90_GLOBAL,"SOUTH-NORTH_GRID_DIMENSION",je)
  call pmc_nc_check(nf90_open(filename_ic, NF90_NOWRITE, ncid_ic))
  status = nf90_get_att(ncid_ic,NF90_GLOBAL,"BOTTOM-TOP_GRID_DIMENSION",nz)

  is = 1 
  ie = ie - 1 
  nx = ie - is + 1 
  js = 1 
  je = je - 1
  ny = je - js + 1 
  nz = nz - 1 
  nt = 4 

  if (mam3) then
     n_modes = 3
  else
     n_modes = 4
  end if
  n_aero_species = aero_data_n_spec(aero_data)

  allocate(mode_diams(n_modes))
  allocate(mode_std(n_modes))
  allocate(mode_num_concs(n_modes))
  allocate(mode_vol_fracs(n_modes,n_aero_species))
  allocate(mode_source(n_modes))
  allocate(aero_spec_name(n_modes))

  ! Zero everything out
  mode_num_concs = 0.0d0
  mode_vol_fracs = 0.0d0

  ! Set values
  if (mam3) then
     aero_spec_name(1) = "aitken"
     mode_diams(1) = 0.0260d-6
     mode_std(1) =  1.6d0
     mode_source(1) = 1 

     aero_spec_name(2) = "accumulation"
     mode_diams(2) =  0.11d-6
     mode_std(2) =  1.8d0
     mode_source(2) = 2 

     aero_spec_name(3) = "coarse"
     mode_diams(3) =  2.0d-6
     mode_std(3) =  1.8d0
     mode_source(3) = 3 
  else
     aero_spec_name(1) = "so4"
     mode_diams(1) =  69.5 * 2.0d0 
     mode_std(1) =  2.03d0
     mode_vol_fracs(1,1) = 1.0
     mode_source(1) = 1

     aero_spec_name(2) = 'no3'
     mode_diams(2) = 69.5 * 2.0d0
     mode_std(2) = 2.03d0
     mode_vol_fracs(2,2) = 62.0d0 / 80.0d0 ! no3
     mode_vol_fracs(2,4) = 18.0d0 / 80.0d0 
     mode_source(2) = 2

     aero_spec_name(3) = 'bc'
     mode_diams(3) = 11.8 * 2.0d0 
     mode_std(3) = 2.0d0
     mode_vol_fracs(3,19) = 1.0
     mode_source(3) = 3

     aero_spec_name(4) = 'oc'
     mode_diams(4) = 21.2 * 2.0d0
     mode_std(4) = 2.2d0
     mode_vol_fracs(4,18) = 1.0
     mode_source(4) = 4

    ! Convert from nanometers to meters
    mode_diams = mode_diams / 1e9
  end if

  allocate(aero_values(n_modes,n_aero_species,is:ie,js:je,nz))
  allocate(density(is:ie,js:je,nz))

  call load_data_aero(aero_values, n_modes, n_aero_species, nx, &
       ny, nz, aero_spec_name, density, ncid_ic, aero_data)

  n_proc = pmc_mpi_size()
  rem = mod(nx, n_proc)
  if (pmc_mpi_rank() < rem) then
      is_local = pmc_mpi_rank() * ((nx / n_proc) + 1) + 1
      ie_local =( pmc_mpi_rank() + 1)  * ((nx / n_proc) + 1)
  else
      is_local = pmc_mpi_rank() * (nx / n_proc) + rem + 1
      ie_local =( pmc_mpi_rank() + 1)  * (nx / n_proc) + rem
  end if

  print*, 'MPI rank:', pmc_mpi_rank(), 'start: ', is_local, 'finish: ',ie_local 

  do i = is_local,ie_local
     do j = js,je
!        do k = 1,nz
           call create_ics(i, j, nz, n_modes, n_aero_species, &
                aero_values(:,:,i,j,:), aero_data, mode_diams, &
                mode_std, mode_vol_fracs, mode_source, file_prefix)
!        end do
     end do
  end do

  deallocate(mode_diams)
  deallocate(mode_std)
  deallocate(mode_num_concs)
  deallocate(mode_vol_fracs)
  deallocate(aero_spec_name)
  deallocate(aero_values)

  call pmc_mpi_finalize()
 
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subroutine to output to a file the aerosol distributions.
  subroutine aero_dist_output_netcdf(aero_dists, ncid)

    !> Array of aerosol distributions.
    type(aero_dist_t), intent(in) :: aero_dists
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

  !> Compute a number concentration based off a log-normal distribution with
  !! a given diameter and standard deviation
  real(kind=dp) function get_num_conc(mass, diam, std, species_density)

    !> Total mass.
    real(kind=dp) :: mass
    !> Median diamter
    real(kind=dp) :: diam
    !> Standard deviation.
    real(kind=dp) :: std
    !> Density of aerosol species.
    real(kind=dp) :: species_density

    real(kind=dp) :: tmp

    tmp = (species_density*const%pi/ 6.0d0)*diam**3.0 * exp(4.5d0 * log(std)**2.0d0)

    get_num_conc = mass / tmp

    ! Fraction of the number concentration that is within the MOSAIC size range
    !tmp = .5 * (erf(log(Dmax / diam) / (sqrt(2.0d0) * log (std))) - &
    !     erf(log(Dmin / diam) / (sqrt(2.0d0)*log(std))))
    !print*, diam, tmp
    !get_num_conc = get_num_conc / tmp

  end function get_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a boundary condition NetCDF for a given grid cell.
  subroutine create_ics(i, j, nz, num_modes, num_aero_species, values, &
       aero_data, mode_diams, mode_std, mode_vol_fracs, mode_source, &
       file_prefix)

    !> Index i for grid cell.
    integer, intent(in) :: i
    !> Index j for grid cell.
    integer, intent(in) :: j
    !> Number of vertical levels.
    integer, intent(in) :: nz 
    !> Number of aerosol modes.
    integer, intent(in) :: num_modes
    !> Number of aerosol species.
    integer, intent(in) :: num_aero_species
    !> Aerosol data
    type(aero_data_t) :: aero_data
    !> Median diameters of each mode.
    real(kind=dp), intent(in) :: mode_diams(num_modes)
    !> Standard deviation of each mode.
    real(kind=dp), intent(in) :: mode_std(num_modes)
    !> Volume fractions of each mode.
    real(kind=dp), intent(in) :: mode_vol_fracs(num_modes, 20)
    !> Source number.
    integer, intent(in) :: mode_source(num_modes)
    !> Aerosol mass concentrations.
    real(kind=dp), intent(in) :: values(num_modes, num_aero_species, nz)
    !> Filename prefix for output.
    character(len=300) :: file_prefix

    character(len=100) :: prefix, suffix
    character(len=300) :: filename
    type(aero_dist_t) :: aero_dists
    integer :: i_time
    integer :: start(1), count(1)
    integer :: start2(2), count2(2)
    integer :: check_dim_size
    integer :: dimid_n_times
    character(len=NF90_MAX_NAME) :: check_name
    integer :: ncid, ncid_group
    real(kind=dp) :: so4_num_conc, no3_num_conc, nh4_ratio, so4_ratio, &
         no3_ratio
    integer :: i_spec
    real(kind=dp) :: total_num_conc
    real(kind=dp), allocatable :: vol_frac(:)
    integer :: k 
    character(len=50) :: group_name

    suffix = '.nc'

    write(filename, '(a,a,i3.3,a,i3.3,a)') trim(file_prefix), '_', &
         i,'_',j,trim(suffix)
    call pmc_nc_open_write(filename, ncid)

    do k = 1,nz
       call pmc_nc_check_msg(nf90_redef(ncid),'in define mode for level')
       write(group_name,'(a,i3.3)') 'level_', k 
       call pmc_nc_check_msg(nf90_def_grp(ncid, group_name, ncid_group), &
            'creating level group')
       call pmc_nc_check_msg(nf90_enddef(ncid),'end define mode for level')

       allocate(aero_dists%mode(num_modes))
       do i_mode = 1,num_modes
          total_num_conc = 0.0d0
          aero_dists%mode(i_mode)%name = 'string'
          aero_dists%mode(i_mode)%type =  1
          aero_dists%mode(i_mode)%char_radius = &
               mode_diams(i_mode) / 2.0
          aero_dists%mode(i_mode)%log10_std_dev_radius = &
               dlog10(mode_std(i_mode))
          do i_spec = 1,num_aero_species
             total_num_conc = total_num_conc + get_num_conc( &
                  values(i_mode,i_spec,k), mode_diams(i_mode), mode_std(i_mode), &
                  aero_data%density(i_spec))
          end do
          aero_dists%mode(i_mode)%num_conc = total_num_conc
          vol_frac = values(i_mode,:,k) / aero_data%density
          aero_dists%mode(i_mode)%vol_frac = vol_frac / sum(vol_frac)
          ! FIXME: Source number
          aero_dists%mode(i_mode)%source = mode_source(i_mode)
       end do
       call aero_dist_output_netcdf(aero_dists,ncid_group)
       deallocate(aero_dists%mode)
    end do
    call pmc_nc_check(nf90_close(ncid))

  end subroutine create_ics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load aerosol data from WRF boundary condition file.
  subroutine load_data_aero(mass_conc, num_mode, num_spec, nx, ny, nz, &
       spec_name, density, ncid, aero_data)

    !> Mass concentration.
    real(kind=dp), intent(inout), dimension(num_mode,num_spec,nx,ny,nz) :: &
         mass_conc
    !> Number of modes.
    integer, intent(in) :: num_mode
    !> Number of species
    integer, intent(in) :: num_spec
    !> Number of grid cells in west-east direction.
    integer, intent(in) :: nx
    !> Number of grid cells in the south-north direction.
    integer, intent(in) :: ny
    !> Number of grid cells in the vertical direction.
    integer, intent(in) :: nz
    !> Density of air.
    real(kind=dp), intent(in), dimension(nx,ny,nz) :: density
    !> NetCDF file.
    integer, intent(in) :: ncid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
   
    character(len=100) :: var_name
    character(len=100), dimension(num_mode) :: spec_name
    integer, parameter :: num_bin = 8
    real(kind=dp), allocatable, dimension(:,:,:,:) :: temp_species
    integer :: i_mode, j_bin, i_spec
    integer :: ind
    integer :: dimid
    integer :: bdy_width
    integer :: i_time, i, spec_index, i_mode_pmc
    
    ! MAM3 species: so4
    integer, parameter :: num_mode_mam3 = 3 
    integer, parameter :: num_spec_mam3 = 6 
    logical, dimension(num_mode_mam3,num_spec_mam3) :: mode_contains_species
    character(len=3), dimension(num_spec_mam3) :: aero_names_mam3 = &
       ["so4", "ncl", "bc ", "pom", "soa", "dst"]
    character(len=4), dimension(num_spec_mam3) :: aero_names_pmc = &
       ["SO4 ","Na  ", "BC  ", "OC  ", "API1", "OIN "]
    integer, dimension(num_mode_mam3) :: pmc_to_mam3_map

    real(kind=dp), parameter :: mw_na = 22.990d0
    real(kind=dp), parameter :: mw_cl = 35.450d0

    real(kind=dp), dimension(num_spec_mam3) :: aero_factors = &
       [96.0d0/115.0d0, mw_na / (mw_na + mw_cl), 1.0d0, 1.0d0, 1.0d0, 1.0d0]


 
    mode_contains_species = .false.

    mode_contains_species(1,1) = .true.
    mode_contains_species(1,2) = .true.
    mode_contains_species(1,5) = .true.
    mode_contains_species(2,:) = .true.
    mode_contains_species(3,1) = .true.
    mode_contains_species(3,2) = .true.
    mode_contains_species(3,6) = .true.

    pmc_to_mam3_map(1) = 2
    pmc_to_mam3_map(2) = 1
    pmc_to_mam3_map(3) = 3

    if (pmc_mpi_rank() == 0) then
       print*, 'loading aerosol data'
    end if

    mass_conc = 0.0d0

    do i_mode = 1,num_mode_mam3
       i_mode_pmc = pmc_to_mam3_map(i_mode)
       do i_spec = 1,num_spec_mam3
          if (mode_contains_species(pmc_to_mam3_map(i_mode), i_spec)) then
             ! in Reverse: Time, bdy_width, bottom_top, south_north
             write(var_name,'(A,A,I1.1)') trim(aero_names_mam3(i_spec)),"_a",i_mode
             if (pmc_mpi_rank() == 0) then
                print*, 'reading aerosol species ', var_name
             end if
             call pmc_nc_read_real_4d(ncid, temp_species, var_name, .true.)
             ! Find the partmc species
             spec_index = aero_data_spec_by_name(aero_data, trim(aero_names_pmc(i_spec)))
             mass_conc(i_mode_pmc,spec_index,:,:,:) = mass_conc(i_mode_pmc,spec_index,:,:,:) &
                  + temp_species(:,:,:,1) * aero_factors(i_spec)
             if (aero_names_mam3(i_spec) == "so4") then
                spec_index = aero_data_spec_by_name(aero_data, "NH4")
                mass_conc(i_mode_pmc,spec_index,:,:,:) = &
                     mass_conc(i_mode_pmc,spec_index,:,:,:) &
                     + temp_species(:,:,:,1) * (1.0d0 - aero_factors(i_spec))
             else if (aero_names_mam3(i_spec) == "ncl") then
                spec_index = aero_data_spec_by_name(aero_data, "Cl")
                mass_conc(i_mode_pmc,spec_index,:,:,:) = &
                     mass_conc(i_mode_pmc,spec_index,:,:,:) &
                     + temp_species(:,:,:,1) * (1.0d0 - aero_factors(i_spec))
             end if
          end if
       end do
    end do

    ! Convert from WRF-Chem units of micrograms per kg of dry air to kg per m^3
    do i_mode = 1,num_mode
    do i_spec = 1,num_spec
    do i = 1,nx
    do j = 1,ny
    do k = 1,nz
        mass_conc(i_mode,i_spec,i,j,k) = mass_conc(i_mode,i_spec,i,j,k) / 1.0d9
    end do
    end do
    end do
    end do
    end do

  end subroutine load_data_aero

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
        & <outfile prefix>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program make_ics
