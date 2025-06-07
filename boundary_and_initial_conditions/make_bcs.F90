program make_bcs

  use netcdf
  use pmc_netcdf
  use pmc_aero_dist
  use pmc_aero_mode
  use pmc_aero_data
  use pmc_spec_file
  use pmc_mpi
  use mpi
  use ic_bc_helper

  implicit none
 
  !> Filename of NetCDF file to open.
  character(len=100) :: filename_bc, filename_ic
  !> NetCDF file ID, in data mode.
  integer :: ncid_bc, ncid_ic

  character(len=100), allocatable, dimension(:) :: aero_mode_name
  real(kind=dp), allocatable, dimension(:,:,:,:,:) :: mass_conc
  real(kind=dp), allocatable, dimension(:,:,:) :: density
  character(len=100) :: name
  integer :: varid, dimid_times
  integer :: status
  integer :: dim
  integer :: n_modes
  integer :: nx, ny, nz, nt
  integer :: i,j,k
  integer :: i_mode
  integer :: n_spec
  type(aero_data_t) :: aero_data
  integer :: is, ie, js, je
  integer :: is_local, ie_local, js_local, je_local
  real(kind=dp), allocatable, dimension(:) :: aero_bc_rate, &
       aero_bc_time
  real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std
  integer, allocatable, dimension(:) :: mode_source
  real(kind=dp), allocatable, dimension(:,:) :: mode_vol_fracs
  character(len=100) :: suffix
  integer :: ncid
  type(spec_file_t) :: sub_file
  character(len=300) :: file_path, filename, file_prefix, command
  character(len=300) :: output_file_path, output_file_prefix
  integer :: ks, ke, n_proc, rem_x, rem_y

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
  write(filename_bc,'(A,A)') trim(file_path),"wrfbdy_d01"
  call pmc_nc_check(nf90_open(filename_bc, NF90_NOWRITE, ncid_bc))
  ! What is our grid - we should get this from the netcdf file
  write(filename_ic,'(A,A)') trim(file_path),"wrfbdy_d01"
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
  nt = 8
  ks = 1
  ke = nz - 1
  n_proc = pmc_mpi_size()
  rem_x = mod(nx, n_proc)
  rem_y = mod(ny, n_proc)

  if (pmc_mpi_rank() < rem_x) then
      is_local = pmc_mpi_rank() * ((nx / n_proc) + 1) + 1
      ie_local =( pmc_mpi_rank() + 1)  * ((nx / n_proc) + 1)
  else
      is_local = pmc_mpi_rank() * (nx / n_proc) + rem_x + 1
      ie_local =( pmc_mpi_rank() + 1)  * (nx / n_proc) + rem_x
  end if

  if (pmc_mpi_rank() < rem_y) then
      js_local = pmc_mpi_rank() * ((ny / n_proc) + 1) + 1
      je_local =( pmc_mpi_rank() + 1)  * ((ny / n_proc) + 1)
  else
      js_local = pmc_mpi_rank() * (ny / n_proc) + rem_y + 1
      je_local =( pmc_mpi_rank() + 1)  * (ny / n_proc) + rem_y
  end if

!  print*, 'MPI rank: ', pmc_mpi_rank(), 'start:', , 'finish:', ke

  allocate(aero_bc_time(nt))
  allocate(aero_bc_rate(nt))

  ! We can actually load this from MOZART
  do i = 1,nt
     aero_bc_time(i) = 6*(i-1)*3600.0d0
     aero_bc_rate(i) = 1.0d0
  end do

  n_modes = 3

  n_spec = aero_data_n_spec(aero_data)

  call get_mode_parameters(aero_mode_name, mode_diams, mode_std, &
       mode_vol_fracs, mode_source, aero_data)

  ! do XS
  i = is
  allocate(mass_conc(n_modes,n_spec,js:je,nz,nt))

  call pmc_mpi_barrier()

  call load_data_aero(mass_conc, 'BXS', n_modes, n_spec, nt, &
       ny, nz, ncid_bc)

  do j = js_local,je_local
     call create_bcs(i, j, nz, aero_data, n_modes, mass_conc(:,:,j,:,:), &
          n_spec, nt, mode_diams, &
          mode_std, mode_vol_fracs, mode_source, aero_bc_rate, aero_bc_time, &
          file_prefix)
  end do

  ! do XE
  i = ie
 
  call load_data_aero(mass_conc, 'BXE', n_modes, n_spec, nt, &
       ny, nz, ncid_bc)

  do j = js_local,je_local
     call create_bcs(i, j, nz, aero_data, n_modes, mass_conc(:,:,j,:,:), &
          n_spec, nt, mode_diams, &
          mode_std, mode_vol_fracs, mode_source, aero_bc_rate, aero_bc_time, &
          file_prefix)
  end do

  deallocate(mass_conc)

  ! do YS
  j = js

  allocate(mass_conc(n_modes,n_spec,is:ie,nz,nt))

  call load_data_aero(mass_conc, 'BYS', n_modes, n_spec, nt, &
       nx, nz, ncid_bc)

  do i = is_local,ie_local
     call create_bcs(i, j, nz, aero_data, n_modes, mass_conc(:,:,i,:,:), &
          n_spec, nt, mode_diams, &
          mode_std, mode_vol_fracs, mode_source, aero_bc_rate, aero_bc_time, &
          file_prefix)
  end do

  ! do YE
  j = je

  call load_data_aero(mass_conc, 'BYE', n_modes, n_spec, nt, &
       nx, nz, ncid_bc)
  do i = is_local,ie_local 
     call create_bcs(i, j, nz, aero_data, n_modes, mass_conc(:,:,i,:,:), &
          n_spec, nt, mode_diams, &
          mode_std, mode_vol_fracs, mode_source, aero_bc_rate, aero_bc_time, &
          file_prefix)
  end do

  deallocate(mass_conc)

  call pmc_mpi_finalize()
 
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a number concentration based off a log-normal distribution with
  !! a given diameter and standard deviation
  real(kind=dp) function get_num_conc(mass, diam, std, density)

    !> Total mass.
    real(kind=dp) :: mass
    !> Median diamter
    real(kind=dp) :: diam
    !> Standard deviation.
    real(kind=dp) :: std
    !>
    real(kind=dp) :: density

    integer :: k
    real(kind=dp) :: tmp
    real(kind=dp) :: Dmax, Dmin

    k = 3
    tmp = density*(const%pi/ 6.0d0)*diam**3.0 * exp(k**2.0d0 / 2.0d0 &
         * log(std)**2.0d0)
    get_num_conc = mass / tmp

  end function get_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a boundary condition NetCDF for a given grid cell.
  subroutine create_bcs(i, j, nz, aero_data, num_modes, values, &
       num_aero_species, num_times, mode_diams, mode_std,  mode_vol_fracs, &
       mode_source, aero_bc_rate, aero_bc_time, file_prefix)

    !> Index i for grid cell.
    integer, intent(in) :: i
    !> Index j for grid cell.
    integer, intent(in) :: j
    !> Number of vertical layers.
    integer, intent(in) :: nz 
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number of aerosol modes.
    integer, intent(in) :: num_modes
    !> Number of timesteps.
    integer, intent(in) :: num_times
    !> Median diameters of each mode.
    real(kind=dp), intent(in) :: mode_diams(num_modes)
    !> Standard deviation of each mode.
    real(kind=dp), intent(in) :: mode_std(num_modes)
    !> Volume fractions of each mode.
    real(kind=dp), intent(in) :: mode_vol_fracs(num_modes, &
         aero_data_n_spec(aero_data))
    !> Source number.
    integer, intent(in) :: mode_source(num_modes)
    !> Number concentrations.
    real(kind=dp), intent(in) :: values(num_modes, num_aero_species, nz, &
         num_times)
    !> Boundary condition scale factor for aerosols.
    real(kind=dp), intent(in) :: aero_bc_rate(num_times)
    !> Boundary condition update time for aerosols.
    real(kind=dp), intent(in) :: aero_bc_time(num_times)
    !> Boundary condition filename input prefix.
    character(len=300), intent(in) :: file_prefix

    character(len=100) :: suffix, filename
    type(aero_dist_t), allocatable :: aero_dists(:)
    integer :: i_time, k
    integer :: start2(2), count2(2)
    integer :: check_dim_size
    integer :: dimid_n_times, dimid_nz, dimid_n_specs, dimid_n_modes
    character(len=NF90_MAX_NAME) :: check_name
    integer :: ncid
    real(kind=dp) :: so4_num_conc, no3_num_conc, nh4_ratio, so4_ratio, &
         no3_ratio
    real(kind=dp) :: total_num_conc
    integer :: i_spec, num_aero_species
    real(kind=dp), allocatable :: char_radius(:,:,:), vol_frac(:,:,:,:), &
         vol_frac_std(:,:,:,:), log10_std_dev_radius(:,:,:), num_conc(:,:,:)
    integer, allocatable :: mode_type(:,:,:), source(:,:,:)
    character(len=AERO_MODE_NAME_LEN) :: mode_name
    integer :: n_spec

    n_spec = aero_data_n_spec(aero_data)

    allocate(mode_type(num_times,nz,num_modes))
    allocate(char_radius(num_times,nz,num_modes))
    allocate(log10_std_dev_radius(num_times,nz,num_modes))
    allocate(num_conc(num_times,nz,num_modes))
    allocate(vol_frac(num_times,nz,num_modes,n_spec))
    allocate(vol_frac_std(num_times,nz,num_modes,n_spec))
    allocate(source(num_times,nz,num_modes))

    do k = 1,nz
    do i_time = 1,num_times
    do i_mode = 1,num_modes
       mode_name = 'string'
       mode_type(i_time,k,i_mode) = 1
       char_radius(i_time,k,i_mode) = mode_diams(i_mode) / 2.0
       log10_std_dev_radius(i_time,k,i_mode) = dlog10(mode_std(i_mode))
       total_num_conc = 0.0d0
       do i_spec = 1,num_aero_species
         total_num_conc = total_num_conc + get_num_conc( &
            values(i_mode,i_spec,k,i_time), mode_diams(i_mode), &
            mode_std(i_mode), aero_data%density(i_spec))
       end do
       num_conc(i_time,k,i_mode) = total_num_conc
       vol_frac(i_time,k,i_mode,:) = values(i_mode,:,k,i_time) / &
            sum(values(i_mode,:,k,i_time))
       vol_frac_std(i_time,k,i_mode,:) = 0.0d0
       source(i_time,k,i_mode) = mode_source(i_mode)
    end do
    end do
    end do

    suffix = '.nc'
    ! Output to the NetCDF file
    write(filename, '(a,a,i3.3,a,i3.3,a)') &
         trim(file_prefix),'_',i,'_',j,trim(suffix)
    call pmc_nc_open_write(filename, ncid)
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "n_times", &
         nt, dimid_n_times))
    call pmc_nc_check(nf90_def_dim(ncid, "n_aero_modes", n_modes, &
         dimid_n_modes))
    call pmc_nc_check(nf90_def_dim(ncid, "n_aero_specs", n_spec, &
         dimid_n_specs))
    call pmc_nc_check(nf90_def_dim(ncid, "nz", nz, dimid_nz))
    call pmc_nc_check(nf90_enddef(ncid))

    call pmc_nc_write_real_1d(ncid, aero_bc_rate, &
          "aero_bc_rate_scale", (/ dimid_n_times /), &
          unit="(1)", &
          long_name="Aerosol bc rate", &
          description="Aerosol boundary condition rate scales at set-points")
    call pmc_nc_write_real_1d(ncid, aero_bc_time, &
          "aero_bc_time", (/ dimid_n_times /), &
          unit="s", &
          long_name="Aerosol boundary condition update time", &
          description="Aerosol boundary condition set-points times (s).")
    ! Output the boundary condition data
    call pmc_nc_write_integer_3d(ncid, mode_type, "mode_type", &
         (/ dimid_n_times, dimid_nz, dimid_n_modes /))
    call pmc_nc_write_real_3d(ncid, char_radius, "char_radius", &
         (/ dimid_n_times, dimid_nz, dimid_n_modes /))
    call pmc_nc_write_real_3d(ncid, log10_std_dev_radius, &
         "log10_std_dev_radius", (/ dimid_n_times, dimid_nz, dimid_n_modes /))
    call pmc_nc_write_real_3d(ncid, num_conc, "num_conc", &
         (/ dimid_n_times, dimid_nz, dimid_n_modes /))
    call pmc_nc_write_real_4d(ncid, vol_frac, "vol_frac", &
         (/ dimid_n_times, dimid_nz, dimid_n_modes, dimid_n_specs /))
    call pmc_nc_write_real_4d(ncid, vol_frac_std, "vol_frac_std", &
         (/ dimid_n_times, dimid_nz, dimid_n_modes, dimid_n_specs /))
    call pmc_nc_write_integer_3d(ncid, source, "source", &
         (/ dimid_n_times, dimid_nz, dimid_n_modes /))

    call pmc_nc_check(nf90_close(ncid))

  end subroutine create_bcs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load aerosol data from WRF boundary condition file.
  subroutine load_data_aero(mass_conc,boundary, num_mode, num_spec, num_time, &
       n_dim1, nz, ncid)

    !> Mass concentration.
    real(kind=dp), intent(out), dimension( &
         num_mode,num_spec,n_dim1,nz,num_time) :: mass_conc
    !> Boundary of interest.
    character(len=3), intent(in) :: boundary
    !> Number of modes.
    integer :: num_mode
    !>
    integer :: num_spec
    !> Number of times.
    integer :: num_time
    !> Horizontal size of the boundary edge.
    integer :: n_dim1
    !> Vertical size of the boundary edge.
    integer :: nz
    !> NetCDF file.
    integer, intent(in) :: ncid
   
    character(len=100) :: var_name
    real(kind=dp), allocatable, dimension(:,:,:,:) :: temp_species
    integer :: i_mode, j_bin
    integer :: ind
    integer :: dimid
    integer :: bdy_width
    integer :: i_time, i
    integer :: i_spec, i_mode_pmc, spec_index

    if (pmc_mpi_rank() == 0) then
       print*, 'loading aerosol data for boundary ', boundary
    end if

    mass_conc = 0.0d0

    if (boundary(3:) .eq. "S") then
       status = NF90_INQ_DIMID(ncid, 'bdy_width', dimid)
       status = nf90_inquire_dimension(ncid, dimid, len = bdy_width)
       ind = 1
    else
       ind = 1
    end if

    do i_mode = 1,num_mode_mam3
       i_mode_pmc = pmc_to_mam3_map(i_mode)
       do i_spec = 1,num_spec_mam3
          if (mode_contains_species(i_mode, i_spec)) then
             ! in Reverse: Time, bdy_width, bottom_top, south_north
             write(var_name,'(A,A,I1.1,A,A)') trim(aero_names_mam3(i_spec)), &
                  "_a",i_mode,"_",boundary

             if (pmc_mpi_rank() == 0) then
                print*, 'reading aerosol species ', trim(var_name)
             endif
             call pmc_nc_read_real_4d(ncid, temp_species, var_name, .true.)
             ! Find the partmc species
             spec_index = aero_data_spec_by_name(aero_data, &
                  trim(aero_names_pmc(i_spec)))
             mass_conc(i_mode_pmc,spec_index,:,:,:) = &
                  mass_conc(i_mode_pmc,spec_index,:,:,:) &
                  + temp_species(:,:,ind,:) * aero_factors(i_spec)
             if (aero_names_mam3(i_spec) == "so4") then
                spec_index = aero_data_spec_by_name(aero_data, "NH4")
                mass_conc(i_mode_pmc,spec_index,:,:,:) = &
                     mass_conc(i_mode_pmc,spec_index,:,:,:) &
                     + temp_species(:,:,ind,:) * (1.0d0 - aero_factors(i_spec))
             else if (aero_names_mam3(i_spec) == "ncl") then
                spec_index = aero_data_spec_by_name(aero_data, "Cl")
                mass_conc(i_mode_pmc,spec_index,:,:,:) = &
                     mass_conc(i_mode_pmc,spec_index,:,:,:) &
                     + temp_species(:,:,ind,:) * (1.0d0 - aero_factors(i_spec))
             end if
          end if
       end do
    end do

    ! Convert from WRF-Chem units of micrograms per kg of dry air to kg per kg
    ! of dry air. WRF-PartMC will convert to kg per m-3.
    mass_conc =  mass_conc / 1d9

  end subroutine load_data_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the usage text to stderr.
  subroutine print_usage()

    write(*,*) 'Usage: ./create_bcs <scenario file path> <output file path> &
         & <output file prefix>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program make_bcs
