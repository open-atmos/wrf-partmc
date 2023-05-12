program compute_emission_statistics

  use netcdf
  use pmc_netcdf
  use pmc_aero_dist
  use pmc_aero_mode
  use pmc_aero_data
  use pmc_spec_file
  use pmc_mpi
  use pmc_stats
  use mpi
  use pmc_scenario
  use pmc_aero_state

  implicit none

  !> Filename of NetCDF file to open.
  character(len=200) :: filename_grid
  character(len=200) :: filename_emissions
  !> NetCDF file ID, in data mode.
  real(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: aero_emissions, &
       aero_emissions_patch
  character(len=100) :: name
  integer :: varid, dimid_times, send_size_a, send_size_g, ierr
  integer :: varid_aero_source
  integer :: status
  integer :: mpi_status(MPI_STATUS_SIZE)
  integer :: var_size(4)
  integer :: temp_size(2)
  integer :: dim
  integer :: t
  integer :: n_source_classes
  integer :: check_dim_size, dimid_n_modes, n_modes, dimid_aero_sources
  character(len=NF90_MAX_NAME) :: check_name
  real(kind=dp) :: total_mass
  real(kind=dp) :: diam_i, diam_j
  real(kind=dp) :: std_i, std_j
  real(kind=dp) :: num_conc_i, num_conc_j
  integer :: x, y, z, nx, ny, nz, nt, s
  integer :: i,j,k, i_source, n_species, i_mode, i_type
  integer :: i_rank, buffer_size, position
  type(aero_data_t) :: aero_data
  character(len=9), allocatable :: partmc_species(:)
  character(len=9), allocatable :: smoke_species(:)
  character(len=9), allocatable :: wrf_species(:)

  ! output structures
  type(aero_dist_t), allocatable :: aero_emission(:)
  real(kind=dp), allocatable, dimension(:) :: aero_emission_rate, &
       aero_emission_time
  integer :: n_mode_counter
  real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std, &
       mode_num_concs
  integer, allocatable, dimension(:) :: mode_source
  integer, allocatable, dimension(:) :: mode_weight_class
  real(kind=dp), allocatable, dimension(:,:) :: mode_vol_fracs, &
       mode_vol_frac_std
  real(kind=dp), allocatable, dimension(:) :: vol_frac
  character(len=100) :: prefix, filename, case_name, dir_path
  character(len=3), parameter  :: suffix = '.nc'
  integer :: ncid
  character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
  type(spec_file_t) :: sub_file

  character(len=300) :: spec_name
  character(len=300) :: grid_name
  character(len=100), allocatable, dimension(:) :: mode_names
  character(len=100) :: name_i, name_j
  integer :: start_date, n_days_per_file, n_days, i_day, i_date
  real(kind=dp) :: gas_scale_factor, aerosol_scale_factor
  real(kind=dp) :: dx

  integer :: i_send_s, i_send_e
  integer :: rem, is_local, ie_local, n_proc, i_proc, root
  integer :: recv_size_a, recv_size_g
  integer, allocatable, dimension(:) :: displs_a, send_counts_a
  integer, allocatable, dimension(:) :: displs_g, send_counts_g
  integer, allocatable, dimension(:) :: aero_source_centers
  integer :: stride, block_length
  integer :: grid_emission
  integer :: n_smoke_species
  integer :: n_emit
  real(kind=dp) :: delta_t, p
  character(len=(AERO_SOURCE_NAME_LEN)*100) :: aero_source_names
  character(len=100) :: unit
  type(scenario_t) :: scenario
  type(aero_state_t) :: aero_state
  type(aero_dist_t) :: emissions, aero_dist_init
  real(kind=dp) :: emission_rate_scale
    character(len=AERO_SOURCE_NAME_LEN) :: weight_class
    character(len=AERO_SOURCE_NAME_LEN), allocatable :: weight_name(:)
    integer :: num_ic_modes, num_bc_modes, num_emit_modes
    integer, allocatable, dimension(:) :: emission_weight_classes
    character(len=200) :: file
  integer :: i_emit, dummy
  type(bin_grid_t) :: bin_grid
  real(kind=dp), allocatable :: chi_values(:,:,:), chi_ccn_values(:,:,:)
  real(kind=dp), allocatable :: d_gamma_ccn(:,:,:), d_gamma(:,:,:)
  real(kind=dp), allocatable :: d_alpha_ccn(:,:,:), d_alpha(:,:,:)
  character, allocatable :: buffer(:)
  real(kind=dp), allocatable :: temp_array(:,:,:)
  integer :: is_remote, ie_remote
  character, allocatable :: buffer_chi_ccn(:), buffer_chi(:)
  character, allocatable :: buffer_d_alpha(:), buffer_d_alpha_ccn(:)
  character, allocatable :: buffer_d_gamma(:), buffer_d_gamma_ccn(:)
  call pmc_mpi_init()

  ! We need aero_data
  sub_filename = 'aero_data.dat'
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_data(sub_file, aero_data)
  call spec_file_close(sub_file)
  call fractal_set_spherical(aero_data%fractal)
  ! FIXME: do not hardcode this
  nx = 169
  ny = 159
  nt = 25
  allocate(d_alpha(nx,ny,nt))
  allocate(d_gamma(nx,ny,nt))
  allocate(chi_values(nx,ny,nt))
  allocate(d_alpha_ccn(nx,ny,nt))
  allocate(d_gamma_ccn(nx,ny,nt))
  allocate(chi_ccn_values(nx,ny,nt))
  delta_t = 60.0d0
  emission_rate_scale = 1.0d0
  prefix = "cares_emissions_new/aero_emit_dists"
  i = 1
  j = 1
  k = 1
  ! Need to know how many weight classes
  write(file, '(a,a,i3.3,a,i3.3,a,i3.3,a)') &
         trim(prefix),'_',i,'_',j,'_',k,'.nc'
  call pmc_nc_open_read(file, ncid)
  status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
  call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, &
         check_name, check_dim_size))
  num_emit_modes = check_dim_size
  call pmc_nc_read_integer_1d(ncid, emission_weight_classes, &
       "source_weight_class", .true.)
  call pmc_nc_check(nf90_close(ncid))
  call ensure_string_array_size(weight_name, num_emit_modes)
  do i_emit = 1,num_emit_modes
     write(weight_class,'(I3)') num_emit_modes !emission_weight_classes(i_emit)
     weight_name(i_emit) = trim(weight_class)
  end do
  do i_emit = 1,num_emit_modes
     dummy = aero_data_weight_class_by_name(aero_data, weight_name(i_emit))
  end do

  call aero_state_zero(aero_state)
  call aero_state_set_weight(aero_state, aero_data, &
       AERO_STATE_WEIGHT_FLAT_SPECIFIED)
  call aero_state_set_n_part_ideal(aero_state, 10000.0d0)
  call aero_state_sort(aero_state, aero_data)

  n_proc = pmc_mpi_size()
  rem = mod(nx, n_proc)
  if (pmc_mpi_rank() == 0) then
     is_local = pmc_mpi_rank() * (nx / n_proc) + 1
  else
     is_local = pmc_mpi_rank() * (nx / n_proc) + rem + 1
  end if
  ie_local =( pmc_mpi_rank() + 1)  * (nx / n_proc) + rem

  print*, pmc_mpi_rank(), is_local, ie_local

  do i = is_local,ie_local
     do j = 1,ny
        do k = 1,1
        ! Load the data
           call read_in_emissions(scenario, i, j, k, aero_data, prefix)
           do t = 1,nt
           emissions = scenario%aero_emission(t)
!           call aero_dist_interp_1d(scenario%aero_emission, &
!                scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
!                (t-1)*3600.0d0, emissions, emission_rate_scale)
           p = 1.0d0 * delta_t
           call aero_state_add_aero_dist_sample(aero_state, aero_data, &
                emissions, p, 1.0d0, 0.0d0, .true., .true., &
                n_emit)
           if (aero_state_total_particles(aero_state) > 0) then
           call aero_state_mixing_state_metrics(aero_state, &
                aero_data, d_alpha_ccn(i,j,t), d_gamma_ccn(i,j,t), &
                chi_ccn_values(i,j,t),  &
                exclude=["H2O"], group=["BC ", "OC ","OIN"])
           call aero_state_mixing_state_metrics(aero_state, &
                aero_data, d_alpha(i,j,t), d_gamma(i,j,t), chi_values(i,j,t), &
                exclude=["H2O"])
           else
              d_alpha(i,j,t) = -1.0d0
              d_gamma(i,j,t) = -1.0d0
              chi_values(i,j,t) = -1.0d0
              d_alpha_ccn(i,j,t) = -1.0d0
              d_gamma_ccn(i,j,t) = -1.0d0
              chi_ccn_values(i,j,t) = -1.0d0
           end if
           call aero_state_zero(aero_state)
           end do
        end do
     end do
  end do

  call pmc_mpi_barrier()
  ! Pack all our local arrays
  if (pmc_mpi_rank() == 0) then
     do i_rank = 1, pmc_mpi_size() - 1
        is_remote = i_rank * (nx / n_proc) + rem + 1
        ie_remote = (i_rank + 1)  * (nx / n_proc) + rem
        allocate(temp_array(is_remote:ie_remote,ny,nt))
        buffer_size =  pmc_mpi_pack_size_real_array_3d(temp_array)
        allocate(buffer(buffer_size))
        call MPI_recv(buffer, buffer_size, MPI_PACKED, i_rank, 0, &
             MPI_COMM_WORLD, mpi_status, ierr)
        position = 0
        call pmc_mpi_unpack_real_array_3d(buffer, position, temp_array)
        chi_values(is_remote:ie_remote,:,:) = temp_array
        call MPI_recv(buffer, buffer_size, MPI_PACKED, i_rank, 0, &
             MPI_COMM_WORLD, mpi_status, ierr)
        position = 0
        call pmc_mpi_unpack_real_array_3d(buffer, position, temp_array)
        d_gamma(is_remote:ie_remote,:,:) = temp_array
        call MPI_recv(buffer, buffer_size, MPI_PACKED, i_rank, 0, &
             MPI_COMM_WORLD, mpi_status, ierr)
        position = 0
        call pmc_mpi_unpack_real_array_3d(buffer, position, temp_array)
        d_alpha(is_remote:ie_remote,:,:) = temp_array
        position = 0
        call pmc_mpi_unpack_real_array_3d(buffer, position, temp_array)
        chi_ccn_values(is_remote:ie_remote,:,:) = temp_array
        call MPI_recv(buffer, buffer_size, MPI_PACKED, i_rank, 0, & 
             MPI_COMM_WORLD, mpi_status, ierr)
        position = 0
        call pmc_mpi_unpack_real_array_3d(buffer, position, temp_array)
        d_gamma_ccn(is_remote:ie_remote,:,:) = temp_array
        call MPI_recv(buffer, buffer_size, MPI_PACKED, i_rank, 0, &
             MPI_COMM_WORLD, mpi_status, ierr)
        position = 0
        call pmc_mpi_unpack_real_array_3d(buffer, position, temp_array)
        d_alpha_ccn(is_remote:ie_remote,:,:) = temp_array

        deallocate(temp_array)
     end do
  else
     temp_array = chi_values(is_local:ie_local,:,:)
     buffer_size = pmc_mpi_pack_size_real_array_3d(temp_array)
     allocate(buffer_chi(buffer_size))
     allocate(buffer_d_alpha(buffer_size))
     allocate(buffer_d_gamma(buffer_size))
     allocate(buffer_chi_ccn(buffer_size))
     allocate(buffer_d_alpha_ccn(buffer_size))
     allocate(buffer_d_gamma_ccn(buffer_size))

     position = 0
     call pmc_mpi_pack_real_array_3d(buffer_chi, position, &
          temp_array)
     position = 0
     temp_array = d_gamma(is_local:ie_local,:,:)
     call pmc_mpi_pack_real_array_3d(buffer_d_gamma, position, &
          temp_array)
     position = 0
     temp_array = d_alpha(is_local:ie_local,:,:)
     call pmc_mpi_pack_real_array_3d(buffer_d_alpha, position, &
          temp_array)

     position = 0
     temp_array =  chi_ccn_values(is_local:ie_local,:,:)
     call pmc_mpi_pack_real_array_3d(buffer_chi_ccn, position, &
          temp_array)
     position = 0
     temp_array = d_gamma_ccn(is_local:ie_local,:,:)
     call pmc_mpi_pack_real_array_3d(buffer_d_gamma_ccn, position, &
          temp_array)
     position = 0
     temp_array = d_alpha_ccn(is_local:ie_local,:,:)
     call pmc_mpi_pack_real_array_3d(buffer_d_alpha_ccn, position, &
          temp_array)

     call MPI_send(buffer_chi, buffer_size, MPI_PACKED, 0, pmc_mpi_rank(), &
             MPI_COMM_WORLD, ierr)
     call MPI_send(buffer_d_gamma, buffer_size, MPI_PACKED, 0, pmc_mpi_rank(), &
             MPI_COMM_WORLD, ierr)
     call MPI_send(buffer_d_alpha, buffer_size, MPI_PACKED, 0, pmc_mpi_rank(), &
             MPI_COMM_WORLD, ierr)
     call MPI_send(buffer_chi_ccn, buffer_size, MPI_PACKED, 0, pmc_mpi_rank(), &
             MPI_COMM_WORLD, ierr)
     call MPI_send(buffer_d_gamma_ccn, buffer_size, MPI_PACKED, 0, pmc_mpi_rank(), &
             MPI_COMM_WORLD, ierr)
     call MPI_send(buffer_d_alpha_ccn, buffer_size, MPI_PACKED, 0, pmc_mpi_rank(), &
             MPI_COMM_WORLD, ierr)
  end if

  call pmc_mpi_barrier()

  if (pmc_mpi_rank() == 0) then
  call pmc_nc_open_write("chi_values.nc", ncid)
  call pmc_nc_write_real_3d(ncid, d_alpha, "d_alpha", dim_name_1="nx", &
       dim_name_2="ny", dim_name_3="nt")
  call pmc_nc_write_real_3d(ncid, d_gamma, "d_gamma", dim_name_1="nx", &
       dim_name_2="ny", dim_name_3="nt")
  call pmc_nc_write_real_3d(ncid, chi_values, "chi", dim_name_1="nx", &
       dim_name_2="ny", dim_name_3="nt")
  call pmc_nc_write_real_3d(ncid, d_alpha_ccn, "d_alpha_ccn", dim_name_1="nx", &
       dim_name_2="ny", dim_name_3="nt")
  call pmc_nc_write_real_3d(ncid, d_gamma_ccn, "d_gamma_ccn", dim_name_1="nx", &
       dim_name_2="ny", dim_name_3="nt")
  call pmc_nc_write_real_3d(ncid, chi_ccn_values, "chi_ccn", dim_name_1="nx", &
       dim_name_2="ny", dim_name_3="nt")

  call pmc_nc_check(nf90_close(ncid)) 

  end if

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read in gas and aerosol emissions for a grid cell.
  subroutine read_in_emissions(scenario, i, j, k, aero_data, prefix)

    !> Scenario data.
    type(scenario_t),intent(inout) :: scenario
    !> East-west index of grid cell.
    integer,intent(in) :: i
    !> North-south index of grid cell.
    integer,intent(in) :: j
    !> Top-bottom index of grid cell.
    integer,intent(in) :: k
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !>
    character(len=*), intent(in) :: prefix

    character(len=200) :: file
    character(len=AERO_MODE_NAME_LEN) :: mode_name, weight_class
    integer :: n_time, ncid
    character(len=NF90_MAX_NAME) :: check_name
    integer :: status, check_dim_size
    integer :: dimid_n_times, dimid_n_modes, dimid_n_aero_specs
    integer :: i_loop
    integer :: n_modes, n_aero_specs, n_gas_specs

    real(kind=dp), allocatable, dimension(:) :: char_radius, std
    integer, allocatable, dimension(:) :: source
    real(kind=dp), allocatable, dimension(:,:) :: num_conc
    real(kind=dp), allocatable, dimension(:,:,:) :: vol_frac, vol_frac_std
    integer :: i_time, i_mode
    integer :: dummy
    ! Source information
    character(len=1000) :: name
    integer, parameter :: MAX_SOURCES = 1000
    integer :: dimid_aero_source, n_source, varid_aero_source, i_source, i_str
    character(len=((AERO_SOURCE_NAME_LEN + 2) * MAX_SOURCES)) &
         :: aero_source_names
    character(len=AERO_SOURCE_NAME_LEN), allocatable :: source_name(:)
    integer, allocatable, dimension(:) :: emission_weight_classes

    write(file, '(a,a,i3.3,a,i3.3,a,i3.3,a)') &
         trim(prefix),'_',i,'_',j,'_',k,'.nc'

    call pmc_nc_open_read(file, ncid)

    status = nf90_inq_dimid(ncid, "n_times", dimid_n_times)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_times, check_name, &
            check_dim_size))
    n_time = check_dim_size

    if (allocated(scenario%aero_emission)) deallocate(scenario%aero_emission)
    if (allocated(scenario%aero_emission_time)) &
         deallocate(scenario%aero_emission_time)
    if (allocated(scenario%aero_emission_rate_scale)) &
         deallocate(scenario%aero_emission_rate_scale)

    allocate(scenario%aero_emission_time(n_time))
    allocate(scenario%aero_emission_rate_scale(n_time))
    allocate(scenario%aero_emission(n_time))

    call pmc_nc_read_real_1d(ncid, scenario%aero_emission_time, &
         'aero_emission_time', .true.)

    status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, check_name, &
            check_dim_size))
    n_modes = check_dim_size
    status = nf90_inq_dimid(ncid, "n_aero_specs", dimid_n_aero_specs)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_aero_specs, &
         check_name, check_dim_size))
    n_aero_specs = check_dim_size

    allocate(char_radius(n_modes))
    call pmc_nc_read_real_1d(ncid, char_radius, 'char_radius', .true.)
    allocate(std(n_modes))
    call pmc_nc_read_real_1d(ncid, std, 'log10_std_dev_radius', .true.)

    allocate(source(n_modes))
    call pmc_nc_read_integer_1d(ncid, source, 'source_id', .true.)
    allocate(num_conc(n_modes, n_time))
    call pmc_nc_read_real_2d(ncid, num_conc, 'num_conc', .true.)
    allocate(vol_frac(n_aero_specs,n_modes,n_time))
    call pmc_nc_read_real_3d(ncid, vol_frac, 'vol_frac', .true.)
    allocate(vol_frac_std(n_aero_specs,n_modes,n_time))
    vol_frac_std = 0.0d0

    ! Read in mode names
    call pmc_nc_check(nf90_inq_varid(ncid, "aero_source", &
         varid_aero_source))
    call pmc_nc_check(nf90_get_att(ncid, varid_aero_source, "names", &
         aero_source_names))
    status = nf90_inq_dimid(ncid, "n_aero_sources", dimid_aero_source)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_source, name, n_source))

    call ensure_string_array_size(source_name, n_source)
    do i_source = 1,n_source
       i_str = 1
       do while ((aero_source_names(i_str:i_str) /= " ") &
            .and. (aero_source_names(i_str:i_str) /= ","))
          i_str = i_str + 1
       end do
       call assert(840982478, i_str > 1)
       source_name(i_source) = aero_source_names(1:(i_str-1))
       aero_source_names = aero_source_names((i_str+1):)
    end do

    call pmc_nc_read_integer_1d(ncid, emission_weight_classes, &
         "source_weight_class", .true.)

    do i_time = 1, n_time
    scenario%aero_emission_rate_scale(i_time) = 1.0d0
    if (allocated(scenario%aero_emission(i_time)%mode)) &
         deallocate(scenario%aero_emission(i_time)%mode)
    allocate(scenario%aero_emission(i_time)%mode(n_modes))

    do i_mode = 1, n_modes
       scenario%aero_emission(i_time)%mode(i_mode)%name =  &
            trim(source_name(source(i_mode)))
       dummy = aero_data_source_by_name(aero_data, &
            source_name(source(i_mode)))
       scenario%aero_emission(i_time)%mode(i_mode)%type = &
            AERO_MODE_TYPE_LOG_NORMAL
       scenario%aero_emission(i_time)%mode(i_mode)%char_radius = &
            char_radius(i_mode)
       scenario%aero_emission(i_time)%mode(i_mode)%log10_std_dev_radius = &
            std(i_mode)
       scenario%aero_emission(i_time)%mode(i_mode)%num_conc = &
            num_conc(i_mode, i_time)
       scenario%aero_emission(i_time)%mode(i_mode)%vol_frac = &
            vol_frac(1:20,i_mode,i_time)
       scenario%aero_emission(i_time)%mode(i_mode)%vol_frac_std = &
            vol_frac_std(1:20,i_mode,i_time)
       ! FIXME: Fix this source information
       scenario%aero_emission(i_time)%mode(i_mode)%source = dummy
       !
       write(weight_class,'(I3)') i_mode !emission_weight_classes(i_mode)
       scenario%aero_emission(i_time)%mode(i_mode)%weight_class = &
            aero_data_weight_class_by_name(aero_data, weight_class)
    end do
    end do

    call pmc_nc_check(nf90_close(ncid))

    deallocate(char_radius)
    deallocate(std)
    deallocate(source)
    deallocate(num_conc)
    deallocate(vol_frac)
    deallocate(vol_frac_std)

  end subroutine read_in_emissions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program compute_emission_statistics 
