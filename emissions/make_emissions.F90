program make_emissions

  use netcdf
  use pmc_netcdf
  use pmc_aero_dist
  use pmc_aero_mode
  use pmc_aero_data
  use pmc_spec_file
  use pmc_gas_state
  use pmc_mpi
  use pmc_stats
  use mpi

  use json_module
  use json_kinds
  use json_parameters, only: unit2str
  use json_string_utilities
  use json_value_module

  implicit none

  type emissions_t
     !> Internally mixed modes.
     type(aero_emission_source_t), allocatable :: mode(:)
  end type emissions_t

  type aero_emission_source_t
      character(len=:),allocatable :: name
      integer :: source_class
      integer :: weight_class
      real(kind=dp), allocatable :: fractions(:,:)
      real(kind=dp), allocatable :: diams(:)
      real(kind=dp), allocatable :: std(:)
  end type aero_emission_source_t

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
  integer :: varid_aero_source
  integer :: status
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
  integer :: i,j,k, i_source, i_mode, iclass, i_type, i_spec
  integer :: n_gas_spec_partmc, n_aero_spec_partmc
  type(aero_data_t) :: aero_data
  type(gas_data_t) :: gas_data
  character(len=9), allocatable :: partmc_species(:)
  character(len=9), allocatable :: smoke_species(:)
  character(len=9), allocatable :: wrf_species(:)

  integer, parameter :: n_species = 5

  ! json
  type(json_file) :: json
  type(json_value),pointer :: p, c, q, r
  type(json_value), pointer :: j_obj
  type(json_value), pointer :: child
  integer :: n_children, n_children_modes
  type(json_core) :: json_obj  ! for manipulating the pointers
  character(len=:),allocatable :: source_name
  real(kind=dp), allocatable :: json_diams(:), json_std(:)
  real(kind=dp), allocatable :: json_fractions(:,:)
  real(kind=dp), allocatable :: fractions_temp(:)
  integer :: source_class, weight_class
  type(emissions_t) :: modes
  logical :: is_found

  ! output structures
  type(aero_dist_t), allocatable :: aero_emission(:)
  real(kind=dp), allocatable, dimension(:) :: aero_emission_rate, &
       aero_emission_time
  real(kind=dp), allocatable, dimension(:) :: gas_emission_rate, &
       gas_emission_time
  integer :: n_mode_counter
  real(kind=dp), allocatable, dimension(:) :: mode_diams, mode_std, &
       mode_num_concs
  integer, allocatable, dimension(:) :: mode_source
  integer, allocatable, dimension(:) :: mode_weight_class
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
  character(len=100), allocatable, dimension(:) :: mode_names
  character(len=100) :: name_i, name_j
  integer :: start_date, n_days_per_file, n_days, i_day, i_date
  logical :: do_point, do_nonpoint, do_bio, do_nonroad, do_mobile, do_rwc, &
      do_ag, do_rail, do_aerosols
  logical :: only_sectional
  real(kind=dp) :: gas_scale_factor, aerosol_scale_factor
  real(kind=dp) :: dx

  integer :: i_send_s, i_send_e
  integer :: rem, is_local, ie_local, n_proc, i_proc, root
  integer :: recv_size_a, recv_size_g
  integer, allocatable, dimension(:) :: displs_a, send_counts_a
  integer, allocatable, dimension(:) :: displs_g, send_counts_g
  integer, allocatable, dimension(:) :: aero_source_centers

!  character(len=200),parameter,dimension(4) :: mobile_sectors = &
!      ["RPH", "RPD", "RPV", "RPP"]
  character(len=200),parameter,dimension(3) :: mobile_sectors = &
      ["RPH", "RPV", "RPD"]

  integer :: stride, block_length
  integer :: grid_emission
  integer :: n_smoke_species

  character(len=(AERO_SOURCE_NAME_LEN)*100) :: aero_source_names
  character(len=100) :: unit

  !
  logical, parameter :: mam3 = .true.
  integer, parameter :: n_emit_species = 5
  integer, parameter :: n_emit_modes = 2
  real(kind=dp), dimension(n_emit_modes) :: diams 
  real(kind=dp), dimension(n_emit_modes) :: stds
  real(kind=dp), dimension(n_emit_modes,n_emit_species) :: factor
  logical, dimension(n_emit_modes,n_emit_species) :: use_species
  real(kind=dp), dimension(n_emit_species) :: aero_spec_density
  real(kind=dp) :: diam, std, num_conc

  call pmc_mpi_init()
  ! aero_species_name = ["PMOTHR", "PSO4  ", "PNO3  ", "POC   ", "PEC   "]

  if (pmc_mpi_rank() == 0) then
     ! only the root process accesses the commandline

     if (command_argument_count() /= 1) then
        call print_usage()
        call die_msg(739173198, "invalid commandline arguments")
     end if

     call get_command_argument(1, spec_name)
  end if

  call pmc_mpi_bcast_string(spec_name)

  call spec_file_open(spec_name, file)

  ! read in aero_data
  call spec_file_read_string(file, 'aerosol_data', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_data(sub_file, aero_data)
  call spec_file_close(sub_file)

  ! read in source names
  call spec_file_read_string(file, 'emission_data', sub_filename)
  call json%initialize()
  call json%load_file(sub_filename)

  call json%info('sources', is_found, n_children=n_children)

  call json%get('sources', p)  ! get a pointer to child
  call json_obj%info(p,n_children=n_children)
  allocate(modes%mode(n_children))
  do i=1,n_children
     call json_obj%get_child(p,i,c) ! get pointer to ith element
     call json_obj%get(c,"source_name", source_name, is_found)
     call json_obj%get(c,"source_class", source_class, is_found)
     call json_obj%get(c,"weight_class", weight_class, is_found)
     call json_obj%get(c,"modes", q)
     call json_obj%info(q,n_children=n_children_modes)

     allocate(json_diams(n_children_modes))
     allocate(json_std(n_children_modes))
     allocate(json_fractions(n_children_modes,n_emit_species))
     do j = 1,n_children_modes
        call json_obj%get_child(q,j,r)
        call json_obj%get(r,"diameter", json_diams(j), is_found)
        call json_obj%get(r,"std", json_std(j), is_found)
        call json_obj%get(r,"fractions", fractions_temp, is_found)
        json_fractions(j,:) = fractions_temp
      end do

     modes%mode(i)%name = source_name
     modes%mode(i)%source_class = source_class
     modes%mode(i)%weight_class = weight_class
     modes%mode(i)%fractions = json_fractions
     modes%mode(i)%diams = json_diams
     modes%mode(i)%std = json_std
     deallocate(json_diams)
     deallocate(json_std)
     deallocate(json_fractions)
  end do

  ! read in gas_data
  call spec_file_read_string(file, 'gas_data', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_gas_data(sub_file, gas_data)
  call spec_file_close(sub_file)

  ! Read in emissions files
  call spec_file_read_string(file, 'grid_name', grid_name)
  call spec_file_read_string(file, 'case_name', case_name)
  call spec_file_read_string(file, 'dir_path', dir_path)
  call spec_file_read_integer(file, 'start_date', start_date)
  call spec_file_read_integer(file, 'n_days_per_file', n_days_per_file)
  call spec_file_read_integer(file, 'n_days', n_days)

  ! Point, nonpt, nonroad, biogenic, mobile
  call spec_file_read_logical(file, 'do_point', do_point)
  call spec_file_read_logical(file, 'do_nonpoint', do_nonpoint)
  call spec_file_read_logical(file, 'do_ag', do_ag)
  call spec_file_read_logical(file, 'do_rwc', do_rwc)
  call spec_file_read_logical(file, 'do_nonroad', do_nonroad)
  call spec_file_read_logical(file, 'do_rail', do_rail)
  call spec_file_read_logical(file, 'do_mobile', do_mobile)
  call spec_file_read_logical(file, 'do_bio', do_bio)

  ! Get grid data
  call spec_file_read_real(file, 'dx', dx)
  call spec_file_read_integer(file, 'nx', nx)
  call spec_file_read_integer(file, 'ny', ny)
  call spec_file_read_integer(file, 'nz', nz)
  call spec_file_read_integer(file, 'nt', nt)

  call spec_file_read_string(file, 'output_prefix', &
       prefix)

  call spec_file_read_logical(file,'only_sectional',only_sectional)

  call species_mapping(smoke_species, partmc_species, wrf_species)

  n_aero_spec_partmc = aero_data_n_spec(aero_data)
  n_gas_spec_partmc = gas_data_n_spec(gas_data)

  n_source_classes = size(modes%mode)

  allocate(aero_source_centers(n_source_classes))

  do i_source = 1,n_source_classes
     aero_source_centers(i_source) = i_source
  end do

  unit='(1)'
  aero_source_names = ""
  do i_source = 1,n_source_classes
     aero_source_names((len_trim(aero_source_names) + 1):) &
          = trim(modes%mode(i_source)%name) !trim(emission_data%name(i_source))
     if (i_source < n_source_classes) then
        aero_source_names((len_trim(aero_source_names) + 1):) = ","
     end if
  end do

!  if (pmc_mpi_rank() == 0) then
     allocate(aero_emissions(nx, ny, nz, nt, n_source_classes, n_species))
     allocate(gas_emissions(nx, ny, nz, nt, n_gas_spec_partmc))
!  end if

  n_smoke_species = size(smoke_species,1)

  ! For memory reasons, only process 0 will read the data.
  if (pmc_mpi_rank() == 0) then
     if (do_point) then
        do_aerosols = .true.
        write(filename_grid, '(a,a,a,a,a,a,a)') trim(dir_path), &
             "source_groups_out.point.", &
             trim(grid_name), '.', trim(case_name), '.ncf'

        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "sginlnts_l.point.", start_date + i_date, '.', &
             i_day, '.', trim(grid_name), '.', trim(case_name),'.ncf'

        write(*,*) 'Reading point emissions', trim(filename_grid), &
             trim(filename_emissions)

        call read_emissions(filename_emissions, filename_grid, & 
             aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, & 
             n_source_classes, nx, ny, nz, nt, gas_data, smoke_species, &
             partmc_species, n_smoke_species, i_date, do_aerosols)
        end do
     end if

     if (do_nonpoint) then
        do_aerosols = .true.
        write(filename_grid, '(a,a,a,a,a,a,a)') trim(dir_path), &
             "source_groups_out.nonpt.", trim(grid_name), '.', &
             trim(case_name), '.ncf'
        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "sginlnts_l.nonpt.", start_date + i_date, '.', i_day, &
             '.',trim(grid_name),'.', trim(case_name),'.ncf'

        write(*,*) 'Reading nonpoint emissions', trim(filename_grid), &
             trim(filename_emissions)

        call read_emissions(filename_emissions, filename_grid, &
             aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, &
             n_source_classes, nx, ny,nz,nt, gas_data,smoke_species, &
             partmc_species, n_smoke_species, i_date, do_aerosols)
        end do
     end if

     if (do_ag) then
        do_aerosols = .false.
        write(filename_grid, '(a,a,a,a,a,a,a)') trim(dir_path), &
             "source_groups_out.ag.", trim(grid_name), '.', &
             trim(case_name), '.ncf'
        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "sginlnts_l.ag.", start_date + i_date, '.', &
             i_day, '.',trim(grid_name),'.', trim(case_name),'.ncf'

        write(*,*) 'Reading nonpoint emissions', trim(filename_grid), &
             trim(filename_emissions)

        call read_emissions(filename_emissions, filename_grid, &
             aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, &
             n_source_classes, nx, ny,nz,nt, gas_data,smoke_species, &
             partmc_species, n_smoke_species, i_date, do_aerosols)
        end do
     end if

     if (do_rwc) then
        do_aerosols = .true.
        write(filename_grid, '(a,a,a,a,a,a,a)') trim(dir_path), &
             "source_groups_out.rwc.", trim(grid_name), '.', &
             trim(case_name), '.ncf'
        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "sginlnts_l.rwc.", start_date + i_date, '.', &
             i_day,'.',trim(grid_name),'.', trim(case_name),'.ncf'

        write(*,*) 'Reading nonpoint emissions', trim(filename_grid), &
             trim(filename_emissions)

        call read_emissions(filename_emissions, filename_grid, &
             aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, &
             n_source_classes, nx, ny,nz,nt, gas_data,smoke_species, &
             partmc_species, n_smoke_species, i_date, do_aerosols)
        end do
     end if

     if (do_nonroad) then
        do_aerosols = .true.
        write(filename_grid, '(a,a,a,a,a,a,a)') trim(dir_path), &
             "source_groups_out.nonroad.", trim(grid_name), '.', &
             trim(case_name), '.ncf'
        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "sginlnts_l.nonroad.", start_date + i_date, '.', i_day, &
             '.', trim(grid_name),'.', trim(case_name), '.ncf'
        write(*,*) 'Reading nonroad emissions', trim(filename_grid), ' ', &
             trim(filename_emissions)
        call read_emissions(filename_emissions, filename_grid, &
             aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, & 
             n_source_classes, nx, ny, nz, nt, gas_data, smoke_species, &
             partmc_species, n_smoke_species, i_date, do_aerosols)
        end do
     end if

     if (do_rail) then
        do_aerosols = .true.
        write(filename_grid, '(a,a,a,a,a,a,a)') trim(dir_path), &
             "source_groups_out.rail.", trim(grid_name), '.', &
             trim(case_name), '.ncf'
        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "sginlnts_l.rail.", start_date + i_date, '.', &
             i_day,'.',trim(grid_name),'.', trim(case_name),'.ncf'

        write(*,*) 'Reading nonpoint emissions', trim(filename_grid), &
             trim(filename_emissions)

        call read_emissions(filename_emissions, filename_grid, &
             aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, &
             n_source_classes, nx, ny,nz,nt, gas_data,smoke_species, &
             partmc_species, n_smoke_species, i_date, do_aerosols)
        end do
     end if

     if (do_mobile) then
        do_aerosols = .true.
        do i_type = 1,size(mobile_sectors)
           write(filename_grid, '(a,a,a,a,a,a,a,a,a)') trim(dir_path), &
                "source_groups_out.", &
                trim(mobile_sectors(i_type)),'.',trim(grid_name), '.', &
                trim(case_name), '.ncf'
           i_day = 1 
           do i_date = 0,n_days-1
           write(filename_emissions, '(a,a,a,a,i8,a,i1,a,a,a,a,a,a)') &
                trim(dir_path), "sginlnts_l.", &
                trim(mobile_sectors(i_type)),'.', start_date + i_date,'.', i_day,'.', &
                trim(grid_name), '.', trim(case_name),'.ncf' 
           write(*,*) 'Reading mobile emissions', trim(filename_grid), &
                trim(filename_emissions)
           call read_emissions(filename_emissions, filename_grid, &
                aero_emissions, gas_emissions, n_gas_spec_partmc, n_species, &
                n_source_classes, nx, ny, nz, nt, gas_data, smoke_species, &
                partmc_species, n_smoke_species, i_date, do_aerosols)
           end do
        end do
     end if

     if (do_bio) then
        i_day = 1
        do i_date = 0,n_days-1
        write(filename_emissions, '(a,a,i8,a,i1,a,a,a,a,a,a)') &
             trim(dir_path), "b3gts_l.", start_date + i_date,'.', i_day, '.', &
             trim(grid_name), '.', trim(case_name),'.ncf'
        write(*,*) 'Reading biogenics', trim(filename_emissions)
        call read_biogenics(filename_emissions, gas_emissions, &
             n_gas_spec_partmc, nx, ny, nz, nt, gas_data, smoke_species, &
             partmc_species, n_smoke_species, i_date)
        end do
     end if
  end if

  ! If we want to zero out any emission classes
  if (pmc_mpi_rank() == 0) then
     !call remove_rare_source_classes(aero_emissions, n_source_classes, nx, ny, &
     !     nz, nt, n_species, emission_data)

     ! We should make some summary statistics
     call create_statistics(gas_emissions, aero_emissions, n_source_classes, &
          nx, ny, nz, nt, n_species, gas_data, gas_data_n_spec(gas_data), dx)
  end if

  call pmc_mpi_barrier()

  n_proc = pmc_mpi_size()
  rem = mod(nx, n_proc)
  if (pmc_mpi_rank() == 0) then
     is_local = pmc_mpi_rank() * (nx / n_proc) + 1
  else
     is_local = pmc_mpi_rank() * (nx / n_proc) + rem + 1
  end if
  ie_local =( pmc_mpi_rank() + 1)  * (nx / n_proc) + rem

  print*, pmc_mpi_rank(), is_local, ie_local

  allocate(aero_emissions_patch(ie_local - is_local + 1, ny, nz, nt, &
       n_source_classes, n_species))
  allocate(gas_emissions_patch(ie_local - is_local + 1, ny, nz, nt, &
       n_gas_spec_partmc))

  if (pmc_mpi_rank() == 0) then
     call make_sectional_emissions(gas_emissions, aero_emissions, nx, ny, nz, &
          nt, n_source_classes, n_species, n_gas_spec_partmc, dx, gas_data)
  end if

  if (only_sectional) then
     call pmc_mpi_finalize()
     stop
  end if

  send_size_a = n_species * n_source_classes * nx * ny * nz * nt
  call MPI_bcast(aero_emissions, send_size_a, MPI_DOUBLE, 0, MPI_COMM_WORLD, &
       ierr)
  send_size_g = nx * ny * nz * nt * n_gas_spec_partmc
  call MPI_bcast(gas_emissions, send_size_g, MPI_DOUBLE, 0, MPI_COMM_WORLD, &
       ierr)

  allocate(aero_emission_time(nt))
  allocate(gas_emission_time(nt))
  ! Hourly emissions
  do i = 1,nt
     aero_emission_time(i) = (i-1)*3600.0d0
     gas_emission_time(i) = (i-1)*3600.0d0
  end do

  ! SMOKE gas units: moles s^-1
  ! WRF-PartMC units: (mol m^{-2} s^{-1})
  gas_scale_factor = 1.0d0 / (dx * dx)
  ! SMOKE aerosol units: g s^-1
  ! WRF-PartMC wants g s^-1 m^-2
  aerosol_scale_factor =  1.0d0 / (dx*dx)

  ! Loop over the grid cells
  do i = is_local, ie_local ! 1,nx
     do j = 1,ny
        do k = 1,nz
           allocate(aero_emission(nt))
           allocate(gas_emission(nt))
           do t = 1,nt
              ! Track the number of modes for this grid cell
              n_mode_counter = 0
              allocate(mode_names(0))
              allocate(mode_diams(0))
              allocate(mode_std(0))
              allocate(mode_num_concs(0))
              allocate(mode_source(0))
              allocate(mode_weight_class(0))
              allocate(mode_vol_fracs(0,20))
              do s = 1,n_source_classes
                 total_mass = sum(aero_emissions(i,j,k,t,s,:)) * &
                     aerosol_scale_factor
                 ! If we want to exclude modes, make this greater than zero
                 if (total_mass >= 0.0d0) then
                    do i_mode = 1,size(modes%mode(s)%diams)
                       diam = modes%mode(s)%diams(i_mode)
                       std = modes%mode(s)%std(i_mode)
                       total_mass = 0.0d0
                       num_conc = 0.0d0
                       do i_spec =1,n_emit_species
                          if (modes%mode(s)%fractions(i_mode,i_spec) > 0.0d0) then
                             num_conc = num_conc + get_num_conc( &
                                  aero_emissions(i,j,k,t,s,i_spec), &
                                  diam, std, aero_spec_density(i_spec)) &
                                  * aerosol_scale_factor &
                                  * modes%mode(s)%fractions(i_mode, i_spec)
                          end if
                       end do
                       call get_vol_frac(aero_emissions(i,j,k,t,s,:), vol_frac, &
                            n_aero_spec_partmc, aero_data, &
                            modes%mode(s)%fractions(i_mode,:))
                       temp_size = shape(mode_vol_fracs)
                       temp_array = mode_vol_fracs
                       deallocate(mode_vol_fracs)
                       allocate(mode_vol_fracs(temp_size(1)+1,temp_size(2)))
                       mode_vol_fracs(1:temp_size(1),:) = temp_array
                       mode_vol_fracs(temp_size(1)+1,:) = vol_frac
                       write(name,'(a,I2)') trim(modes%mode(s)%name), i_mode 
                       if (s == 1 .and. i_mode == 1) then
                          mode_names = [name_i(1:100)]
                       else
                          mode_names = [mode_names, name(1:100)]
                       end if
                       mode_diams = [mode_diams, diam]
                       mode_std = [mode_std, std]
                       mode_num_concs = [mode_num_concs, num_conc]
                       mode_source = [mode_source, modes%mode(s)%source_class]
                       mode_weight_class = [mode_weight_class, &
                            modes%mode(s)%weight_class]
                       n_mode_counter = n_mode_counter + 1
                    end do
                 end if
              end do ! end of source classes

              ! We don't have this information anywhere
              allocate(mode_vol_frac_std(n_mode_counter,20))
              mode_vol_frac_std = 0.0d0
              ! Construct the aero_emissions for this time step 

              ! Allocate the correct size
              allocate(aero_emission(t)%mode(n_mode_counter))
              ! Set the values from the temporary arrays
              do i_mode = 1,n_mode_counter
                 ! FIXME: this needs to be changed
                 aero_emission(t)%mode(i_mode)%name = mode_names(i_mode)(1:100)
                 aero_emission(t)%mode(i_mode)%type = AERO_MODE_TYPE_LOG_NORMAL
                 aero_emission(t)%mode(i_mode)%char_radius = &
                      mode_diams(i_mode) / 2.0
                 aero_emission(t)%mode(i_mode)%log10_std_dev_radius = & 
                      dlog10(mode_std(i_mode))
                 aero_emission(t)%mode(i_mode)%num_conc = &
                      mode_num_concs(i_mode)
                 aero_emission(t)%mode(i_mode)%vol_frac = &
                      mode_vol_fracs(i_mode,:)
                 aero_emission(t)%mode(i_mode)%vol_frac_std = &
                      mode_vol_frac_std(i_mode,:) 
                 aero_emission(t)%mode(i_mode)%source = mode_source(i_mode) 
                 ! FIXME: Add weight class here
              end do
              ! Deallocate single time step variables
              deallocate(mode_diams)
              deallocate(mode_std)
              deallocate(mode_num_concs)
              deallocate(mode_source)
              deallocate(mode_weight_class)
              deallocate(mode_vol_fracs)
              deallocate(mode_vol_frac_std)
              deallocate(mode_names)

              call gas_state_set_size(gas_emission(t), n_gas_spec_partmc)
              gas_emission(t)%mix_rat = gas_emissions(i,j,k,t,:) * &
                   gas_scale_factor

           end do ! end time

           allocate(aero_emission_rate(nt))
           aero_emission_rate = 1.0d0

           ! Output to the NetCDF file
           write(filename, '(a,i3.3,a,i3.3,a,i3.3,a)') trim(prefix),(i), &
                '_',j,'_',k,trim(suffix)
           call pmc_nc_open_write(filename, ncid)
           call pmc_nc_check(nf90_redef(ncid))
           call pmc_nc_check(nf90_def_dim(ncid, "n_times", &
                nt, dimid_times))
           call pmc_nc_check(nf90_enddef(ncid))
           call pmc_nc_write_real_1d(ncid, aero_emission_rate, &
                 "aero_emission_rate_scale", (/ dimid_times /), &
                 unit="(1)", &
                 long_name="Aerosol emission rate", &
                 description="Aerosol emission rate scales at set-points")
           call pmc_nc_write_real_1d(ncid, aero_emission_time, &
                 "aero_emission_time", (/ dimid_times /), &
                 unit="s", &
                 long_name="Aerosol emission time", &
                 description="Aerosol emission set-points times (s).")

           ! Output the emissions data
           call aero_dist_output_netcdf(aero_emission,ncid)

           ! JHC: This will get put into the above subroutine call eventually
           allocate(mode_weight_class(n_source_classes * 2))
           do s = 1,n_source_classes
              mode_weight_class(2*s - 1) = &
                   modes%mode(s)%weight_class
              mode_weight_class(2*s) = &
                   modes%mode(s)%weight_class
           end do
           status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
           call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, check_name, &
                check_dim_size))
           call pmc_nc_write_integer_1d(ncid, mode_weight_class, &
                 "source_weight_class", (/ dimid_n_modes /), &
                 unit="(1)", &
                 long_name="Aerosol weight class", &
                 description="Weight class ID for each aerosol mode.")
           deallocate(mode_weight_class)

           ! Output class names
           call pmc_nc_check(nf90_redef(ncid))
           call pmc_nc_check(nf90_def_dim(ncid, "n_aero_sources", &
                n_source_classes, dimid_aero_sources))
           call pmc_nc_check(nf90_def_var(ncid, "aero_source", NF90_INT, &
                dimid_aero_sources, varid_aero_source))
           call pmc_nc_check(nf90_put_att(ncid, varid_aero_source, "names", &
                aero_source_names))
           call pmc_nc_check(nf90_put_att(ncid, varid_aero_source, "description", &
                "dummy dimension variable (no useful value) - read source names " &
                // "as comma-separated values from the 'names' attribute"))
           call pmc_nc_check(nf90_enddef(ncid))
           call pmc_nc_check(nf90_put_var(ncid, varid_aero_source, &
                aero_source_centers))

           ! Output the gas emissions
           allocate(gas_emission_rate(nt))
           gas_emission_time = aero_emission_time 
           call pmc_nc_write_real_1d(ncid, gas_emission_rate, &
                 "gas_emission_rate_scale", (/ dimid_times /), &
                 unit="(1)", &
                 long_name="Gas emission rate scale factor", &
                 description="Gas emission rate scales at set-points")
           call pmc_nc_write_real_1d(ncid, gas_emission_time, &
                 "gas_emission_time", (/ dimid_times /), &
                 unit="s", &
                 long_name="Gas emission time", &
                 description="Gas emission set-points times (s).")

           call gas_emissions_output_netcdf(gas_emission, ncid)

           call pmc_nc_check(nf90_close(ncid))

           deallocate(aero_emission_rate)
           deallocate(aero_emission)
           deallocate(gas_emission_rate)
           deallocate(gas_emission)
        end do
     end do
  end do

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the shape of the 4D array.
  subroutine get_array_dimensions(ncid, dim_size)

    !> NetCDF file ID.
    integer, intent(in) :: ncid
    !> Dimensions of the input array.
    integer, intent(out) :: dim_size(4)

    integer :: dimids(4)
    integer :: varid, dimid

    ! name isn't important
    ! Find the varid
    status = nf90_inq_varid(ncid, name, varid)
    varid = 2
    ! Inquire for the dimension IDs
    call pmc_nc_check(nf90_inquire_variable(ncid, varid, dimids=dimids))
    ! Loop over all the dimensions to get their lengths
    do dim = 1, 4
       call pmc_nc_check(nf90_inquire_dimension(ncid, dimids(dim), &
            len=dim_size(dim)))
    end do

  end subroutine get_array_dimensions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subroutine to output to a file the aerosol emission.
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

    character(len=(AERO_SOURCE_NAME_LEN)*100) :: aero_source_names

    integer :: start(1), count(1)
    integer :: start2(2), count2(2)
    integer :: start3(3), count3(3)

    integer :: varid_aero_source , dimid_aero_source

    status = nf90_inq_dimid(ncid, "n_times", dimid_n_times)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_times, check_name, &
            check_dim_size))
    n_times = check_dim_size

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

    name='char_radius'
    unit='m'
    long_name = 'characteristic_radius'
    standard_name = 'characteristic_radius'
    description = 'Characteristic radius, with meaning dependent on mode' &
         // ' type (m)'
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
 
    name='source_id'
    unit='(1)'
    long_name = 'Source number'
    standard_name = 'Source number'
    description = 'Source ID number for each emission mode. Maps to names in ' &
         // 'aero_source'
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

    name='num_conc'
    unit='# m^{-2} s^{-1}'
    long_name = 'total number concentration flux'
    standard_name = 'total number concentration flux'
    description = 'Total number concentration flux of mode (#/m^2/s^1).'
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
  real(kind=dp) function get_num_conc(mass, diam, std, aero_spec_density)

    !> Mass concentrations (g s^-1 m^-2).
    real(kind=dp) :: mass
    !> Geometric mean diameter (m).
    real(kind=dp) :: diam
    !> Geometric standard deviation (1).
    real(kind=dp) :: std
    !> Aerosol species density (kg m^-3).
    real(kind=dp) :: aero_spec_density

    integer :: k
    real(kind=dp) :: density, tmp

    k = 3
    density = aero_spec_density * 1000.0d0 ! g / m^3
    tmp = density * (const%pi / 6.0d0) * diam**k * &
         exp((k**2.0d0 / 2.0d0) * log(std)**2.0d0)
    get_num_conc = mass / tmp

  end function get_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine get_vol_frac(species, vol_frac, n_specs, aero_data, factor)

    !> Masses for all species.
    real(kind=dp), dimension(5) :: species
    !> Volume fractions.
    real(kind=dp), allocatable :: vol_frac(:)
    !> Number of aerosol species.
    integer :: n_specs
    !> Aerosol data.
    type(aero_data_t) :: aero_data
    !>
    real(kind=dp), dimension(5) :: factor 

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
          if (factor(i) > 0.0d0) then
             vol_frac(spec_index(i)) = factor(i)*species(i)
          end if
       end do
       vol_frac = vol_frac / aero_data%density
       tot_vol_frac = sum(vol_frac)

       ! convert mass fractions to volume fractions
       vol_frac = vol_frac / tot_vol_frac
       !vol_frac_std = vol_frac_std / aero_data%density
   else ! There isn't anything for this mode so set a dummy value
       vol_frac = 1.0d0
       vol_frac = vol_frac / sum(vol_frac)
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
    name='gas_emission'
    unit='mol m^{-2} s^{-1}'
    long_name = 'gas emissions'
    standard_name = 'gas emissions'
    description = 'gas phase emission rates.'
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

  !> Read emissions from SMOKE.
  subroutine read_emissions(filename_emissions, filename_grid, &
       aero_emissions, gas_emissions, n_gas_species, n_aero_species, &
       n_source_classes, nx, ny, nz, nt, gas_data, smoke_species, &
       partmc_species, n_emission_species, i_date, do_aerosols)

    !> Filename of the emissions data.
    character(len=200), intent(in) :: filename_emissions
    !> Filename of the gridding data.
    character(len=200), intent(in) :: filename_grid
    !> Aerosol emisisons array with source tracking.
    real(kind=dp), dimension(nx,ny,nz,nt,n_source_classes,n_aero_species) :: &
        aero_emissions
    !> Gas emissions array with no source tracking
    real(kind=dp), dimension(nx,ny,nz,nt,n_gas_species) :: &
        gas_emissions
    !> Number of gas species (currently in PartMC).
    integer :: n_gas_species
    !> Number of aerosol species (currently in emissions file).
    integer :: n_aero_species
    !> Number of tracked sources.
    integer :: n_source_classes
    !> East-west dimension size of grid.
    integer :: nx
    !> North-south dimension size of grid.
    integer :: ny
    !> Vertical dimension size of grid.
    integer :: nz
    !> Number of timesteps.
    integer :: nt
    !> Gas data type.
    type(gas_data_t) :: gas_data
    !>
    character(len=9), dimension(n_emission_species), intent(in) :: &
         smoke_species
    !>
    character(len=9), dimension(n_emission_species), intent(in) :: &
         partmc_species
    !>
    integer :: n_emission_species
    integer :: i_date
    logical, intent(in) :: do_aerosols

    real(kind=dp), allocatable, dimension(:,:,:,:,:) :: &
        aero_emissions_source, gas_emissions_source
    character(len=100), parameter, dimension(5) :: &
        aero_species_name = ["PMOTHR", "PSO4  ", "PNO3  ", "POC   ", "PEC   "]
    character(len=100), parameter, dimension(19) :: &
        smoke_gas_species_name = [ &
        "SO2     ", "NO      ", "ALD2    ", &
        "FORM    ", "        ", "NH3     ", &
        "        ", "        ", "        ", &
        "ETH     ", "CO      ", "        ", &
        "OLE     ", "IOLE    ", "TOL     ", &
        "XYL     ", "        ", "        ", &
        "ISOP    "]

    character(len=100), parameter, dimension(19) :: &
        pmc_gas_species_name = [ &
        "SO2     ", "NO      ", "ALD2    ", &
        "HCHO    ", "        ", "NH3     ", &
        "        ", "        ", "        ", &
        "ETH     ", "CO      ", "        ", &
        "OLET    ", "OLEI    ", "TOL     ", &
        "XYL     ", "        ", "        ", &
        "ISOP    "]

    real(kind=dp), allocatable, dimension(:,:,:,:) :: temp_species
    integer, allocatable, dimension(:,:,:,:) :: IGROUP, ROW, COL
    integer :: i_spec, i_spec_partmc, i_class
    integer :: n_sources(1)
    integer :: ncid_grid, ncid_emissions
    integer :: len_aero_species_list(1), len_gas_species_list(1)
    integer :: pmc_gas_index, varid
    character(len=100) :: name, pmc_name
    integer :: t_start, t_end, smoke_start, smoke_end
    
    if (i_date == 0) then
       t_start = 1
       t_end = t_start + 24
       smoke_start = 1
    else
       t_start = i_date * 24 + 2 
       t_end = t_start + 23
       smoke_start = 2
    end if
    smoke_end = 25

    len_aero_species_list = shape(aero_species_name)
    len_gas_species_list = shape(smoke_species)

    call pmc_nc_check( nf90_open(filename_grid, NF90_NOWRITE, &
         ncid_grid))
    call pmc_nc_check( nf90_open(filename_emissions, NF90_NOWRITE, &
         ncid_emissions))
    name = "IGROUP"
    call pmc_nc_read_integer_4d(ncid_grid, IGROUP, name, .true.)
    name = "ROW"
    call pmc_nc_read_integer_4d(ncid_grid, ROW, name, .true.)
    name = "COL"
    call pmc_nc_read_integer_4d(ncid_grid, COL, name, .true.)

    call get_array_dimensions(ncid_emissions, var_size)
  
    ! Allocate the big array - we will put species information in it
    allocate(aero_emissions_source(len_aero_species_list(1),var_size(1), &
         var_size(2), var_size(3), var_size(4)))
    allocate(temp_species(var_size(1), var_size(2), var_size(3), var_size(4)))

    if (do_aerosols) then
    do i_spec = 1,len_aero_species_list(1)
       name = trim(aero_species_name(i_spec))
       call pmc_nc_read_real_4d(ncid_emissions, temp_species, name, .true.)
       aero_emissions_source(i_spec,:,:,:,:) = temp_species
    end do
    end if

    allocate(gas_emissions_source(n_gas_species,var_size(1), var_size(2), &
         var_size(3), var_size(4)))
    do i_spec = 1, n_emission_species
       if (partmc_species(i_spec) .ne. "         ") then
          name = smoke_species(i_spec)
          pmc_name = partmc_species(i_spec)
          pmc_gas_index = gas_data_spec_by_name(gas_data, pmc_name)
          write(*,'(A,A9,1X,A,A9,1X,A,I4)') 'reading SMOKE variable: ', &
              trim(name), 'PartMC variable: ',trim(pmc_name), 'Index = ', &
              pmc_gas_index

          status = nf90_inq_varid(ncid_emissions, name, varid)
          if (status /= NF90_ENOTVAR) then
             call pmc_nc_read_real_4d(ncid_emissions, temp_species, name, &
                  .true.)
             print*, 'emissions total =', sum(temp_species)
             gas_emissions_source(pmc_gas_index,:,:,:,:) = temp_species
          end if
       end if
    end do

    n_sources =  shape(IGROUP(1,:,1,1))

    do i_source = 1, n_sources(1)
       ! FIXME: Add 1 to the index
       ! Get the x and y position
       x = COL(1,i_source,1,1) + 1
       y = ROW(1,i_source,1,1) + 1
       z = 1
       i_class = IGROUP(1,i_source,1,1)
       ! FIXME: Why is there ever a zero? Must be a mistake in SMOKE group file
       if (i_class .eq. 0) cycle

       do i_spec = 1,len_aero_species_list(1)
          aero_emissions(x,y,z,t_start:t_end,i_class,i_spec) = &
               aero_emissions(x,y,z,t_start:t_end,i_class,i_spec) + &
               aero_emissions_source(i_spec,1,i_source,1,smoke_start:smoke_end)
       end do
       do i_spec = 1,n_gas_species
          gas_emissions(x,y,z,t_start:t_end,i_spec) =  &
               gas_emissions(x,y,z,t_start:t_end,i_spec) + &
               gas_emissions_source(i_spec,1,i_source,1,smoke_start:smoke_end)
       end do
    end do

  end subroutine read_emissions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the biogenic species from a file. Only consists of gas emissions.
  subroutine read_biogenics(filename_emissions, gas_emissions, n_gas_species, &
       nx, ny, nz, nt, gas_data, smoke_species, partmc_species, &
       n_emission_species, i_date)

    !> Filename of the emissions data.
    character(len=200), intent(in) :: filename_emissions
    !> Gas emissions array with no source tracking
    real(kind=dp), dimension(nx, ny, nz, nt, n_gas_species) :: &
        gas_emissions
    !> Number of gas species (currently in PartMC).
    integer, intent(in) :: n_gas_species
    !> East-west dimension size of grid.
    integer, intent(in) :: nx
    !> North-south dimension size of grid.
    integer, intent(in) :: ny
    !> Vertical dimension size of grid.
    integer, intent(in) :: nz
    !> Number of timesteps.
    integer, intent(in) :: nt
    !> Gas data type.
    type(gas_data_t), intent(in) :: gas_data
    !> SMOKE speciation.
    character(len=9), dimension(n_emission_species), intent(in) :: &
         smoke_species
    !> PartMC speciation.
    character(len=9), dimension(n_emission_species), intent(in) :: &
         partmc_species
    !> Number of emitted species
    integer, intent(in) :: n_emission_species
    integer :: i_date

    integer :: i_spec, pmc_gas_index
    character(len=100) :: name, pmc_name
    integer :: ncid_emissions
    real(kind=dp), allocatable, dimension(:,:,:,:) :: temp_species
    integer :: t_start, t_end, smoke_start, smoke_end

    if (i_date == 0) then
       t_start = 1
       t_end = t_start + 24
       smoke_start = 1
    else
       t_start = i_date * 24 + 2
       t_end = t_start + 23
       smoke_start = 2
    end if
    smoke_end = 25

    call pmc_nc_check( nf90_open(filename_emissions, NF90_NOWRITE, &
         ncid_emissions))

    do i_spec = 1, n_emission_species
       if (partmc_species(i_spec) .ne. "         ") then
          name = smoke_species(i_spec)
          pmc_name = partmc_species(i_spec)
          pmc_gas_index = gas_data_spec_by_name(gas_data, pmc_name)
          write(*,'(A,A9,1X,A,A9,1X,A,I4)') 'reading SMOKE variable: ', &
              trim(name), &
              'PartMC variable: ',trim(pmc_name), 'Index = ', pmc_gas_index

          status = nf90_inq_varid(ncid_emissions, name, varid)
          if (status /= NF90_ENOTVAR) then
             call pmc_nc_read_real_4d(ncid_emissions, temp_species, name, &
                  .true.)
             gas_emissions(2:nx-1,2:ny-1,:,t_start:t_end,pmc_gas_index) = &
                  gas_emissions(2:nx-1,2:ny-1,:,t_start:t_end,pmc_gas_index) +  &
                  temp_species(:,:,:,smoke_start:smoke_end)
          end if
       end if
    end do

    call pmc_nc_check(nf90_close(ncid_emissions))

  end subroutine read_biogenics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates the binary file needed for conversion program for WRF-CHEM.
  subroutine make_sectional_emissions(gas_emissions, aero_emissions, nx, ny, &
       nz, nt, n_source_classes, n_aero_species, n_gas_species, dx, gas_data)

    !> Gas emissions.
    real(kind=dp), dimension(nx, ny, nz, nt, n_gas_species), intent(in) :: &
         gas_emissions
    !> Aerosol emissions.
    real(kind=dp), dimension(nx, ny, nz, nt, n_source_classes, &
         n_aero_species), intent(in) :: aero_emissions
    !> Number of points in x dimension.
    integer, intent(in) :: nx
    !> Number of points in y dimension.
    integer, intent(in) :: ny
    !> Number of points in the z dimension
    integer, intent(in) :: nz
    !> Number of time steps.
    integer, intent(in) :: nt
    !> Number of aerosol source classes.
    integer, intent(in) :: n_source_classes
    !> Number of aerosol species.
    integer, intent(in) :: n_aero_species
    !> Number of gas species.
    integer, intent(in) :: n_gas_species
    !> Grid spacing.
    real(kind=dp), intent(in) :: dx
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data

    integer, parameter :: NRADM=31
    CHARACTER(len=9), DIMENSION(NRADM)     ::  ENAME
    CHARACTER(len=9), DIMENSION(20) :: PMC
    integer :: IHR, i_spec, index
    integer :: pmc_gas_index
    ! EM3RS must not be double precision 
    real(kind=dp), dimension(nx,nz,ny) :: EM3RS
    integer, parameter :: N_PM25 = 10
    integer :: fac_index
    REAL(kind=dp), parameter, dimension(N_PM25) :: pm_factor = (/ .2d0, .8d0, &
         .2d0, .8d0, .2d0, .8d0, .2d0, .8d0, .2d0, .8d0 /)
    REAL(kind=dp), PARAMETER :: MGPG = 1.0d6     ! ug/g
    character(len=PMC_MAX_FILENAME_LEN) :: out_filename, unit_name
    REAL(kind=dp) :: gas_scale_factor, aerosol_scale_factor
    integer :: dimid_nx, dimid_ny, dimid_nz
    DATA ENAME /    &
      'E_SO2 ', 'E_NO  ','E_NO2 ', 'E_ALD ','E_HCHO','E_ORA2',                     &
      'E_NH3 ','E_HC3 ','E_HC5 ','E_HC8 ',                              &
      'E_ETH ','E_CO  ','E_OL2 ','E_OLT ','E_OLI ','E_TOL ','E_XYL ',   &
      'E_KET ','E_CSL ','E_ISO ','E_PM25I','E_PM25J',                   &
      'E_SO4I','E_SO4J','E_NO3I','E_NO3J','E_ORGI','E_ORGJ','E_ECI',    &
      'E_ECJ','E_PM_10'/
!      'e_so2 ','e_no  ','e_ald ','e_hcho','e_ora2',                     &
!      'e_nh3 ','e_hc3 ','e_hc5 ','e_hc8 ',                              &
!      'e_eth ','e_co  ','e_ol2 ','e_olt ','e_oli ','e_tol ','e_xyl ',   &
!      'e_ket ','e_csl ','e_iso ','e_pm25i','e_pm25j',                   &
!      'e_so4i','e_so4j','e_no3i','e_no3j','e_orgi','e_orgj','e_eci',    &
!      'e_ecj','e_pm10'/
    DATA PMC / &
      'SO2      ', 'NO      ', 'NO2      ', 'ALD2     ', &
      'HCHO     ', '        ', 'NH3      ', &
      '         ', 'PAR     ', '         ', &
      'ETH      ', 'CO      ', '         ', &
      'OLET     ', 'OLEI    ', 'TOL      ', &
      'XYL      ', '        ', '         ', &
      'ISOP     ' / 

    ! Open file and write header information
!    OPEN(19,FILE='wrfem_00to12z_d01',FORM='UNFORMATTED')
!    WRITE(19)NRADM
!    WRITE(19)ename

    ! SMOKE gas units: moles s^-1
    ! WRF-Chem units: mol km^-2 hr^-1
    gas_scale_factor = 3600.0d0 / ((dx * dx) / 1000.0**2)
    ! SMOKE aerosol units: g s^-1
    ! WRF-Chem units: microgram m^-3 m s^-1
    aerosol_scale_factor =  MGPG / (dx*dx)
    ! Loop over time
    do IHR = 1,24
       write(out_filename,'(A,I3.3,A)')  "bulk_emissions_", IHR, ".nc"
       call pmc_nc_open_write(out_filename, ncid)
       call pmc_nc_check(nf90_redef(ncid))
       call pmc_nc_check(nf90_def_dim(ncid, "nx", &
            nx, dimid_nx))
       call pmc_nc_check(nf90_def_dim(ncid, "nz", &
            nz, dimid_nz))
       call pmc_nc_check(nf90_def_dim(ncid, "ny", &
            ny, dimid_ny))
!       call pmc_nc_check(nf90_def_dim(ncid, "nz", &
!            nz, dimid_nz))
       call pmc_nc_check(nf90_enddef(ncid)) 
       ! Collapse gas species
       do i_spec = 1, NRADM
         EM3RS = 0.0
         if (i_spec < NRADM - 10) then
            pmc_gas_index = gas_data_spec_by_name(gas_data, PMC(i_spec))
            print*, ename(i_spec), pmc_gas_index
            if (pmc_gas_index > 0) then
            DO I=1,NX
            DO K=1,NZ
            DO J=1,NY
               EM3RS(I,K,J) = gas_emissions(I,J,K,IHR, pmc_gas_index) * &
                    gas_scale_factor
            ENDDO
            ENDDO
            ENDDO
            end if
            unit_name = 'gas'
            call pmc_nc_write_real_3d(ncid, em3rs, ename(i_spec), dim_name_1='nx', &
                 dim_name_2='nz', dim_name_3='ny', unit=unit_name)
         else if (i_spec < NRADM) then
            index = (i_spec - (NRADM - 10)) / 2 + 1
            fac_index = i_spec - (NRADM - 11)
            write(*,*) i_spec, index, ename(i_spec), pm_factor(fac_index)
            DO I=1,NX
            DO K=1,NZ
            DO J=1,NY
               EM3RS(I,K,J) = (sum(aero_emissions(I,J,K,IHR,:,index)) * &
                    pm_factor(fac_index)) * aerosol_scale_factor
            ENDDO
            ENDDO
            ENDDO
            unit_name = 'aerosol'
            call pmc_nc_write_real_3d(ncid, em3rs, ename(i_spec) , dim_name_1='nx', &
                 dim_name_2='nz', dim_name_3='ny', unit=unit_name)
         else
            EM3RS = 0.0
         end if

       end do
       call pmc_nc_close(ncid)
    end do

  end subroutine make_sectional_emissions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads map between different speciations. 
  subroutine species_mapping(smoke_species, partmc_species, wrf_species)

    !> SMOKE speciation table.
    character(len=9), allocatable, intent(inout) :: smoke_species(:)
    !> PartMC speciation table.
    character(len=9), allocatable, intent(inout) :: partmc_species(:)
    !> WRF-CHEM speciation table.
    character(len=9), allocatable, intent(inout) :: wrf_species(:)

    character(len=9) :: a,b,c
    !> Complete line read.
    character(len=255) :: line
    !> True if at EOF.
    logical :: eof
    type(spec_file_t) :: sub_file
    integer :: n_species
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: species_name(:)
    real(kind=dp), allocatable :: species_data(:,:)
    integer :: i_spec 

    call spec_file_open("species_map.dat", sub_file)
    eof = .false.

    n_species = number_species(sub_file)

    allocate(partmc_species(n_species))
    allocate(smoke_species(n_species))
    allocate(wrf_species(n_species))

    do i_spec = 1, n_species
       call spec_file_read_next_data_line(sub_file, line, eof)
       !FIXME: Length is fixed here
       read(line,'(A9,A9,A9)') a,b,c
       smoke_species(i_spec) = a
       partmc_species(i_spec) = b 
       wrf_species(i_spec) = c    
       print*, 'SMOKE: ', a,'PARTMC: ',b,'WRF-CHEM: ',c
    end do

    call spec_file_close(sub_file)

  end subroutine species_mapping

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the number of species in a mappping file.
  integer function number_species(file)

    !> Spec file.
    type(spec_file_t), intent(in) :: file

    character(len=255) :: line
    integer :: num_lines
    logical :: eof
    integer :: max_lines

    num_lines = 0
    eof = .false.
    max_lines = 0

    call spec_file_read_next_data_line(sub_file, line, eof)
    do while (.not. eof)
       num_lines = num_lines + 1
       if (num_lines > SPEC_FILE_MAX_LIST_LINES) then
          call spec_file_die_msg(450564159, file, &
               'maximum number of lines exceeded')
       end if
       if (max_lines > 0) then
          if (num_lines >= max_lines) then
             eof = .true.
          end if
       end if
       if (.not. eof) then
         call spec_file_read_next_data_line(sub_file, line, eof) 
       end if
    end do

    number_species = num_lines

    rewind(unit=sub_file%unit)

  end function number_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !> Find statistics regarding emission groups over the entire domain.
  !! This identifies the less common sources and allows removal if desired.
  subroutine remove_rare_source_classes(aero_emissions, n_source_classes, nx, &
       ny, nz, nt, n_species) !, emission_data)

    !> Aerosol emissions.
    real(kind=dp), dimension(nx, ny, nz, nt, n_source_classes, n_species), &
      intent(inout) :: aero_emissions
    integer, intent(in) :: n_source_classes
    !> East-west dimension size of grid.
    integer, intent(in) :: nx
    !> North-south dimension size of grid.
    integer, intent(in) :: ny
    !> Vertical dimension size of grid.
    integer, intent(in) :: nz
    !> Number of emission times.
    integer, intent(in) :: nt
    !> Number of emission species.
    integer, intent(in) :: n_species

    integer :: i_source
    real(kind=dp), allocatable, dimension(:) :: total_mass
    integer, allocatable, dimension(:) :: total_cells

    allocate(total_cells(n_source_classes))
    allocate(total_mass(n_source_classes))

    do i_source = 1, n_source_classes
       total_cells(i_source) = get_total_cells(aero_emissions(:,:,:,:,i_source,:), &
            nx,ny,nz,nt,n_species)
       total_mass(i_source) = sum(aero_emissions(:,:,:,:,i_source,:))
!       print*, 'source name:', trim(emission_data%name(i_source)), ' ', &
!          'source number:', i_source, 'mass:', total_mass(i_source), &
!          'number of cells:', total_cells(i_source)
    end do

    ! FIXME: Determine criteria for removal
    ! Eventually we should regroup these groups so they are accounted for.
    do i_source = 1, n_source_classes
       if(total_cells(i_source) < ((nx*ny) * 5.0d-4)) then
          aero_emissions(:,:,:,:,i_source,:) = 0.0d0
       end if
    end do

  end subroutine remove_rare_source_classes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the total number of grid cells that contain an emission source.
  integer function get_total_cells(aero_emissions,nx,ny,nz,nt,n_species)

    !> Aerosol emissions.
    real(kind=dp), dimension(nx, ny, nz, nt, n_species) :: aero_emissions
    !> East-west dimension size of grid.
    integer :: nx
    !> North-south dimension size of grid.
    integer :: ny
    !> Vertical dimension size of grid.
    integer :: nz
    !> Number of emission times.
    integer :: nt
    !> Number of emission species.
    integer :: n_species

    integer :: i,j,k

    get_total_cells = 0

    k = 1

    do i = 1,nx
    do j = 1,ny
       if (sum(aero_emissions(i,j,k,:,:)) > 0) then
          get_total_cells = get_total_cells + 1
       end if
    end do
    end do

  end function get_total_cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs summary statistics for the emissions.
  subroutine create_statistics(gas_emissions, aero_emissions, n_source_classes, &
       nx, ny, nz, nt, n_species, gas_data, n_gas_species, dx)

    !> Gas emissions
    real(kind=dp), dimension(nx, ny, nz, nt, n_gas_species), &
         intent(inout) :: gas_emissions
    !> Aerosol emissions.
    real(kind=dp), dimension(nx, ny, nz, nt, n_source_classes, n_species), &
      intent(inout) :: aero_emissions
    integer, intent(in) :: n_source_classes
    !> East-west dimension size of grid.
    integer, intent(in) :: nx
    !> North-south dimension size of grid.
    integer, intent(in) :: ny
    !> Vertical dimension size of grid.
    integer, intent(in) :: nz
    !> Number of emission times.
    integer, intent(in) :: nt
    !> Number of emission species.
    integer, intent(in) :: n_species
    !>
    type(gas_data_t), intent(in) :: gas_data
    integer, intent(in) :: n_gas_species
    !
    real(kind=dp), intent(in) :: dx

    real(kind=dp), dimension(nt) :: source_values
    integer, dimension(n_source_classes) :: total_cells
    real(kind=dp), dimension(n_species) :: masses, mass_frac
    type(stats_1d_t), dimension(n_source_classes) :: mass_fractions
    type(stats_1d_t), dimension(n_source_classes) :: emission_time_series
    integer :: i,j,k,i_source,i_spec,i_time
    real(kind=dp), allocatable :: mass_frac_array(:,:,:)
    real(kind=dp) :: total_mass, max_val, min_val
    character(len=PMC_MAX_FILENAME_LEN) :: out_filename
    character(len=100) :: var_name
    integer :: ncid, dimid
    integer :: source_count(nx,ny)
    real(kind=dp) :: gas_scale_factor

    k = 1
    i_time = 1
    do i_source = 1,n_source_classes
        total_cells(i_source) = get_total_cells(aero_emissions(:,:,:,:,i_source,:), &
            nx,ny,nz,nt,n_species)
        ! Lets sum everything and divide by the number of total cells?
        do i = 1,nx
        do j = 1,ny
        do i_spec = 1,n_species
            masses = aero_emissions(i,j,k,i_time,i_source,:)
            if (sum(masses) > 0.0d0) then
               mass_frac = masses / sum(masses)
               call stats_1d_add(mass_fractions(i_source), mass_frac)
            end if
        end do
        end do
        end do
    end do
    source_count = 0
    do i_source = 1,n_source_classes
        do i = 1,nx
        do j = 1,ny
        do i_time = 1,nt
            masses = aero_emissions(i,j,k,i_time,i_source,:)
            if (sum(masses) > 0.0d0) then
               call stats_1d_add_entry(emission_time_series(i_source),&
                    sum(masses), &
                    i_time)
            end if
        end do
        if (sum(aero_emissions(i,j,k,1,i_source,:)) > 0) then
           source_count(i,j) = source_count(i,j) + 1
        end if
        end do
        end do
    end do

    out_filename = "emission_statistics.nc"
    call pmc_nc_open_write("emission_statistics.nc", ncid)
    do i_source = 1,n_source_classes
       write(var_name,'(A,I2.2)') "source_", i_source
       if (allocated(mass_fractions(i_source)%mean)) then
            call stats_1d_output_netcdf(mass_fractions(i_source), ncid, var_name, &
                 dim_name="n_species", unit="(1)")
    end if
    end do
    do i_source = 1,n_source_classes
       write(var_name,'(A,I2.2)') "emission_", i_source
       if (allocated(mass_fractions(i_source)%mean)) then
            call stats_1d_output_netcdf(emission_time_series(i_source), ncid, &
                 var_name, &
                 dim_name="n_times", unit="(1)")
    end if
    end do

    call pmc_nc_write_integer_2d(ncid, source_count, "source_count", &
       dim_name_1="nx", dim_name_2="ny", unit="(1)", &
       long_name="source_count", standard_name="source_count", &
       description="number of emission sources in each grid cell")

    ! FIXME: Should come up with a better idea here...
    ! We should have a pmc_nc_write_real_3d
    do i_source = 1,n_source_classes
       write(var_name,'(A,I2.2)') "emission_values_", i_source
       call pmc_nc_write_real_2d(ncid, &
          sum(aero_emissions(:,:,1,12,i_source,:),dim=3), &
          var_name, dim_name_1='nx', &
          dim_name_2='ny', unit='',description='emissions average')
    end do

    ! Output gas data
    call gas_data_output_netcdf(gas_data, ncid)

    ! output the gas emissions 
    k = 1
    i_time = 1 
    gas_scale_factor = 1.0d0 / (dx * dx)
!    do i_spec = 1,gas_data_n_spec(gas_data)
!        ! Skip ones that don't exist as emissions
!        if (sum(gas_emissions(:,:,k,:,i_spec)) > 0.0d0) then
!           write(var_name, '(A)') gas_data%name(i_spec)
!           call pmc_nc_write_real_2d(ncid, & 
!              gas_emissions(:,:,k,i_time,i_spec) * gas_scale_factor, &
!              var_name, dim_name_1='nx', dim_name_2='ny', &
!              unit='mol m^{-2} s^{-1}', description='gas emission rate')
!        end if
!    end do    

    call pmc_nc_close(ncid)

    allocate(mass_frac_array(nx,ny,5))
    do i_time = 1,nt
       write(out_filename,'(a,i3.3,a)') 'emissions_', i_time, '.nc'
       call pmc_nc_open_write(out_filename, ncid)
       do i_source = 1,n_source_classes
          write(var_name,'(A,I2.2)') "emission_values_", i_source
          call pmc_nc_write_real_2d(ncid, &
             sum(aero_emissions(:,:,1,i_time,i_source,:),dim=3), &
             var_name, dim_name_1='nx', &
             dim_name_2='ny', unit='', description='emissions')
       end do
     
       do i_source =1, n_source_classes
       mass_frac_array = -1.0
       do i = 1,nx
       do j = 1,ny
            masses = aero_emissions(i,j,1,i_time,i_source,:)
            if (sum(masses) > 0.0d0) then
               mass_frac_array(i,j,:) = masses / sum(masses)
            end if
       end do
       end do
       do i_spec = 1,n_species
          write(var_name,'(A,I2.2,A,I2.2)') "mass_frac_source_", i_source,"_", i_spec
          call pmc_nc_write_real_2d(ncid, mass_frac_array(:,:,i_spec), &
             var_name, dim_name_1='nx', &
             dim_name_2='ny', unit='', description='emissions')
       end do
       end do

       do i_spec = 1,gas_data_n_spec(gas_data)
        if (sum(gas_emissions(:,:,k,:,i_spec)) > 0.0d0) then
           write(var_name, '(A)') gas_data%name(i_spec)
           call pmc_nc_write_real_2d(ncid, &
              gas_emissions(:,:,k,i_time,i_spec) * gas_scale_factor, &
              var_name, dim_name_1='nx', dim_name_2='ny', &
              unit='mol m^{-2} s^{-1}', description='gas emission rate')
        end if
       end do
       call pmc_nc_close(ncid)
   end do

  end subroutine create_statistics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program make_emissions
