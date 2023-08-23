! Copyright (C) 2011-2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc_init module.

!> Initializes PartMC.

module wrf_pmc_init

  use pmc_aero_data
  use pmc_aero_state
  use pmc_aero_dist
  use pmc_aero_weight
  use pmc_spec_file
  use pmc_gas_data
  use pmc_gas_state
  use pmc_scenario
  use pmc_env_state
  use pmc_spec_line
  use pmc_util
  use pmc_rand
  use pmc_mpi

  ! Temporary til some subroutines are moved from partmc_driver
  use wrf_pmc_driver

  ! wrf modules
  use module_model_constants
  use module_wrf_error
  use module_domain_type
  use module_domain
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Maximum length of a file path.
  integer, parameter :: INPUT_FILE_PATH_NAME_LEN = 500

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes PartMC data structures on the WRF domain based on input files
  subroutine init_wrf_partmc(grid, scenario, env_states, aero_data, &
       aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
       pmc_ks, pmc_ke, nx, ny, nz, global_nx, global_ny, global_nz, &
       config_flags)

    !> PartMC east-west start of domain.
    integer, intent(in) :: pmc_is
    !> PartMC east-west end of domain.
    integer, intent(in) :: pmc_ie
    !> PartMC north-south start of domain.
    integer, intent(in) :: pmc_js
    !> PartMC north-south end of domain.
    integer, intent(in) :: pmc_je
    !> PartMC top-bottom start of domain.
    integer, intent(in) :: pmc_ks
    !> PartMC top-bottom end of domain.
    integer, intent(in) :: pmc_ke
    !> WRF domain.
    type(domain), intent(inout) :: grid
    !> Environment data.
    type(scenario_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenario
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data 
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout):: aero_states
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data 
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: gas_states
    !> Number of boxes in east-west direction.
    integer, intent(in) :: nx
    !> Number of boxes in the north-south direction.
    integer, intent(in) :: ny
    !> Number of boxes in the top-bottom direction.
    integer, intent(in) :: nz
    !> East-west domain dimension size of WRF domain.
    integer, intent(in) :: global_nx
    !> North-south domain dimension size of WRF domain.
    integer, intent(in) :: global_ny
    !> Top-bottom domain dimension size of WRF domain.
    integer, intent(in) :: global_nz
    !> WRF namelist configuration.
    type(grid_config_rec_type), intent(in) :: config_flags

    ! Local variables
    integer :: i
    integer :: j
    integer :: k
    integer :: p
    real(kind=dp) :: alpha
    ! UUID
    character(len=PMC_UUID_LEN) :: uuid
    ! Counting particles
    integer :: n_parts_added
    real(kind=dp) :: n_parts
    ! Local PartMC structures
    type(spec_file_t) :: sub_file
    character(LEN=INPUT_FILE_PATH_NAME_LEN) :: filename
    character(LEN=:), allocatable :: name
    ! Various structures to store inputs
    type(gas_data_t) :: input_gas_data
    type(aero_data_t) :: input_aero_data
    type(aero_data_t) :: global_aero_data
    character(len=AERO_MODE_NAME_LEN) :: mode_name
    integer :: n_ic, n_emit, i_ic, i_emit

    ! MPI variables
    character, allocatable :: buffer(:)
    character, allocatable :: aero_data_buffer(:)
    integer :: aero_data_size
    integer :: position
    integer :: max_buffer_size
    integer :: buffer_size
    integer :: wrf_i,wrf_j,wrf_k, i_mode, dummy
    real(kind=dp) :: time_ic, time_bc, t1, t2
    integer :: i_start,i_end,j_start,j_end,k_start,k_end
    integer :: n_class

    ! Initialize the random number generator
    call pmc_srand(grid%random_seed, pmc_mpi_rank())

    print*,'Initializing:',pmc_mpi_rank(),pmc_is,pmc_ie,pmc_js,pmc_je,pmc_ks,pmc_ke

    call wrf_message('In the partmc initialization subroutine')

    ! Processor 0 reads in all the input data
    if (pmc_mpi_rank() == 0) then
       ! Load gas_data inputs
       call wrf_message('PartMC_init: Loading gas_data from input files')
       filename = 'gas_data.dat'
       call spec_file_open(filename, sub_file)
       call spec_file_read_gas_data(sub_file, input_gas_data)
       call spec_file_close(sub_file)
       ! Load aero_data inputs
       call wrf_message('PartMC_init: Loading aero_data from input files')
       filename = 'aero_data.dat'
       call spec_file_open(filename,sub_file)
       call spec_file_read_aero_data(sub_file, input_aero_data)
       call fractal_set_spherical(input_aero_data%fractal)
       call spec_file_close(sub_file)

       call get_sources_and_weights(input_aero_data, grid%partmc_ics, &
            grid%partmc_bcs, grid%partmc_emissions, &
            config_flags%periodic_x .and. config_flags%periodic_y)

       ! Get the maximum pack size
       max_buffer_size = pmc_mpi_pack_size_gas_data(input_gas_data) &
            +  pmc_mpi_pack_size_aero_data(input_aero_data)
       allocate(buffer(max_buffer_size))
       ! Pack the buffer
       position = 0
       call pmc_mpi_pack_gas_data(buffer, position, input_gas_data)
       call pmc_mpi_pack_aero_data(buffer, position, input_aero_data)
       buffer_size = position
    end if

    ! Broadcast the buffer size
    call pmc_mpi_bcast_integer(max_buffer_size)

    ! Allocate the buffer on other nodes
    if (pmc_mpi_rank() /= 0) then
       allocate(buffer(max_buffer_size))
    end if

    ! Broadcast the buffer
    call pmc_mpi_bcast_packed(buffer)
    ! Unpack the buffer
    if (pmc_mpi_rank() /= 0) then
       position = 0
       call pmc_mpi_unpack_gas_data(buffer, position, input_gas_data)
       call pmc_mpi_unpack_aero_data(buffer, position, input_aero_data)
    end if
    deallocate(buffer)

    call wrf_message('PartMC_init: Setting initial information for env_states')

    do i = pmc_is, pmc_ie
    do k = pmc_ks, pmc_ke
    do j = pmc_js, pmc_je
        wrf_i = i
        wrf_j = j
        wrf_k = k
        env_states(i,k,j)%ix = wrf_i
        env_states(i,k,j)%iy = wrf_j
        env_states(i,k,j)%iz = wrf_k
        ! Set latitude/longitude
        env_states(i,k,j)%latitude = grid%xlat(wrf_i,wrf_j)
        env_states(i,k,j)%longitude = grid%xlong(wrf_i,wrf_j)
        ! Set an altitude/height in m
        ! WRF initially doesn't have z values so we must use geopotential here
        ! in a more creative way
        ! (Top - bottom)/g is the total height
        env_states(i,k,j)%height = (grid%phb(wrf_i,wrf_k+1,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k+1,wrf_j)-(grid%phb(wrf_i,wrf_k,wrf_j) + &
             grid%ph_1(wrf_i,wrf_k,wrf_j)))/9.8
        ! Find the bottom and add half the height for the center point
        ! (Top - bottom)/2g + (bottom)/g
        env_states(i,k,j)%altitude = 0.5d0*((grid%phb(wrf_i,wrf_k+1,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k+1,wrf_j)-(grid%phb(wrf_i,wrf_k,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k,wrf_j)))/9.8) + &
             ((grid%phb(wrf_i,wrf_k,wrf_j)+grid%ph_1(wrf_i,wrf_k,wrf_j))/9.8)
        env_states(i,k,j)%z_min = (grid%phb(wrf_i,wrf_k,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k,wrf_j))/9.8
        env_states(i,k,j)%z_max = (grid%phb(wrf_i,wrf_k+1,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k+1,wrf_j))/9.8
        ! FIXME: Density or inverse density
        env_states(i,k,j)%rrho = grid%alt(wrf_i,wrf_k,wrf_j)
        ! Set time and date information
        env_states(i,k,j)%start_day = grid%julian
        env_states(i,k,j)%start_time = grid%start_hour * 3600.0d0 + &
             grid%start_minute * 60.0d0 + grid%start_second
    end do
    end do
    end do

    call wrf_message('PartMC_init: Setting the initial meteorological &
         conditions for env_states')

    call init_wrf_to_partmc(grid, env_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
         pmc_ks, pmc_ke)

    call wrf_message('PartMC_init: Setting aerosol data for WRF domain')

    aero_data = input_aero_data

    n_parts = grid%num_particles

    call wrf_message('PartMC_init: Setting gas data for the WRF domain')

    gas_data = input_gas_data

    call pmc_mpi_barrier()

    call wrf_message('PartMC_init: Setting initial aerosol and gas states for &
         WRF domain')

    do j = pmc_js, pmc_je
    do k = pmc_ks, pmc_ke
    do i = pmc_is, pmc_ie
       call aero_state_zero(aero_states(i,k,j))
       call aero_state_set_weight(aero_states(i,k,j), aero_data, &
             AERO_STATE_WEIGHT_FLAT_SPECIFIED)
       call aero_state_set_n_part_ideal(aero_states(i,k,j), n_parts)
       call gas_state_set_size(gas_states(i,k,j), &
            gas_data_n_spec(gas_data))
    end do
    end do
    end do

    call wrf_message('PartMC_init: Setting initial conditions')

    t1 = MPI_Wtime()

    if (config_flags%periodic_x .and. config_flags%periodic_y) then
      j_start = max(pmc_js,1)
      j_end = min(pmc_je,global_ny)
      i_start = max(pmc_is,1)
      i_end = min(pmc_ie,global_nx)
    else if (config_flags%do_restart) then
      j_start = max(pmc_js,1)
      j_end = min(pmc_je,global_ny)
      i_start = max(pmc_is,1)
      i_end = min(pmc_ie,global_nx)
    else
       j_start = max(pmc_js,2)
       j_end = min(pmc_je,global_ny-1)
       i_start = max(pmc_is,2)
       i_end = min(pmc_ie,global_nx-1)
    end if
    k_start = pmc_ks
    k_end = pmc_ke

    if (.not. config_flags%do_restart) then
       do j = j_start,j_end
       do i = i_start,i_end
          call init_read_in_ics(grid, aero_states(i,:,j), gas_states(i,:,j), &
               i, j, nz, aero_data, gas_data, env_states(i,:,j), &
               grid%partmc_ics)
       end do
       end do
    else
       do j = j_start,j_end
       do i = i_start,i_end
          call init_read_in_restart(aero_states(i,:,j), gas_states(i,:,j), i, &
               j, pmc_ke, aero_data, gas_data, config_flags%partmc_restart_prefix, &
               config_flags%partmc_restart_index)
       end do
       end do
    end if

    time_ic = MPI_Wtime() - t1
    call pmc_mpi_barrier()

    call wrf_message('PartMC_init: Setting scenario for the WRF domain')

    if (grid%do_emission) then
       do j = j_start,j_end
       do k = pmc_ks,pmc_ke
       do i = i_start,i_end
          ! Emissions for the lowest kemit layers.
          if (k <= config_flags%kemit) then
            call init_read_in_emissions(scenario(i,k,j), i, j, k, aero_data, &
                 gas_data, grid%partmc_emissions)
          else
             call init_zero_emissions(scenario(i,k,j), aero_data, &
                  gas_data)
          end if
       end do
       end do
       end do
    else
       do j = j_start,j_end
       do k = pmc_ks,pmc_ke
       do i = i_start,i_end
       call init_zero_emissions(scenario(i,k,j), aero_data, &
            gas_data)
       end do
       end do
       end do
    end if

    ! Read in boundary conditions
    ! West boundary
    t1 = MPI_Wtime()
    if (.not. (config_flags%periodic_x .and. config_flags%periodic_y)) then
       print*, 'partmc_init: west boundary', pmc_is
       if (pmc_is == 1) then
          i = 1
          do j = max(pmc_js,1), min(pmc_je,global_ny)
             do k = k_start,k_end !pmc_ks, pmc_ke
                call init_read_in_bcs(grid, scenario(i,k,j), i, j, k, aero_data, &
                     aero_states(i,k,j), gas_states(i,k,j), env_states(i,k,j), &
                     grid%partmc_bcs, config_flags%do_restart)
             end do
          end do
       end if
       print*, 'partmc_init: east boundary', pmc_ie, global_nx
       ! East boundary
       if (pmc_ie == global_nx) then
          i = global_nx
          do j = max(pmc_js,1), min(pmc_je,global_ny)
             do k = k_start,k_end !pmc_ks, pmc_ke
              call init_read_in_bcs(grid, scenario(i,k,j), i, j, k, aero_data, &
                   aero_states(i,k,j), gas_states(i,k,j), env_states(i,k,j), &
                   grid%partmc_bcs, config_flags%do_restart)
             end do
          end do
       end if

       ! South boundary
       print*, 'partmc_init: south boundary', pmc_js
       if (pmc_js == 1) then
          j = 1
          ! Loop changed to avoid corners that are also east/west boundaries
          do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
             do k = k_start,k_end !pmc_ks, pmc_ke
                call init_read_in_bcs(grid, scenario(i,k,j), i, j, k, aero_data, &
                     aero_states(i,k,j), gas_states(i,k,j), env_states(i,k,j), &
                     grid%partmc_bcs, config_flags%do_restart)
             end do
          end do
       end if
       ! North boundary
       print*, 'partmc_init: north boundary', pmc_je, global_ny
       if (pmc_je == global_ny) then
          j = global_ny
          ! Loop changed to avoid corners that are also east/west boundaries
          do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
             do k = k_start,k_end !pmc_ks, pmc_ke
                 call init_read_in_bcs(grid, scenario(i,k,j), i, j, k, aero_data, &
                      aero_states(i,k,j), gas_states(i,k,j), env_states(i,k,j), &
                      grid%partmc_bcs, config_flags%do_restart)
             end do
          end do
       end if
    end if
    time_bc = MPI_Wtime() - t1

    ! UUID
    if (pmc_mpi_rank() == 0) then
      call uuid4_str(uuid)
    end if

    call pmc_mpi_bcast_string(uuid)
    grid%uuid = uuid


    ! Allocate the probability arrays for all cells including boundary
    do i = pmc_is, pmc_ie
    do k = pmc_ks, pmc_ke
    do j = pmc_js, pmc_je
       n_class = aero_weight_array_n_class(aero_states(i,k,j)%awa)
       allocate(env_states(i,k,j)%prob_advection(-1:1,-1:1,-1:1,n_class))
       allocate(env_states(i,k,j)%prob_diffusion(-1:1,-1:1,-1:1,n_class))
       allocate(env_states(i,k,j)%prob_vert_diffusion(pmc_ks:pmc_ke,n_class))
    end do
    end do
    end do

    call pmc_mpi_barrier()

    print*, 'timing initialization - ICs:', time_ic

    ! Put initial gas mixing ratios from chem array into partmc
    call partmc_from_wrf(grid, gas_states, gas_data, pmc_is, pmc_ie, &
       pmc_js, pmc_je, pmc_ks, pmc_ke)

    ! Put initial gas mixing ratios into chem array
!    call partmc_to_wrf(grid, aero_states, aero_data, gas_states, gas_data, &
!         pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke)


    !!! Uncomment to generically add scalars
    !do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
    !do k = pmc_ks, pmc_ke
    !do j = max(pmc_js,2), min(pmc_je,global_ny-1)
    !   grid%chem(i,k,j,p_NUM_CONC_a39) = 5.0e9
    !end do
    !end do
    !end do
    !if (pmc_js == 1) then
    !j = 1
    !do i = pmc_is, pmc_ie
    !do k = pmc_ks, pmc_ke
    !   grid%chem(i,k,j,p_NUM_CONC_a40) = 5.0e9
    !end do
    !end do
    !end if
    !if (pmc_je == global_ny) then
    !j = global_ny
    !do i = pmc_is, pmc_ie
    !do k = pmc_ks, pmc_ke
    !   grid%chem(i,k,j,p_NUM_CONC_a40) = 5.0e9
    !end do
    !end do
    !end if
    !if (pmc_is == 1) then
    !i = pmc_is
    !do j = pmc_js, pmc_je
    !do k = pmc_ks, pmc_ke
    !   grid%chem(i,k,j,p_NUM_CONC_a40) = 5.0e9
    !end do
    !end do
    !end if
    !if  (pmc_ie == global_nx) then
    !i = global_nx
    !do j = pmc_js, pmc_je
    !do k = pmc_ks, pmc_ke
    !   grid%chem(i,k,j,p_NUM_CONC_a40) = 5.0e9
    !end do
    !end do
    !end if

    !do i = pmc_is,pmc_ie
    !do k = pmc_ks,pmc_ke
    !do j = pmc_js,pmc_je
    !    if (j < global_ny/2) then
    !       grid%chem(i,k,j,p_NUM_CONC_a37) = 5.0e9
    !    else
    !       grid%chem(i,k,j,p_NUM_CONC_a38) = 5.0e9
    !    end if
    !end do
    !end do
    !end do

    !do i = pmc_is,pmc_ie
    !do k = pmc_ks,pmc_ke
    !do j = pmc_js,pmc_je
    !    grid%chem(i,k,j,p_NUM_CONC_a36) = 5.0e9
    !end do
    !end do
    !end do


  end subroutine init_wrf_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of WRF variables into partmc env_state variables
  !> Nearly the same as wrf_to_partmc
  subroutine init_wrf_to_partmc(grid, env_states, pmc_is, pmc_ie, pmc_js, &
       pmc_je, pmc_ks, pmc_ke)

    !> PartMC east-west start of domain.
    integer, intent(in) :: pmc_is
    !> PartMC east-west end of domain.
    integer, intent(in) :: pmc_ie
    !> PartMC north-south start of domain.
    integer, intent(in) :: pmc_js
    !> PartMC north-south end of domain.
    integer, intent(in) :: pmc_je
    !> PartMC top-bottom start of domain.
    integer, intent(in) :: pmc_ks
    !> PartMC top-bottom end of domain.
    integer, intent(in) :: pmc_ke
    !> WRF domain.
    type(domain), intent(in) :: grid
    !> Full domain of environmental states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states

    ! Local variables
    integer :: i,j,k
    integer :: wrf_i, wrf_j, wrf_k
    real(kind=dp) :: pressure
    real(kind=dp) :: potential_temp,temperature
    real(kind=dp) :: rh

    call wrf_message('wrf_to_partmc_init: starting initialization of &
         pressure, temperature and rh')

    do i = pmc_is, pmc_ie
    do k = pmc_ks, pmc_ke
    do j = pmc_js, pmc_je
       wrf_i = env_states(i,k,j)%ix
       wrf_j = env_states(i,k,j)%iy
       wrf_k = env_states(i,k,j)%iz
       ! Calculate the pressure
       pressure = grid%p(wrf_i,wrf_k,wrf_j)+grid%pb(wrf_i,wrf_k,wrf_j)
       ! Calculate the temperature
       potential_temp = grid%t_2(wrf_i,wrf_k,wrf_j)+t0
       temperature  = potential_temp*(pressure/p0)**(rcp)
       ! Set the value in the env_state
       env_states(i,k,j)%pressure = pressure
       env_states(i,k,j)%temp = temperature
       env_states(i,k,j)%rel_humid = rel_hum( &
            real(grid%moist(wrf_i,wrf_k,wrf_j,p_qv),dp),pressure,temperature)
    end do
    end do
    end do

    call wrf_message('wrf_to_partmc_init: finished initialization of &
         pressure, temperature and rh')

  end subroutine init_wrf_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Cleans up all the PartMC structures
  subroutine partmc_cleanup(scenarios, env_states, aero_data, aero_states, &
       gas_data,gas_states, pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke)

    !> PartMC east-west start of domain.
    integer, intent(in) :: pmc_is
    !> PartMC east-west end of domain.
    integer, intent(in) :: pmc_ie
    !> PartMC north-south start of domain.
    integer, intent(in) :: pmc_js
    !> PartMC north-south end of domain.
    integer, intent(in) :: pmc_je
    !> PartMC top-bottom start of domain.
    integer, intent(in) :: pmc_ks
    !> PartMC top-bottom end of domain.
    integer, intent(in) :: pmc_ke
    !> Scenario data.
    type(scenario_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenarios
    !> Environmental states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: aero_states
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data 
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: gas_states

    ! Nothing to clean up anymore
    call pmc_mpi_barrier()

  end subroutine partmc_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read in gas and aerosol emissions for a grid cell.
  subroutine init_read_in_emissions(scenario, i, j, k, aero_data, gas_data, &
       prefix)

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
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data
    !> File prefix.
    character(len=*), intent(in) :: prefix

    character(len=INPUT_FILE_PATH_NAME_LEN) :: file
    character(len=AERO_MODE_NAME_LEN) :: mode_name, weight_class
    integer :: n_time, ncid
    character(len=NF90_MAX_NAME) :: check_name
    integer :: status, check_dim_size
    integer :: dimid_n_times, dimid_n_modes, dimid_n_aero_specs, &
         dimid_n_gas_specs
    integer :: i_loop
    integer :: n_modes, n_aero_specs, n_gas_specs

    real(kind=dp), allocatable, dimension(:) :: char_radius, std
    integer, allocatable, dimension(:) :: dist_type, source
    real(kind=dp), allocatable, dimension(:,:) :: num_conc, gas_mix_rats
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

    ! FIXME: Add error checking to assure that n_spec in the input file
    ! is equal to the n_spec in aero_data

    write(file, '(a,a,i3.3,a,i3.3,a,i3.3,a)') &
         trim(prefix),'_',i,'_',j,'_',k,'.nc'

    !write(*,'(a,a)') 'reading emissions file ', trim(file)

    call pmc_nc_open_read(file, ncid)

    ! Find the number of times
    status = nf90_inq_dimid(ncid, "n_times", dimid_n_times)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_times, check_name, &
            check_dim_size))
    n_time = check_dim_size

    if (allocated(scenario%aero_emission)) deallocate(scenario%aero_emission)
    if (allocated(scenario%aero_emission_time)) &
         deallocate(scenario%aero_emission_time)
    if (allocated(scenario%aero_emission_rate_scale)) &
         deallocate(scenario%aero_emission_rate_scale)
    if (allocated(scenario%gas_emission_time)) &
         deallocate(scenario%gas_emission_time)
    if (allocated(scenario%gas_emission_rate_scale)) &
         deallocate(scenario%gas_emission_rate_scale)
    if (allocated(scenario%gas_emission)) deallocate(scenario%gas_emission)

    allocate(scenario%aero_emission_time(n_time))
    allocate(scenario%aero_emission_rate_scale(n_time))
    allocate(scenario%aero_emission(n_time))
    allocate(scenario%gas_emission_time(n_time))
    allocate(scenario%gas_emission_rate_scale(n_time))
    allocate(scenario%gas_emission(n_time))

    ! We want to set the following
    !  - aero_emission_time (real)
    !  - aero_emission_rate (real)
    !  - aero_emission (aero_dist_t)

    call pmc_nc_read_real_1d(ncid, scenario%aero_emission_time, 'aero_emission_time', .true.)

    !
    status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, check_name, &
            check_dim_size))
    n_modes = check_dim_size
    status = nf90_inq_dimid(ncid, "n_aero_specs", dimid_n_aero_specs)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_aero_specs, & 
         check_name, check_dim_size))
    n_aero_specs = check_dim_size

!    status = nf90_inq_dimid(ncid, "n_gas_specs", dimid_n_gas_specs)
!    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_gas_specs, &
!         check_name, check_dim_size))
!    n_gas_specs = check_dim_size

    !print*, 'number of modes =', n_modes, 'number of specs = ', n_aero_specs, &
    !     'number of times = ', n_time 
 
    !character(len=AERO_MODE_NAME_LEN) :: name
    ! Mode type (given by module constants).
    !integer :: type
    ! Characteristic radius, with meaning dependent on mode type (m).
    !real(kind=dp) :: char_radius
    ! Log base 10 of geometric standard deviation of radius, (m).
    !real(kind=dp) :: log10_std_dev_radius
    ! Total number concentration of mode (#/m^3).
    !real(kind=dp) :: num_conc
    ! Species fractions by volume [length \c aero_data%%n_spec] (1).
    !real(kind=dp), pointer :: vol_frac(:)
    ! Source number.
    !integer :: source

    ! Values that are n_modes in size
    ! FIXME: We've currently just kept the number of modes constant for all
    ! time for a cell and for all cells. This is why there isn't a time
    ! dependence for each grid cell.
    !   - name
    !   - type
    !   - char_radius
    !   - log10_std_dev_radius
    !   - source

    ! FIXME: We've omitted name for now since characters are annoying.

    ! FIXME: currently missing from input, trivial to add
    allocate(dist_type(n_modes))
    !call pmc_nc_read_integer_1d(ncid, dist_type, 'type', .true.)

    allocate(char_radius(n_modes))
    call pmc_nc_read_real_1d(ncid, char_radius, 'char_radius', .true.)

    allocate(std(n_modes))
    call pmc_nc_read_real_1d(ncid, std, 'log10_std_dev_radius', .true.)

    allocate(source(n_modes))
    call pmc_nc_read_integer_1d(ncid, source, 'source_id', .true.)

    ! Values that are n_modes and n_times in size
    !   - number concentration
    allocate(num_conc(n_modes, n_time))
    call pmc_nc_read_real_2d(ncid, num_conc, 'num_conc', .true.)

    ! Values that are n_modes, n_times and n_specs in size
    !   - vol_frac
    allocate(vol_frac(n_aero_specs,n_modes,n_time))
    call pmc_nc_read_real_3d(ncid, vol_frac, 'vol_frac', .true.)

    !FIXME: Read in
    allocate(vol_frac_std(n_aero_specs,n_modes,n_time))
    vol_frac_std = 0.0d0

    !FIXME: read in gas emissions
    status = nf90_inq_dimid(ncid, "n_gas_specs", dimid_n_gas_specs)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_gas_specs, &
            check_name, check_dim_size))
    n_gas_specs = check_dim_size
    allocate(gas_mix_rats(n_gas_specs, n_time))

    call pmc_nc_read_real_1d(ncid, scenario%gas_emission_time, &
         'gas_emission_time', .true.)
    call pmc_nc_read_real_2D(ncid, gas_mix_rats, 'gas_emission', .true.)

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

    ! Loop to store the data
    do i_time = 1, n_time
    ! FIXME: This doesn't need to be hardcoded - read from a NetCDF file
    scenario%aero_emission_rate_scale(i_time) = 1.0d0
    scenario%gas_emission_rate_scale(i_time) = 1.0d0
    if (allocated(scenario%aero_emission(i_time)%mode)) &
         deallocate(scenario%aero_emission(i_time)%mode) 
    allocate(scenario%aero_emission(i_time)%mode(n_modes))

    do i_mode = 1, n_modes
       ! FIXME: This should be from the netcdf file as well but strings are
       ! somewhat difficult
       !write(mode_name,'(a,i2.2)') 'emit_mode_', i_mode
       scenario%aero_emission(i_time)%mode(i_mode)%name =  &
            trim(source_name(source(i_mode)))
       ! Check to see if it is in aero_data
       dummy = aero_data_source_by_name(aero_data, &
            source_name(source(i_mode)))
       ! FIXME: This should be from the netcdf file but we are always doing
       ! log-normal at least for now
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
       scenario%aero_emission(i_time)%mode(i_mode)%source = dummy !source(i_mode)
       !
       if (char_radius(i_mode) < 0.01e-6) then
          write(weight_class,'(I3)') emission_weight_classes(i_mode)
       else
          write(weight_class,'(I3,a)') emission_weight_classes(i_mode), 'accumulation'
       end if
       scenario%aero_emission(i_time)%mode(i_mode)%weight_class = &
            aero_data_weight_class_by_name(aero_data, trim(weight_class))
    end do
    ! Gases
    call gas_state_set_size(scenario%gas_emission(i_time), &
         gas_data_n_spec(gas_data))
    scenario%gas_emission(i_time)%mix_rat = gas_mix_rats(:,i_time)
    end do

    call pmc_nc_check(nf90_close(ncid))

    ! deallocate
    deallocate(dist_type)
    deallocate(char_radius)
    deallocate(std)
    deallocate(source)
    deallocate(num_conc)
    deallocate(vol_frac)
    deallocate(vol_frac_std)

  end subroutine init_read_in_emissions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read in aerosol and gas boundary conditions.
  subroutine init_read_in_bcs(grid, scenario, i, j, k, aero_data, aero_state, &
       gas_state, env_state, prefix, restart)

    !> WRF domain.
    type(domain), intent(inout) :: grid
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
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environmental state.
    type(env_state_t), intent(in) :: env_state
    !> File prefix.
    character(len=*), intent(in) :: prefix
    !> If this is a restart, we don't sample
    logical, intent(in) :: restart

    character(len=INPUT_FILE_PATH_NAME_LEN) :: file
    character(len=AERO_MODE_NAME_LEN) :: mode_name, weight_class
    integer :: ncid
    character(len=NF90_MAX_NAME) :: check_name
    integer :: status, check_dim_size
    integer :: dimid_n_times, dimid_n_modes
    integer :: dimid_n_aero_specs, dimid_n_gas_specs
    integer :: i_loop
    integer :: n_times, n_modes, n_aero_specs, n_gas_specs

    real(kind=dp), allocatable, dimension(:) :: char_radius, std
    integer, allocatable, dimension(:) :: dist_type, source
    real(kind=dp), allocatable, dimension(:,:) :: num_conc, gas_mix_rats
    real(kind=dp), allocatable, dimension(:,:,:) :: vol_frac, vol_frac_std
    integer :: i_time, i_mode
    integer :: dummy, counter
    integer :: n_part_added
    integer :: n_extra_modes
    real(kind=dp) :: density
    logical, parameter :: extra_modes = .false.

    n_extra_modes = 0
    if (extra_modes) n_extra_modes = 2

    write(file, '(a,a,i3.3,a,i3.3,a,i3.3,a)') &
         trim(prefix),'_',i,'_',j,'_',k,'.nc'

    !write(*,'(a,a)') 'reading BC file ', trim(file)

    call pmc_nc_open_read(file, ncid)

    ! Find the number of times
    status = nf90_inq_dimid(ncid, "n_times", dimid_n_times)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_times, check_name, &
            check_dim_size))
    n_times = check_dim_size

    if (allocated(scenario%aero_background)) &
         deallocate(scenario%aero_background)
    if (allocated(scenario%aero_dilution_time)) &
         deallocate(scenario%aero_dilution_time)
    if (allocated(scenario%aero_dilution_rate)) &
         deallocate(scenario%aero_dilution_rate)
    if (allocated(scenario%gas_background)) &
         deallocate(scenario%gas_background)
    if (allocated(scenario%gas_dilution_time)) &
         deallocate(scenario%gas_dilution_time)
    if (allocated(scenario%gas_dilution_rate)) &
         deallocate(scenario%gas_dilution_rate)

    allocate(scenario%aero_background(n_times))
    allocate(scenario%aero_dilution_time(n_times))
    allocate(scenario%aero_dilution_rate(n_times))
    allocate(scenario%gas_background(n_times))
    allocate(scenario%gas_dilution_time(n_times))
    allocate(scenario%gas_dilution_rate(n_times))

    call pmc_nc_read_real_1d(ncid, scenario%aero_dilution_time, &
         'aero_bc_time', .true.)

    status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, check_name, &
            check_dim_size))
    n_modes = check_dim_size
    status = nf90_inq_dimid(ncid, "n_aero_specs", dimid_n_aero_specs)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_aero_specs, &
            check_name, check_dim_size))
    n_aero_specs = check_dim_size

    allocate(dist_type(n_modes))
    !call pmc_nc_read_integer_1d(ncid, dist_type, 'type', .true.)
    allocate(char_radius(n_modes))
    call pmc_nc_read_real_1d(ncid, char_radius, 'char_radius', .true.)
    allocate(std(n_modes))
    call pmc_nc_read_real_1d(ncid, std, 'log10_std_dev_radius', .true.)

    allocate(source(n_modes))
    call pmc_nc_read_integer_1d(ncid, source, 'source', .true.)

    allocate(num_conc(n_modes, n_times))
    call pmc_nc_read_real_2d(ncid, num_conc, 'num_conc', .true.)

    allocate(vol_frac(n_aero_specs,n_modes,n_times))
    call pmc_nc_read_real_3d(ncid, vol_frac, 'vol_frac', .true.)

    ! FIXME:
    allocate(vol_frac_std(n_aero_specs,n_modes,n_times))
    vol_frac_std = 0.0d0

    ! Loop to store the data
    do i_time = 1,n_times
       ! FIXME: This doesn't need to be hardcoded - read from a NetCDF file
        scenario%aero_dilution_rate(i_time) = 1.0d0
        allocate(scenario%aero_background(i_time)%mode(n_modes+n_extra_modes))
        do i_mode = 1, n_modes
           ! FIXME: Disable special boundary condition mode
           !write(mode_name,'(a,i2.2)') 'bc_mode_', i_mode
           write(mode_name,'(a,i2.2)') 'bc_mode_', i_mode
           scenario%aero_background(i_time)%mode(i_mode)%name = trim(mode_name)
           ! Check to see if it is in aero_data
           dummy = aero_data_source_by_name(aero_data, mode_name)
           scenario%aero_background(i_time)%mode(i_mode)%type = &
                AERO_MODE_TYPE_LOG_NORMAL
           scenario%aero_background(i_time)%mode(i_mode)%char_radius = &
                char_radius(i_mode)
           scenario%aero_background(i_time)%mode(i_mode)%log10_std_dev_radius = &
                std(i_mode)
           scenario%aero_background(i_time)%mode(i_mode)%num_conc = &
                num_conc(i_mode, i_time)
           scenario%aero_background(i_time)%mode(i_mode)%vol_frac = &
                vol_frac(1:20,i_mode,i_time)
           scenario%aero_background(i_time)%mode(i_mode)%vol_frac_std = &
                vol_frac_std(1:20,i_mode,i_time)
           ! FIXME: Fix this source information
           scenario%aero_background(i_time)%mode(i_mode)%source = dummy !source(i_mode)
           write(weight_class,'(a,i3.3)') "BC", i_mode
           scenario%aero_background(i_time)%mode(i_mode)%weight_class = &
                aero_data_weight_class_by_name(aero_data, trim(weight_class))
       end do

       ! Hack to add some more modes for now
       if (extra_modes) then
       counter = 1
       do i_mode = n_modes+1, n_modes+n_extra_modes
       write(mode_name,'(a,i2.2)') 'ic_mode_new', counter
       scenario%aero_background(i_time)%mode(i_mode)%name = trim(mode_name)
       ! Check to see if it is in aero_data
       dummy = aero_data_source_by_name(aero_data, mode_name)
       scenario%aero_background(i_time)%mode(i_mode)%type = &
            AERO_MODE_TYPE_LOG_NORMAL
       scenario%aero_background(i_time)%mode(i_mode)%char_radius = 2.0d-8
       scenario%aero_background(i_time)%mode(i_mode)%log10_std_dev_radius = .3d0
       scenario%aero_background(i_time)%mode(i_mode)%num_conc = 0.0d0
       if ((counter == 1 .and. j < 9) .or.  &
          (counter == 2 .and. j >= 9))  then
          scenario%aero_background(i_time)%mode(i_mode)%num_conc = 5.0d9 / &
               grid%alt(i,k,j)
       endif
       scenario%aero_background(i_time)%mode(i_mode)%vol_frac = &
            vol_frac(1:20,1,i_time)
       scenario%aero_background(i_time)%mode(i_mode)%source = dummy !source(i_mode)
       scenario%aero_background(i_time)%mode(i_mode)%vol_frac_std = &
           vol_frac_std(1:20,1,i_time)

       write(weight_class,'(a,i3.3)') "IC_new", counter
       scenario%aero_background(i_time)%mode(i_mode)%weight_class = &
            aero_data_weight_class_by_name(aero_data, weight_class)
       counter = counter + 1
       end do
       end if

       ! Gases
!       scenario%gas_dilution_rate(i_time) = 1.0d0
!       call gas_state_set_size(scenario%gas_background(i_time), &
!            n_gas_specs)
!       scenario%gas_background(i_time)%mix_rat = gas_mix_rats(:,i_time)
    end do

    call pmc_nc_check(nf90_close(ncid))

    ! FIXME:
    ! Sample particles with first boundary condition
    density = 1.0d0 / grid%alt(i,k,j) !env_state_air_den(env_state)

    if (.not. restart) then
       call aero_state_add_aero_dist_sample(aero_state, aero_data, &
            scenario%aero_background(1), density, 1.0d0, 0d0, .true., .true., &
            n_part_added)
    end if

    ! deallocate
    deallocate(dist_type)
    deallocate(char_radius)
    deallocate(std)
    deallocate(source)
    deallocate(num_conc)
    deallocate(vol_frac)
    deallocate(vol_frac_std)
!    deallocate(gas_mix_rats)

  end subroutine init_read_in_bcs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reads in aerosol and gas initial conditions for a grid cell.
  subroutine init_read_in_ics(grid, aero_state, gas_state, i, j, nz, aero_data, &
       gas_data, env_state, prefix)

    integer,intent(in) :: nz 
    !> WRF domain.
    type(domain), intent(in) :: grid
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state(nz)
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state(nz)
    !> East-west index of grid cell.
    integer,intent(in) :: i
    !> North-south index of grid cell.
    integer,intent(in) :: j
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data
    !> Environmental states.
    type(env_state_t), intent(in) :: env_state(nz)
    !> Filename prefix.
    character(len=*), intent(in) :: prefix

    character(len=INPUT_FILE_PATH_NAME_LEN) :: file, group
    character(len=AERO_MODE_NAME_LEN) :: mode_name
    character(len=AERO_SOURCE_NAME_LEN) :: weight_class
    integer :: n_time, ncid, ncid_group
    character(len=NF90_MAX_NAME) :: check_name
    integer :: status, check_dim_size
    integer :: dimid_n_times, dimid_n_modes
    integer :: dimid_n_aero_specs, dimid_n_gas_specs
    integer :: i_loop
    integer :: n_modes, n_aero_specs, n_gas_specs
    type(aero_dist_t) :: aero_dist_init
    real(kind=dp), allocatable, dimension(:) :: char_radius, std
    integer, allocatable, dimension(:) :: dist_type, source
    real(kind=dp), allocatable, dimension(:) :: num_conc
    real(kind=dp), allocatable, dimension(:,:) :: vol_frac, vol_frac_std
    integer :: i_time, i_mode, k
    integer :: dummy
    integer :: n_part_added
    real(kind=dp) :: density

    write(file, '(a,a,i3.3,a,i3.3,a)') trim(prefix),'_',i,'_',j,'.nc'
    call pmc_nc_open_read(file, ncid)

    do k = 1,nz
       write(group, '(a,i3.3)') 'level_', k
       status = nf90_inq_ncid(ncid, group, ncid_group)

       ! Find the number of modes and species from netcdf file
       status = nf90_inq_dimid(ncid_group, "n_modes", dimid_n_modes)
       call pmc_nc_check(nf90_Inquire_Dimension(ncid_group, dimid_n_modes, check_name, &
            check_dim_size))
       n_modes = check_dim_size
       status = nf90_inq_dimid(ncid_group, "n_aero_specs", dimid_n_aero_specs)
       call pmc_nc_check(nf90_Inquire_Dimension(ncid_group, dimid_n_aero_specs, &
            check_name, check_dim_size))
       n_aero_specs = check_dim_size

       allocate(dist_type(n_modes))
       !call pmc_nc_read_integer_1d(ncid_group, dist_type, 'type', .true.)
       allocate(char_radius(n_modes))
       call pmc_nc_read_real_1d(ncid_group, char_radius, 'char_radius', .true.)
       allocate(std(n_modes))
       call pmc_nc_read_real_1d(ncid_group, std, 'log10_std_dev_radius', .true.)

       allocate(source(n_modes))
       call pmc_nc_read_integer_1d(ncid_group, source, 'source', .true.)

       allocate(num_conc(n_modes))
       call pmc_nc_read_real_1d(ncid_group, num_conc, 'num_conc', .true.)

       allocate(vol_frac(n_aero_specs,n_modes))
       call pmc_nc_read_real_2d(ncid_group, vol_frac, 'vol_frac', .true.)    

       ! FIXME: Read this in?
       allocate(vol_frac_std(n_aero_specs,n_modes))
       vol_frac_std = 0.0d0
       allocate(aero_dist_init%mode(n_modes))
       ! Recreate the aerosol distribution to sample
       do i_mode = 1, n_modes
          write(mode_name,'(a,i2.2)') 'ic_mode_', i_mode
          aero_dist_init%mode(i_mode)%name = trim(mode_name)
          ! Check to see if it is in aero_data
          dummy = aero_data_source_by_name(aero_data, mode_name)
          aero_dist_init%mode(i_mode)%type = &
               AERO_MODE_TYPE_LOG_NORMAL
          aero_dist_init%mode(i_mode)%char_radius = &
               char_radius(i_mode)
          aero_dist_init%mode(i_mode)%log10_std_dev_radius = &
               std(i_mode)
          aero_dist_init%mode(i_mode)%num_conc = &
               num_conc(i_mode)
          aero_dist_init%mode(i_mode)%vol_frac = &
               vol_frac(1:n_aero_specs,i_mode)
          ! FIXME: This should more flexible
          aero_dist_init%mode(i_mode)%source = dummy !source(i_mode)
          ! FIXME: we don't have this information
          aero_dist_init%mode(i_mode)%vol_frac_std = vol_frac_std(:,i_mode)

          write(weight_class,'(a,i3.3)') "IC", i_mode

          aero_dist_init%mode(i_mode)%weight_class = & 
               aero_data_weight_class_by_name(aero_data, trim(weight_class))
       end do

       ! FIXME: Halving and doubling from grid
       density = 1.0d0 / grid%alt(i,k,j) !env_state_air_den(env_state)
       call aero_state_add_aero_dist_sample(aero_state(k), aero_data, &
            aero_dist_init, density, 1.0d0, 0d0, .true., .true., &
            n_part_added)

       deallocate(dist_type)
       deallocate(char_radius)
       deallocate(std)
       deallocate(source)
       deallocate(num_conc)
       deallocate(vol_frac)
       deallocate(vol_frac_std)
       deallocate(aero_dist_init%mode)
    end do

    call pmc_nc_check(nf90_close(ncid))

  end subroutine init_read_in_ics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_read_in_restart(aero_states, gas_states, i, j, nz, &
       aero_data, gas_data, partmc_restart_prefix, partmc_restart_index)

    !> Aerosol state to reconstruct.
    type(aero_state_t), dimension(nz), intent(inout) :: aero_states
    !> Gas state to reconstruct.
    type(gas_state_t), dimension(nz), intent(inout) :: gas_states
    !> East-west grid cell index.
    integer, intent(in) :: i
    !> North-south grid cell index.
    integer, intent(in) :: j
    !> Number of grid cells in column.
    integer, intent(in) :: nz
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Prefix for restart files.
    character(len=*), intent(in) :: partmc_restart_prefix
    !> Time step index of restart files.
    integer, intent(in) :: partmc_restart_index

    character(LEN=INPUT_FILE_PATH_NAME_LEN) :: filename
    integer :: ncid, status
    integer :: k, i_part, i_part_g, i_comp_g
    integer :: i_group, i_class, n_group, n_class
    integer :: total_particles, total_components
    real(kind=dp), allocatable :: gas_mixing_ratios(:,:), &
         aero_particle_mass(:,:)
    integer :: dimid_gas_species, dimid_num_levels, dimid_aero_particle
    integer :: dimid_aero_components, dimid_aero_species
    integer, allocatable :: n_parts(:)
    integer, allocatable :: part_start_index(:)
    real(kind=dp), allocatable :: weight_exponent(:,:,:), &
         weight_magnitude(:,:,:)
    integer, allocatable :: weight_type(:,:,:)
    integer, allocatable :: aero_particle_weight_group(:)
    integer, allocatable :: aero_particle_weight_class(:)
    real(kind=dp), allocatable :: aero_absorb_cross_sect(:,:)
    real(kind=dp), allocatable :: aero_scatter_cross_sect(:,:)
    real(kind=dp), allocatable :: aero_asymmetry(:,:)
    real(kind=dp), allocatable :: aero_refract_shell_real(:,:)
    real(kind=dp), allocatable :: aero_refract_shell_imag(:,:)
    real(kind=dp), allocatable :: aero_refract_core_real(:,:)
    real(kind=dp), allocatable :: aero_refract_core_imag(:,:)
    real(kind=dp), allocatable :: aero_core_vol(:)
    integer, allocatable :: aero_water_hyst_leg(:)
    real(kind=dp), allocatable :: aero_num_conc(:)
    integer(kind=8), allocatable :: aero_id(:)
    real(kind=dp), allocatable :: aero_least_create_time(:)
    real(kind=dp), allocatable :: aero_greatest_create_time(:)
    integer, allocatable :: aero_component_particle_num(:)
    integer, allocatable :: aero_component_source_num(:)
    integer, allocatable :: aero_component_len(:)
    integer, allocatable :: aero_component_start_ind(:)
    real(kind=dp), allocatable :: aero_component_create_time(:)
    integer :: array_position, i_comp, n_part
    type(aero_particle_t) :: aero_particle

    write(filename, '(a,a,i3.3,a,i3.3,a,i8.8,a)') &
         trim(partmc_restart_prefix),'_',i,'_',j,'_', &
         partmc_restart_index, '.nc'

    call pmc_nc_open_read(filename, ncid)

    call pmc_nc_read_real_2d(ncid, gas_mixing_ratios, &
         "gas_mixing_ratio")

    call pmc_nc_read_integer_1d(ncid, n_parts, "n_parts")
    call pmc_nc_read_integer_1d(ncid, part_start_index, &
         "part_start_index")
    call pmc_nc_read_real_2d(ncid, aero_particle_mass, &
         "aero_particle_mass")
    call pmc_nc_read_integer_1d(ncid, aero_component_len, &
         "aero_component_len") 
    call pmc_nc_read_integer_1d(ncid, aero_component_start_ind, &
         "aero_component_start_ind")
    call pmc_nc_read_integer_1d(ncid, aero_component_particle_num, &
         "aero_component_particle_num")
    call pmc_nc_read_integer_1d(ncid, aero_component_source_num, &
         "aero_component_source_num")
    call pmc_nc_read_real_1d(ncid, aero_component_create_time, &
         "aero_component_create_time")
    call pmc_nc_read_real_1d(ncid, aero_least_create_time, &
         "aero_least_create_time")
    call pmc_nc_read_real_1d(ncid, aero_greatest_create_time, &
         "aero_greatest_create_time")
    call pmc_nc_read_integer_1d(ncid, aero_particle_weight_group, &
         "aero_particle_weight_group")
    call pmc_nc_read_integer_1d(ncid, aero_particle_weight_class, &
         "aero_particle_weight_class")
    call pmc_nc_read_integer_1d(ncid, aero_water_hyst_leg, &
         "aero_water_hyst_leg")
    call pmc_nc_read_real_1d(ncid, aero_num_conc, &
         "aero_num_conc")
    call pmc_nc_read_integer64_1d(ncid, aero_id, "aero_id")

    ! Optional optical
    call pmc_nc_read_real_2d(ncid, aero_absorb_cross_sect, &
         "aero_absorb_cross_sect", must_be_present=.false.)
    call pmc_nc_read_real_2d(ncid, aero_scatter_cross_sect, &
         "aero_scatter_cross_sect", must_be_present=.false.)
    call pmc_nc_read_real_2d(ncid, aero_asymmetry, &
         "aero_asymmetry", must_be_present=.false.)
    call pmc_nc_read_real_2d(ncid, aero_refract_shell_real, &
         "aero_refract_shell_real", must_be_present=.false.)
    call pmc_nc_read_real_2d(ncid, aero_refract_shell_imag, &
         "aero_refract_shell_imag", must_be_present=.false.)
    call pmc_nc_read_real_2d(ncid, aero_refract_core_real, &
         "aero_refract_core_real", must_be_present=.false.)
    call pmc_nc_read_real_2d(ncid, aero_refract_core_imag, &
         "aero_refract_core_imag", must_be_present=.false.)
    call pmc_nc_read_real_1d(ncid, aero_core_vol, &
         "aero_core_vol", must_be_present=.false.)

    do k = 1,nz
       gas_states(k)%mix_rat = gas_mixing_ratios(k,:)
    end do

    call pmc_nc_read_integer_3d(ncid, weight_type, &
         "weight_type")
    call pmc_nc_read_real_3d(ncid, weight_magnitude, &
         "weight_magnitude")
    call pmc_nc_read_real_3d(ncid, weight_exponent, &
         "weight_exponent")

    n_group = size(weight_magnitude,2)
    n_class = size(weight_magnitude,3)
    do k = 1,nz
    call aero_weight_array_set_sizes(aero_states(k)%awa, n_group, n_class)
    do i_group = 1,n_group
    do i_class = 1,n_class
       aero_states(k)%awa%weight(i_group,i_class)%type = &
            weight_type(k,i_group,i_class)
       aero_states(k)%awa%weight(i_group,i_class)%magnitude = &
            weight_magnitude(k,i_group,i_class)
       aero_states(k)%awa%weight(i_group,i_class)%exponent = &
            weight_exponent(k,i_group,i_class)
    end do
    end do
    end do

    i_comp_g = 1
    n_part = size(aero_num_conc)
    do k = 1,nz
    do i_part = part_start_index(k),part_start_index(k) + n_parts(k) - 1
       call aero_particle_zero(aero_particle, aero_data)

       aero_particle%vol = aero_particle_mass(i_part, :) / aero_data%density

       if (allocated(aero_particle%component)) &
            deallocate(aero_particle%component)
       allocate(aero_particle%component(aero_component_len(i_part)))
       do i_comp = 1, aero_component_len(i_part)
            aero_particle%component(i_comp)%source_id = &
                 aero_component_source_num(i_comp_g)
            aero_particle%component(i_comp)%create_time = &
                 aero_component_create_time(i_comp_g)
            i_comp_g = i_comp_g + 1
       end do
       aero_particle%weight_group = aero_particle_weight_group(i_part)
       aero_particle%weight_class = aero_particle_weight_class(i_part)
       if (size(aero_absorb_cross_sect,1) == n_part) then
          aero_particle%absorb_cross_sect = aero_absorb_cross_sect(i_part,:)
       end if
       if (size(aero_scatter_cross_sect,1) == n_part) then
          aero_particle%scatter_cross_sect = aero_scatter_cross_sect(i_part,:)
       end if
       if (size(aero_asymmetry,1) == n_part) then
            aero_particle%asymmetry = aero_asymmetry(i_part,:)
       end if
       if ((size(aero_refract_shell_real,1) == n_part) &
            .and. (size(aero_refract_shell_imag) == n_part)) then
          aero_particle%refract_shell = &
               cmplx(aero_refract_shell_real(i_part,:), &
               aero_refract_shell_imag(i_part,:), kind=dc)
       end if
       if ((size(aero_refract_core_real,1) == n_part) &
            .and. (size(aero_refract_core_imag) == n_part)) then
          aero_particle%refract_core = cmplx(aero_refract_core_real(i_part,:), &
               aero_refract_core_imag(i_part,:), kind=dc)
       end if
       if (size(aero_core_vol) == n_part) then
          aero_particle%core_vol = aero_core_vol(i_part)
       end if
       aero_particle%water_hyst_leg = aero_water_hyst_leg(i_part)
       aero_particle%id = aero_id(i_part)
       aero_particle%least_create_time = aero_least_create_time(i_part)
       aero_particle%greatest_create_time = aero_greatest_create_time(i_part)
       call aero_state_add_particle(aero_states(k), aero_particle, aero_data)
    end do
    end do

    do k = 1,nz
       call aero_info_array_zero(aero_states(k)%aero_info_array)
       call aero_state_sort(aero_states(k), aero_data)
    end do

    call pmc_nc_read_integer64(ncid, next_id, "next_id")

    call pmc_nc_close(ncid)

  end subroutine init_read_in_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the number of sources for aerosols from input files.
  subroutine get_sources_and_weights(aero_data, prefix_ics, &
       prefix_bcs, prefix_emissions, periodic_bcs)

    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> File prefix of initial conditions.
    character(len=*), intent(in) :: prefix_ics
    !> File prefix of boundary conditions.
    character(len=*), intent(in) :: prefix_bcs
    !> File prefix of emissions.
    character(len=*), intent(in) :: prefix_emissions
    !> Whether or not boundary conditions are periodic.
    logical, intent(in) :: periodic_bcs 

    integer :: ncid, ncid_group
    integer :: status
    integer :: dimid_n_modes, check_dim_size
    character(len=NF90_MAX_NAME) :: check_name, name
    integer :: i, j, k
    character(len=INPUT_FILE_PATH_NAME_LEN) :: file, group
    integer, parameter :: MAX_SOURCES = 1000
    character(len=((AERO_SOURCE_NAME_LEN + 2) * MAX_SOURCES)) &
         :: aero_source_names
    character(len=AERO_SOURCE_NAME_LEN), allocatable :: source_name(:)

    integer :: dimid_aero_source, n_source, varid_aero_source, i_source, i_str
    integer :: dummy, i_ic, i_bc, i_emit
    integer :: i_mode
    character(len=AERO_MODE_NAME_LEN) :: mode_name
    character(len=AERO_SOURCE_NAME_LEN) :: weight_class
    character(len=AERO_SOURCE_NAME_LEN), allocatable :: weight_name(:)
    integer :: num_ic_modes, num_bc_modes, num_emit_modes
    integer, allocatable, dimension(:) :: emission_weight_classes
    integer :: i_ss

    i = 1
    j = 1
    k = 1

    write(file, '(a,a,i3.3,a,i3.3,a)') &
         trim(prefix_ics),'_',i,'_',j,'.nc'

    call pmc_nc_open_read(file, ncid)

    write(group, '(a,i3.3)') 'level_', k
    status = nf90_inq_ncid(ncid, group, ncid_group)
    status = nf90_inq_dimid(ncid_group, "n_modes", dimid_n_modes)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid_group, dimid_n_modes, &
         check_name, check_dim_size))
    num_ic_modes = check_dim_size
    ! FIXME: We should read in weight classes for ICs
    do i_mode = 1,num_ic_modes
    write(weight_class,'(a,i3.3)') "IC", i_mode
    dummy = aero_data_weight_class_by_name(aero_data, weight_class)
    end do
    call pmc_nc_check(nf90_close(ncid))

    if (.not. periodic_bcs) then
    write(file, '(a,a,i3.3,a,i3.3,a,i3.3,a)') &
         trim(prefix_bcs),'_',i,'_',j,'_',k,'.nc'

    call pmc_nc_open_read(file, ncid)
    status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, &
         check_name, check_dim_size))
    num_bc_modes = check_dim_size
    ! FIXME: We should read in weight classes for BCs
    do i_mode = 1,num_bc_modes
    write(weight_class,'(a,i3.3)') "BC", i_mode
    dummy = aero_data_weight_class_by_name(aero_data, weight_class)
    end do
    call pmc_nc_check(nf90_close(ncid))
    end if

    write(file, '(a,a,i3.3,a,i3.3,a,i3.3,a)') &
         trim(prefix_emissions),'_',i,'_',j,'_',k,'.nc'

    call pmc_nc_open_read(file, ncid)
    status = nf90_inq_dimid(ncid, "n_modes", dimid_n_modes)
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_n_modes, &
         check_name, check_dim_size))
    num_emit_modes = check_dim_size

    ! Read in strings
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

    ! FIXME: This should be strings?
    call pmc_nc_read_integer_1d(ncid, emission_weight_classes, &
         "source_weight_class", .true.)
    call pmc_nc_check(nf90_close(ncid))
    call ensure_string_array_size(weight_name, num_emit_modes)
    do i_emit = 1,num_emit_modes
        write(weight_class,'(I3)') emission_weight_classes(i_emit)
        weight_name(i_emit) = trim(weight_class)
    end do

    do i_ic = 1, num_ic_modes
       write(mode_name,'(a,i2.2)') 'ic_mode_', i_ic
       dummy = aero_data_source_by_name(aero_data, mode_name)
    end do
    do i_source= 1, n_source
       dummy = aero_data_source_by_name(aero_data, source_name(i_source))
    end do
    do i_emit = 1,num_emit_modes
       dummy = aero_data_weight_class_by_name(aero_data, weight_name(i_emit))
    end do

    ! Hack to add another weight class for now
    do i_emit = 1,num_emit_modes
       write(weight_class,'(I3,a)') emission_weight_classes(i_emit), 'accumulation'
       dummy = aero_data_weight_class_by_name(aero_data, weight_class)
    end do

    do i_bc = 1, num_bc_modes
       write(mode_name,'(a,i2.2)') 'bc_mode_', i_bc
       dummy = aero_data_source_by_name(aero_data, mode_name)
    end do

    ! Added sea salt
    do i_ss = 1,2
    write(mode_name,'(a,i2.2)') 'sea_salt_', i_ss
    dummy = aero_data_weight_class_by_name(aero_data, mode_name)
    dummy = aero_data_source_by_name(aero_data, mode_name)
    end do

  end subroutine get_sources_and_weights

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets emissions to be empty.
  subroutine init_zero_emissions(scenario, aero_data, gas_data)

    !> Scenario data.
    type(scenario_t), intent(inout) :: scenario
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data

    if (allocated(scenario%aero_emission_time)) deallocate(scenario%aero_emission_time)
    allocate(scenario%aero_emission_time(0))
    if (allocated(scenario%aero_emission_rate_scale)) deallocate(scenario%aero_emission_rate_scale)
    allocate(scenario%aero_emission_rate_scale(0))
    if (allocated(scenario%aero_emission)) deallocate(scenario%aero_emission)
    allocate(scenario%aero_emission(0))
    if (allocated(scenario%gas_emission_time)) deallocate(scenario%gas_emission_time)
    allocate(scenario%gas_emission_time(0))
    if (allocated(scenario%gas_emission_rate_scale)) deallocate(scenario%gas_emission_rate_scale)
    allocate(scenario%gas_emission_rate_scale(0))
    if (allocated(scenario%gas_emission)) deallocate(scenario%gas_emission)
    allocate(scenario%gas_emission(0))

    scenario%aero_emission_rate_scale = 0.0d0
    scenario%gas_emission_rate_scale = 0.0d0

  end subroutine init_zero_emissions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_partmc_to_wrf(grid, gas_states, pmc_is, pmc_ie, pmc_js, &
       pmc_je, pmc_ks, pmc_ke)

    !> PartMC east-west start of domain.
    integer, intent(in) :: pmc_is
    !> PartMC east-west end of domain.
    integer, intent(in) :: pmc_ie
    !> PartMC north-south start of domain.
    integer, intent(in) :: pmc_js
    !> PartMC north-south end of domain.
    integer, intent(in) :: pmc_je
    !> PartMC top-bottom start of domain.
    integer, intent(in) :: pmc_ks
    !> PartMC top-bottom end of domain.
    integer, intent(in) :: pmc_ke
    !> WRF domain.
    type(domain), intent(inout) :: grid
    !> Full domain of gas_states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(in) :: gas_states

    integer :: i, j, k, i_spec

    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       !FIXME: Should pass gas_data
       do i_spec = 1,77
       grid%chem(i,k,j,1+i_spec) = real(gas_states(i,k,j)%mix_rat(i_spec)) &
           / 1000.0d0
       end do
    end do
    end do
    end do

  end subroutine init_partmc_to_wrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes PartMC data structures on the WRF domain based on input files
  subroutine init_wrf_partmc_idealized(grid, scenario, env_states, &
       aero_data, aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, &
       pmc_je, pmc_ks, pmc_ke, nx, ny, nz, global_nx, global_ny, global_nz, &
       config_flags)

    !> PartMC east-west start of domain.
    integer, intent(in) :: pmc_is
    !> PartMC east-west end of domain.
    integer, intent(in) :: pmc_ie
    !> PartMC north-south start of domain.
    integer, intent(in) :: pmc_js
    !> PartMC north-south end of domain.
    integer, intent(in) :: pmc_je
    !> PartMC top-bottom start of domain.
    integer, intent(in) :: pmc_ks
    !> PartMC top-bottom end of domain.
    integer, intent(in) :: pmc_ke
    !> WRF domain.
    type(domain), intent(inout) :: grid
    !> Environment data.
    type(scenario_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenario
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout):: aero_states
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: gas_states
    !> Number of boxes in east-west direction.
    integer, intent(in) :: nx
    !> Number of boxes in the north-south direction.
    integer, intent(in) :: ny
    !> Number of boxes in the top-bottom direction.
    integer, intent(in) :: nz
    !> East-west domain dimension size of WRF domain.
    integer, intent(in) :: global_nx
    !> North-south domain dimension size of WRF domain.
    integer, intent(in) :: global_ny
    !> Top-bottom domain dimension size of WRF domain.
    integer, intent(in) :: global_nz
    !> WRF namelist configuration.
    type(grid_config_rec_type), intent(in) :: config_flags

    ! Local variables
    integer :: i, j, k
    integer :: p
    real(kind=dp) :: alpha
    ! UUID
    character(len=PMC_UUID_LEN) :: uuid
    ! Counting particles
    integer :: n_parts_added
    real(kind=dp) :: n_parts
    ! Local PartMC structures
    type(spec_file_t) :: sub_file
    character(LEN=INPUT_FILE_PATH_NAME_LEN) :: filename
    character(LEN=:), allocatable :: name
    ! Various structures to store inputs
    type(gas_data_t) :: input_gas_data
    type(aero_data_t) :: input_aero_data
    type(aero_data_t) :: global_aero_data
    character(len=AERO_MODE_NAME_LEN) :: mode_name, weight_class
    integer :: n_ic, n_emit, i_ic, i_emit
    integer :: n_modes

    real(kind=dp), allocatable, dimension(:) :: vol_frac, vol_frac_std
    type(aero_dist_t) :: aero_dist_init

    ! MPI variables
    character, allocatable :: buffer(:)
    character, allocatable :: aero_data_buffer(:)
    integer :: aero_data_size
    integer :: position
    integer :: max_buffer_size
    integer :: buffer_size
    integer :: wrf_i,wrf_j,wrf_k, i_mode, dummy
    real(kind=dp) :: time_ic, time_bc, t1, t2

    real(kind=dp) :: xpos, z
    real(kind=dp) :: xrad, yrad, r, r0, x0, y0, z0

    integer :: q, ix, iy, iz
    real(kind=dp) :: a,b,c,d,e,f,x_pos, y_pos, z_pos
    real(kind=dp) :: x_center, y_center
    integer :: rx, ry, rz
    real(kind=dp) :: q_value
    real(kind=dp) :: gauss_weights(5)
    real(kind=dp) :: gauss_points(5)
    integer :: pmc_o3
    integer :: n_class, i_spec
   
    integer, parameter :: n_gas_emit = 17 
    integer, dimension(n_gas_emit) :: gas_index
    character(len=5), dimension(n_gas_emit) :: gas_emit_names = &
        (/ "SO2  ", "NO2  ", "NO   ", "NH3  ", "CO   ", "ALD2 ", &
           "HCHO ", "ETH  ", "OLEI ", "OLET ", "TOL  ", "XYL  ", &
           "AONE ", "PAR  ", "ISOP ", "CH3OH", "ANOL "/) 
    real(kind=dp), dimension(n_gas_emit), parameter :: gas_emit_values = &
        (/4d-9, 3d-9, 6d-8, 9d-9, 8d-7, 2d-9, &
          4d-9, 2d-8, 6d-9, 6d-9, 6d-9, 6d-9, &
          8d-10, 2d-7, 2.4d-10, 2.4d-10, 5d-9 /)

    ! Initialize the random number generator
    call pmc_srand(grid%random_seed, pmc_mpi_rank())

    print*,'Initializing:',pmc_mpi_rank(),pmc_is,pmc_ie,pmc_js,pmc_je,pmc_ks,pmc_ke

    call wrf_message('In the partmc initialization subroutine')

    n_modes = 1

    ! Processor 0 reads in all the input data
    if (pmc_mpi_rank() == 0) then
       ! Load gas_data inputs
       call wrf_message('PartMC_init: Loading gas_data from input files')
       filename = 'gas_data.dat'
       call spec_file_open(filename, sub_file)
       call spec_file_read_gas_data(sub_file, input_gas_data)
       call spec_file_close(sub_file)
       ! Load aero_data inputs
       call wrf_message('PartMC_init: Loading aero_data from input files')
       filename = 'aero_data.dat'
       call spec_file_open(filename,sub_file)
       call spec_file_read_aero_data(sub_file, input_aero_data)
       call fractal_set_spherical(input_aero_data%fractal)
       call spec_file_close(sub_file)

       write(mode_name,'(a)') 'source_A'
       dummy = aero_data_weight_class_by_name(input_aero_data, mode_name)
       dummy = aero_data_source_by_name(input_aero_data, mode_name)
       ! Define a second mode
       if (n_modes > 1) then
       write(mode_name,'(a)') 'source_B'
       dummy = aero_data_weight_class_by_name(input_aero_data, mode_name)
       dummy = aero_data_source_by_name(input_aero_data, mode_name)
       end if

       ! Get the maximum pack size
       max_buffer_size = pmc_mpi_pack_size_gas_data(input_gas_data) &
            +  pmc_mpi_pack_size_aero_data(input_aero_data)
       allocate(buffer(max_buffer_size))
       ! Pack the buffer
       position = 0
       call pmc_mpi_pack_gas_data(buffer, position, input_gas_data)
       call pmc_mpi_pack_aero_data(buffer, position, input_aero_data)
       buffer_size = position
    end if

    ! Broadcast the buffer size
    call pmc_mpi_bcast_integer(max_buffer_size)

    ! Allocate the buffer on other nodes
    if (pmc_mpi_rank() /= 0) then
       allocate(buffer(max_buffer_size))
    end if

    ! Broadcast the buffer
    call pmc_mpi_bcast_packed(buffer)
    ! Unpack the buffer
    if (pmc_mpi_rank() /= 0) then
       position = 0
       call pmc_mpi_unpack_gas_data(buffer, position, input_gas_data)
       call pmc_mpi_unpack_aero_data(buffer, position, input_aero_data)
    end if
    deallocate(buffer)

    call wrf_message('PartMC_init: Setting initial information for env_states')

    do i = pmc_is, pmc_ie
    do k = pmc_ks, pmc_ke
    do j = pmc_js, pmc_je
        wrf_i = i
        wrf_j = j
        wrf_k = k
        env_states(i,k,j)%ix = wrf_i
        env_states(i,k,j)%iy = wrf_j
        env_states(i,k,j)%iz = wrf_k
        ! Set latitude/longitude
        env_states(i,k,j)%latitude = grid%xlat(wrf_i,wrf_j)
        env_states(i,k,j)%longitude = grid%xlong(wrf_i,wrf_j)
        ! Set an altitude/height in m
        ! WRF initially doesn't have z values so we must use geopotential here
        ! in a more creative way
        ! (Top - bottom)/g is the total height
        env_states(i,k,j)%height = (grid%phb(wrf_i,wrf_k+1,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k+1,wrf_j)-(grid%phb(wrf_i,wrf_k,wrf_j) + &
             grid%ph_1(wrf_i,wrf_k,wrf_j)))/9.8
        ! Find the bottom and add hlf the height for the center point
        ! (Top - bottom)/2g + (bottom)/g
        env_states(i,k,j)%altitude = 0.5d0*((grid%phb(wrf_i,wrf_k+1,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k+1,wrf_j)-(grid%phb(wrf_i,wrf_k,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k,wrf_j)))/9.8) + &
             ((grid%phb(wrf_i,wrf_k,wrf_j)+grid%ph_1(wrf_i,wrf_k,wrf_j))/9.8)
        env_states(i,k,j)%z_min = (grid%phb(wrf_i,wrf_k,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k,wrf_j))/9.8
        env_states(i,k,j)%z_max = (grid%phb(wrf_i,wrf_k+1,wrf_j)+ &
             grid%ph_1(wrf_i,wrf_k+1,wrf_j))/9.8
        ! FIXME: Density or inverse density
        env_states(i,k,j)%rrho = grid%alt(wrf_i,wrf_k,wrf_j)
        ! Set time and date information
        env_states(i,k,j)%start_day = grid%julian
        env_states(i,k,j)%start_time = grid%start_hour * 3600.0d0 + &
             grid%start_minute * 60.0d0 + grid%start_second
    end do
    end do
    end do

    call wrf_message('PartMC_init: Setting the initial meteorological &
         conditions for env_states')

    call init_wrf_to_partmc(grid, env_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
         pmc_ks, pmc_ke)

    call wrf_message('PartMC_init: Setting aerosol data for WRF domain')

    aero_data = input_aero_data

    n_parts = grid%num_particles

    call wrf_message('PartMC_init: Setting gas data for the WRF domain')

    gas_data = input_gas_data

    call pmc_mpi_barrier()

    call wrf_message('PartMC_init: Setting initial aerosol and gas states for &
         WRF domain')

    !write(mode_name,'(a)') 'sea_salt'
    !dummy = aero_data_weight_class_by_name(aero_data, mode_name)
    !dummy = aero_data_source_by_name(aero_data, mode_name)

    do j = pmc_js, pmc_je
    do k = pmc_ks, pmc_ke
    do i = pmc_is, pmc_ie
       call aero_state_zero(aero_states(i,k,j))
       call aero_state_set_weight(aero_states(i,k,j), aero_data, &
             AERO_STATE_WEIGHT_FLAT_SPECIFIED)
       call aero_state_set_n_part_ideal(aero_states(i,k,j), n_parts)
       call gas_state_set_size(gas_states(i,k,j), &
            gas_data_n_spec(gas_data))
    end do
    end do
    end do

    call wrf_message('PartMC_init: Setting initial conditions')
    t1 = MPI_Wtime()
    allocate(aero_dist_init%mode(n_modes))
    ! Recreate the aerosol distribution to sample
    write(mode_name,'(a)') 'source_A'
    i_mode = 1
    aero_dist_init%mode(i_mode)%name = mode_name
    dummy = aero_data_source_by_name(aero_data, mode_name)
    aero_dist_init%mode(i_mode)%type = &
         AERO_MODE_TYPE_LOG_NORMAL
    aero_dist_init%mode(i_mode)%char_radius = 2d-8 
    aero_dist_init%mode(i_mode)%log10_std_dev_radius = 0.25 
    allocate(vol_frac(aero_data_n_spec(aero_data)))
    vol_frac = 0.0d0
    vol_frac(1) = 1.0d0
    aero_dist_init%mode(i_mode)%vol_frac = vol_frac
    aero_dist_init%mode(i_mode)%source = dummy
    allocate(vol_frac_std(aero_data_n_spec(aero_data)))
    vol_frac_std = 0.0d0
    aero_dist_init%mode(i_mode)%vol_frac_std = vol_frac_std
    aero_dist_init%mode(i_mode)%weight_class = &
         aero_data_weight_class_by_name(aero_data, mode_name)

    if (n_modes > 1) then
    write(mode_name,'(a)') 'source_B'
    i_mode = 2 
    aero_dist_init%mode(i_mode)%name = mode_name
    dummy = aero_data_source_by_name(aero_data, mode_name)
    aero_dist_init%mode(i_mode)%type = &
         AERO_MODE_TYPE_LOG_NORMAL
    aero_dist_init%mode(i_mode)%char_radius = 2d-8
    aero_dist_init%mode(i_mode)%log10_std_dev_radius = 0.25
    vol_frac = 0.0d0
    vol_frac(1) = 1.0d0
    aero_dist_init%mode(i_mode)%vol_frac = vol_frac
    aero_dist_init%mode(i_mode)%source = dummy
    vol_frac_std = 0.0d0
    aero_dist_init%mode(i_mode)%vol_frac_std = vol_frac_std
    aero_dist_init%mode(i_mode)%weight_class = &
         aero_data_weight_class_by_name(aero_data, mode_name)
    end if

    gauss_points(1) = - (1.0/3) * sqrt(5 + 2*sqrt(10.0/7.0))
    gauss_points(2) = - (1.0/3) * sqrt(5 - 2*sqrt(10.0/7.0))
    gauss_points(3) = 0.0
    gauss_points(4) = (1.0/3) * sqrt(5 - 2*sqrt(10.0/7.0))
    gauss_points(5) = (1.0/3) * sqrt(5 + 2*sqrt(10.0/7.0))

    gauss_weights(1) = (322.0 -  13.0*sqrt(70.0)) / 900
    gauss_weights(2) = (322.0 +  13.0*sqrt(70.0)) / 900
    gauss_weights(3) = 128.0 / 225.0
    gauss_weights(4) = (322.0 +  13.0*sqrt(70.0)) / 900
    gauss_weights(5) = (322.0 - 13.0*sqrt(70.0)) / 900

    if (config_flags%do_uniform) then
       do j = pmc_js,pmc_je
       do k = pmc_ks,pmc_ke
       do i = pmc_is,pmc_ie
          b = config_flags%dx*(float(i)) / 80000.0
          a = config_flags%dx*(float(i)-1) / 80000.0
          do i_mode = 1,2
             q_value = 0.0d0
             if (config_flags%do_zero_background) then
                do q = 1,5
                   xpos = (.5*(b-a) * gauss_points(q) + .5*(a+b))
                   z = abs(xpos - .5)
                   q_value = q_value + gauss_weights(q) * 10d9 * (1.0 /  &
                     (1.0 + exp(80.0 * (z - 0.15))))
                end do
                aero_dist_init%mode(i_mode)%num_conc = q_value  / 2.0d0
                if (aero_dist_init%mode(i_mode)%num_conc  < .1 * 1d8) then
                   aero_dist_init%mode(i_mode)%num_conc  = 0.0d0
                end if
             else if (config_flags%do_constant) then
                aero_dist_init%mode(i_mode)%num_conc = 1.0d9
             else
                do q = 1,5
                   xpos = (.5*(b-a) * gauss_points(q) + .5*(a+b))
                   z = abs(xpos - .5)
                   q_value = q_value + gauss_weights(q) * (1d9 + 9d9 * (1.0 /  &
                     (1.0 + exp(80.0 * (z - 0.15)))))
                end do
                aero_dist_init%mode(i_mode)%num_conc =  q_value / 2.0d0
             end if
          end do
          call aero_state_add_aero_dist_sample(aero_states(i,k,j), aero_data, &
               aero_dist_init, 1.0/real(grid%alt(i,k,j),kind=dp), 1.0d0, 0d0, &
               .true., .true.)
       end do
       end do
       end do
    else if (config_flags%do_rotational) then

       do i_spec = 1,n_gas_emit
          gas_index(i_spec) = gas_data_spec_by_name(gas_data, &
               trim(gas_emit_names(i_spec)))
       end do
       r0 = 6.0d0
       x_center = real(global_nx,kind=dp) / 2 !50.0d0
       y_center = 1.5*(real(global_ny,kind=dp) / 2) !75.0d0
       do j = pmc_js,pmc_je
       do k = pmc_ks,pmc_ke
       do i = pmc_is,pmc_ie
          b = float(I)
          a = float(I-1)
          d = float(J)
          c = float(J-1)
          if (grid%do_emission) then
             allocate(scenario(i,k,j)%aero_emission_time(1))
             scenario(i,k,j)%aero_emission_time(1) = 0.0d0
             allocate(scenario(i,k,j)%aero_emission_rate_scale(1))
             allocate(scenario(i,k,j)%aero_emission(1))
             allocate(scenario(i,k,j)%gas_emission_time(1))
             scenario(i,k,j)%gas_emission_time(1) = 0.0d0
             allocate(scenario(i,k,j)%gas_emission_rate_scale(1))
             allocate(scenario(i,k,j)%gas_emission(1))
             allocate(scenario(i,k,j)%gas_emission(1)%mix_rat(gas_data_n_spec(gas_data)))
             scenario(i,k,j)%aero_emission_rate_scale = 1.0d0
             scenario(i,k,j)%gas_emission_rate_scale = 1.0d0
             allocate(scenario(i,k,j)%aero_emission(1)%mode(n_modes))
          end if
          do i_mode = 1,n_modes
             q_value = 0.0d0
             do ix = 1,5
             do iy = 1,5
                x_pos = (.5*(b-a) * gauss_points(ix) + .5*(a+b))
                y_pos = (.5*(d-c) * gauss_points(iy) + .5*(d+c))
                xrad = sqrt((x_pos - x_center)**2 + (y_pos - y_center)**2)
                if (config_flags%do_zero_background) then
                   q_value = q_value + gauss_weights(ix)*gauss_weights(iy) * &
                        (10d9*exp(-(xrad/r0)**2))
                else
                   q_value = q_value + gauss_weights(ix)*gauss_weights(iy) * &
                        (1d9 + 9d9*exp(-(xrad/r0)**2))
                end if
             end do
             end do
             if (i .ne. 1 .and. i .ne. global_nx .and. j .ne. 1 .and. j .ne. global_ny) then
                aero_dist_init%mode(i_mode)%num_conc = max(q_value / 4.0, 1d-15)
             else
                aero_dist_init%mode(i_mode)%num_conc = 1d-15
             endif
          end do

          call aero_state_add_aero_dist_sample(aero_states(i,k,j), & 
               aero_data, aero_dist_init, 1.0/real(grid%alt(i,k,j),kind=dp), 1.0d0, &
                0d0, .true., .true.)
          if (i == 1 .or. i == global_nx .or. j == 1 .or. j == global_ny) then
             if (.not. allocated(scenario(i,k,j)%aero_background)) then
                allocate(scenario(i,k,j)%aero_background(1))
                allocate(scenario(i,k,j)%aero_background(1)%mode(n_modes))
                allocate(scenario(i,k,j)%aero_dilution_time(1))
                allocate(scenario(i,k,j)%aero_dilution_rate(1))
                allocate(scenario(i,k,j)%gas_background(1))
                allocate(scenario(i,k,j)%gas_dilution_time(1))
                allocate(scenario(i,k,j)%gas_dilution_rate(1))
             end if
             do i_mode = 1,n_modes
                scenario(i,k,j)%aero_background(1)%mode(i_mode) = &
                     aero_dist_init%mode(i_mode)
             end do 
          end if

          if (grid%do_emission) then
             scenario(i,k,j)%aero_emission(1)%mode  = &
                        aero_dist_init%mode
             if (aero_dist_init%mode(1)%num_conc > 2d9) then
                ! Scale the emissions to be a fraction
                scenario(i,k,j)%aero_emission_rate_scale = 10.0d0
             else
                scenario(i,k,j)%aero_emission_rate_scale = 10.0d0
             end if
             scenario(i,k,j)%gas_emission(1)%mix_rat = 0.0d0
             do i_spec = 1,n_gas_emit
                scenario(i,k,j)%gas_emission(1)%mix_rat(gas_index(i_spec)) = gas_emit_values(i_spec) 
             end do
          end if
       end do
       end do
       end do

    else if (config_flags%do_seabreeze) then
       r0 = 6.0d0
       do j = pmc_js,pmc_je
       do k = pmc_ks,pmc_ke
       do i = pmc_is,pmc_ie
          b = float(I)
          a = float(I-1)
          d = float(J)
          c = float(J-1)
          do i_mode = 1,2
             q_value = 0.0d0
             do ix = 1,5
             do iy = 1,5
                x_pos = (.5*(b-a) * gauss_points(ix) + .5*(a+b))
                xrad = sqrt((x_pos - 100.0)**2)
                if (config_flags%do_zero_background) then
                   q_value = q_value + gauss_weights(ix) * &
                        (10d9*exp(-(xrad/r0)**2))
                else
                   q_value = q_value + gauss_weights(ix) * &
                        (1d9 + 9d9*exp(-(xrad/r0)**2))
                end if
             end do
             end do
             aero_dist_init%mode(i_mode)%num_conc = q_value / 4.0
          end do
          call aero_state_add_aero_dist_sample(aero_states(i,k,j), &
               aero_data, aero_dist_init, 1.0d0, 1.0d0, 0d0, .true., .true.)
          if (i == 1 .or. i == global_nx .or. j == 1 .or. j == global_ny) then
             if (.not. allocated(scenario(i,k,j)%aero_background)) then
                allocate(scenario(i,k,j)%aero_background(1))
                allocate(scenario(i,k,j)%aero_background(1)%mode(n_modes))
                allocate(scenario(i,k,j)%aero_dilution_time(1))
                allocate(scenario(i,k,j)%aero_dilution_rate(1))
                allocate(scenario(i,k,j)%gas_background(1))
                allocate(scenario(i,k,j)%gas_dilution_time(1))
                allocate(scenario(i,k,j)%gas_dilution_rate(1))
             end if
          end if
       end do
       end do
       end do
    else if (config_flags%do_ideal_init_cond .and. &
        (config_flags%periodic_x .and. config_flags%periodic_y)) then
       x0 = global_nx / 2
       y0 = global_ny / 2
       z0 = 0
       rx = 6.0
       ry = 6.0
       rz = 4.0
       pmc_o3 = gas_data_spec_by_name(gas_data, "O3")
       do j = pmc_js,pmc_je
       do k = pmc_ks,pmc_ke
       do i = pmc_is,pmc_ie
          b = float(I)
          a = float(I-1)
          d = float(J)
          c = float(J-1)
          f = float(K)
          e = float(K-1)
          q_value = 0.0d0
          do ix = 1,5
          do iy = 1,5
          do iz = 1,5
             x_pos = (.5*(b-a) * gauss_points(ix) + .5*(a+b))
             y_pos = (.5*(d-c) * gauss_points(iy) + .5*(d+c))
             z_pos = (.5*(f-e) * gauss_points(iz) + .5*(f+e))
             xrad = sqrt(((x_pos - x0)/rx)**2 + ((y_pos - y0)/ry)**2 + &
                  ((z_pos - z0)/rz)**2)
             ! FIXME: Decide if we should go with exp(-r**2) here
             q_value = q_value + gauss_weights(ix) * gauss_weights(iy) * &
                  gauss_weights(iz) * (10d9*exp(-xrad))
          end do
          end do
          end do

          do i_mode = 1,n_modes
             aero_dist_init%mode(i_mode)%num_conc = q_value / 8.0d0
          end do
          gas_states(i,k,j)%mix_rat(pmc_o3) = q_value / 8.0d0  / 1.0d8
          call aero_state_add_aero_dist_sample(aero_states(i,k,j), &
               aero_data, aero_dist_init, 1.0 &
               / real(grid%alt(i,k,j), kind=dp), 1.0d0, 0.0d0, .true., &
               .true.)

       end do
       end do
       end do
    else if (config_flags%do_ideal_init_cond) then
       x0 = 75
       y0 = 75
       z0 = 0
       rx = 6.0
       ry = 6.0
       rz = 4.0
       pmc_o3 = gas_data_spec_by_name(gas_data, "O3")
       do j = pmc_js,pmc_je
       do k = pmc_ks,pmc_ke
       do i = pmc_is,pmc_ie
          b = float(I)
          a = float(I-1)
          d = float(J)
          c = float(J-1)
          f = float(K)
          e = float(K-1)
          q_value = 0.0d0
          do ix = 1,5
          do iy = 1,5
          do iz = 1,5
             x_pos = (.5*(b-a) * gauss_points(ix) + .5*(a+b))
             y_pos = (.5*(d-c) * gauss_points(iy) + .5*(d+c))
             z_pos = (.5*(f-e) * gauss_points(iz) + .5*(f+e))
             xrad = sqrt(((x_pos - x0)/rx)**2 + ((y_pos - y0)/ry)**2 + &
                  ((z_pos - z0)/rz)**2)
             q_value = q_value + gauss_weights(ix) * gauss_weights(iy) * &
                  gauss_weights(iz) * (10d9*exp(-(xrad)))
          end do
          end do
          end do
          do i_mode = 1,n_modes
             aero_dist_init%mode(i_mode)%num_conc = max(q_value / 8.0d0, 1d-15)
          end do

          if (.not. config_flags%do_restart) then
             ! ppb will be converted to ppm
             gas_states(i,k,j)%mix_rat(pmc_o3) = 1000 * max(q_value / 8.0d0, 1d-15)
             call aero_state_add_aero_dist_sample(aero_states(i,k,j), &
                  aero_data, aero_dist_init, 1.0/real(grid%alt(i,k,j), kind=dp), 1.0d0, &
                  0d0, .true., .true.)
          end if
          if (i == 1 .or. i == global_nx .or. j == 1 .or. j == global_ny) then
             if (.not. allocated(scenario(i,k,j)%aero_background)) then
                allocate(scenario(i,k,j)%aero_background(1))
                allocate(scenario(i,k,j)%aero_background(1)%mode(n_modes))
                allocate(scenario(i,k,j)%aero_dilution_time(1))
                allocate(scenario(i,k,j)%aero_dilution_rate(1))
                allocate(scenario(i,k,j)%gas_background(1))
                allocate(scenario(i,k,j)%gas_dilution_time(1))
                allocate(scenario(i,k,j)%gas_dilution_rate(1))
             end if
             scenario(i,k,j)%aero_background(1)%mode(1) = aero_dist_init%mode(1)
             scenario(i,k,j)%aero_background(1)%mode(1)%num_conc = q_value / 8.0
             if (n_modes > 1) then
                scenario(i,k,j)%aero_background(1)%mode(2) = aero_dist_init%mode(2)
                scenario(i,k,j)%aero_background(1)%mode(2)%num_conc = q_value / 8.0
             end if
          end if


          grid%chem_bxs(1,1,1,1) = 1.0d0
       end do
       end do
       end do

       ! If we are a restart
       if (config_flags%do_restart) then
          do j = pmc_js,pmc_je
          do i = pmc_is,pmc_ie 
             call init_read_in_restart(aero_states(i,:,j), gas_states(i,:,j), i, j, &
                  pmc_ke, aero_data, gas_data, config_flags%partmc_restart_prefix, &
                  config_flags%partmc_restart_index)
          end do
          end do
       end if

    else
       print*, 'not implemented yet'
    end if

    time_ic = MPI_Wtime() - t1
    call pmc_mpi_barrier()

    call wrf_message('PartMC_init: Setting scenario for the WRF domain')

    ! UUID
    if (pmc_mpi_rank() == 0) then
      call uuid4_str(uuid)
    end if

    call pmc_mpi_bcast_string(uuid)
    grid%uuid = uuid

    ! All test cases have no emissions
    do i = pmc_is, pmc_ie
    do k = pmc_ks, pmc_ke
    do j = pmc_js, pmc_je
       call init_zero_emissions(scenario(i,k,j), aero_data, &
            gas_data)
    end do
    end do
    end do

    ! Allocate the probability arrays for all cells including boundary
    do i = pmc_is, pmc_ie
    do k = pmc_ks, pmc_ke
    do j = pmc_js, pmc_je
       n_class = aero_weight_array_n_class(aero_states(i,k,j)%awa)
       allocate(env_states(i,k,j)%prob_advection(-1:1,-1:1,-1:1,n_class))
       allocate(env_states(i,k,j)%prob_diffusion(-1:1,-1:1,-1:1,n_class))
       allocate(env_states(i,k,j)%prob_vert_diffusion(pmc_ks:pmc_ke,n_class))
    end do
    end do
    end do

    call pmc_mpi_barrier()

    print*, 'timing initialization - ICs:', time_ic, 'BCs:', time_bc

    ! FIXME: We are not being consistent in the passing of information
    ! partmc_from_wrf is needed if WRF has set the gases with its own init 
    ! routines. This is the case with uniform and rotational.
    ! For real cases, the input is still handled by PartMC files for both
    ! the aerosols and gases from the ic/bc files.
    ! The init_cond case sets gases here so partmc_to_wrf is called first.
    if (config_flags%do_ideal_init_cond) then
       call partmc_to_wrf(grid, aero_states, aero_data, gas_states, gas_data, &
            pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke)
    end if

    call partmc_from_wrf(grid, gas_states, gas_data, pmc_is, pmc_ie, pmc_js, &
         pmc_je, pmc_ks, pmc_ke)

  end subroutine init_wrf_partmc_idealized

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module wrf_pmc_init
