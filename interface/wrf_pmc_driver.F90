! Copyright (C) 2011-2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc_driver module

!> Main driver for partmc chemistry integration in WRF
module wrf_pmc_driver

  use pmc_aero_state
  use pmc_gas_state
  use pmc_gas_data
  use pmc_aero_data
  use pmc_mosaic
  use pmc_scenario
  use pmc_env_state
  use pmc_mpi
  use pmc_nucleate
  use pmc_util
  use pmc_coagulation
  use pmc_condense
  use pmc_output
  use pmc_rand

  ! WRF modules
  use module_wrf_error
  use module_domain_type
  use module_domain
  use module_model_constants
  use module_configure, only: p_qv

  USE module_state_description
  USE module_domain, ONLY : domain, domain_clock_get
  USE module_configure, ONLY : grid_config_rec_type

#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a single PartMC-MOSAIC timestep.
  subroutine partmc_timestep(grid,config_flags, scenario, env_states, &
       aero_data, aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, &
       pmc_je, pmc_ks, pmc_ke, global_nx, global_ny, global_nz)

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
    !> Model configuration settings.
    type(grid_config_rec_type), intent(in) :: config_flags
    !> Full domain of aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Full domain of aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout):: aero_states
    !> Full domain of scenario data.
    type(scenario_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenario
    !> Full domain of environmental states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Full domain of gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Full domain of gas_states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: gas_states
    !> Full WRF domain east-west dimension.
    integer :: global_nx
    !> Full WRF domain north-south dimension.
    integer :: global_ny
    !> Full WRF domain top-bottom dimension.
    integer :: global_nz

    ! Local variables
    integer :: i,j,k
    integer :: n_samp
    integer :: n_coag
    integer :: n_part_before
    integer :: n_nuc
    integer :: progress_n_nuc
    integer :: n_emit, n_dil_in,n_dil_out
    integer :: i_gas
    real(kind=dp) :: elapsed_time
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je) :: &
         old_env_states
    real(kind=dp) :: t1, t2, coag_time
    real(kind=dp) :: t_loop_end, t_loop_start
    real(kind=dp) :: emission_time, chem_time, chemistry_call
    integer :: total_particles_before, total_particles_after, len_aero_info
    integer :: allocated_aero_info
    real(kind=dp) :: wind_speed
    integer :: i_start, i_end, j_start, j_end
    character(len=PMC_UUID_LEN) :: uuid
    logical :: do_seasalt
    real :: real_time
    
    uuid =  grid%uuid

    do_seasalt = .true.

#ifdef PMC_DEBUG
    call wrf_message('PartMC_timestep: Updating the meteorology')
#endif
    ! Update the meteorology
    call wrf_to_partmc(grid, env_states, old_env_states, pmc_is, pmc_ie, &
         pmc_js, pmc_je, pmc_ks, pmc_ke)
#ifdef PMC_DEBUG
    call wrf_message('PartMC_timestep: Back from updating meteorology')
    call wrf_message("PartMC_timestep: Starting PartMC timestep")
#endif
    elapsed_time = real( &
         real_time(domain_get_time_since_sim_start(grid)),kind=dp)

    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       env_states(i,k,j)%elapsed_time = elapsed_time
    end do
    end do
    end do
    !print*, 'Current timestep = ',  grid%itimestep

    total_particles_before = 0
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       total_particles_before = total_particles_before + &
            aero_state_total_particles(aero_states(i,k,j))
    end do
    end do
    end do

    t_loop_start = MPI_Wtime()

    ! Timing variables
    coag_time = 0.0d0
    emission_time = 0.0d0
    chem_time = 0.0d0

    chemistry_call = real(grid%partmc_chem_dt,kind=dp)

    if (config_flags%periodic_x .and. config_flags%periodic_y) then
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

    do j = j_start,j_end
    do k = pmc_ks,pmc_ke
    do i = i_start,i_end
       ! FIXME: Add this correctly
       ! Nucleation
       n_part_before = aero_state_total_particles(aero_states(i,k,j))
       !call nucleate(1, env_states(i,k,j), gas_data(i,k,j), aero_data(i,k,j), &
       !     aero_states(i,k,j), gas_states(i,k,j), real(grid%dt,dp))
       n_nuc = aero_state_total_particles(aero_states(i,k,j)) - n_part_before

       ! Debugging
#ifdef PMC_DEBUG
       write(*,*) 'grid cell',i,j,k, 'contains ', &
           aero_state_total_particles(aero_states(i,k,j)), &
           aero_state_total_num_conc(aero_states(i,k,j), aero_data)
#endif

       ! Coagulation - Brownian kernel assumed
       n_coag = 0
       n_samp = 0
       if (mod(elapsed_time, chemistry_call) == 0.0d0) then
          if (grid%do_coagulation) then
             t1 = MPI_Wtime()
             call mc_coag(COAG_KERNEL_TYPE_BROWN, env_states(i,k,j), &
                  aero_data, aero_states(i,k,j), chemistry_call, &
                  n_samp, n_coag)
             t2 = MPI_Wtime()
             coag_time = coag_time + (t2 - t1)
          end if
       end if

       ! Update gas emissions
       t1 = MPI_Wtime()
       call scenario_update_gas_state(scenario(i,k,j), &
            real(grid%dt,dp), env_states(i,k,j), old_env_states(i,k,j), &
            gas_data, gas_states(i,k,j))

       n_emit = 0
       call scenario_update_aero_state(scenario(i,k,j), real(grid%dt,dp), &
            env_states(i,k,j), old_env_states(i,k,j), aero_data, &
            aero_states(i,k,j), n_emit, n_dil_in, n_dil_out, &
            grid%allow_doubling, grid%allow_halving)

       if (grid%partmc_seasalt_param > 0) then
          if (k == 1) then
             if (grid%xland(i,j) > 1.5) then
                wind_speed = (grid%u10(i,j) * grid%u10(i,j) + &
                     grid%v10(i,j) * grid%v10(i,j))**.5
                call seasalt_emissions(wind_speed, aero_states(i,k,j), &
                     aero_data, env_states(i,k,j), real(grid%dt,dp), &
                     grid%allow_doubling, grid%allow_halving, &
                     grid%partmc_seasalt_param)
             end if
          end if
       end if
       t2 = MPI_Wtime()
       emission_time = emission_time + (t2 - t1)

#ifdef PMC_DEBUG
        print*,'n_coag',n_coag,'n_samp',n_samp,'n_nuc',n_nuc,'n_emit', &
             n_emit,'n_parts', aero_state_total_particles(aero_states(i,k,j))
#endif

       if (mod(elapsed_time, chemistry_call) == 0.0d0) then
          if (grid%do_mosaic) then
             t1 = MPI_Wtime()
             call mosaic_timestep(env_states(i,k,j), aero_data, &
                  aero_states(i,k,j), gas_data, gas_states(i,k,j), &
                  grid%do_optical, uuid)
             if (grid%do_optical) then
                call compute_bulk_optical_props(grid, aero_data, &
                     aero_states(i,k,j), env_states(i,k,j))
             end if
          end if
          t2 = MPI_Wtime()
          chem_time = chem_time + (t2 - t1)
       end if

       ! Rebalance the particle population if needed 
       call aero_state_rebalance(aero_states(i,k,j), aero_data, &
            grid%allow_doubling, grid%allow_halving, &
            initial_state_warning=.false.)
       call aero_info_array_zero(aero_states(i,k,j)%aero_info_array)
    end do
    end do
    end do

    t_loop_end = MPI_Wtime()

    call memtrack(1)

#ifdef PMC_DEBUG
    write(*,*) 'Timing - Total chemistry:', t_loop_end - t_loop_start, &
         'coagulation', coag_time, 'emissions', emission_time, 'MOSAIC', &
         chem_time
#endif

#ifdef PMC_DEBUG
    call wrf_message("PartMC_timestep: Finished PartMC timestep")
#endif
#ifdef PMC_DEBUG
    total_particles_after = 0
    len_aero_info = 0
    allocated_aero_info = 0
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       total_particles_after = total_particles_after + &
            aero_state_total_particles(aero_states(i,k,j))
       len_aero_info = len_aero_info + &
           aero_info_array_n_item(aero_states(i,k,j)%aero_info_array)
    end do
    end do
    end do

    print*, 'total particles in chemistry. before:', total_particles_before, &
      ' after: ', total_particles_after, 'aero_info', len_aero_info
#endif

   call pmc_mpi_barrier()

  end subroutine partmc_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	       	Find a better home for these subroutines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test if output should be done based on wrf history interval and output
  !> if it is true 
  subroutine partmc_before_solve_io(grid, config_flags, env_states, &
       aero_data, aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, &
       pmc_je, pmc_ks, pmc_ke, global_nx, global_ny, global_nz, output_index)

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
    !>
    type(grid_config_rec_type), intent(in) :: config_flags
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: env_states
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data 
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: aero_states
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: gas_states
    !> Full WRF domain east-west dimension.
    integer, intent(in) :: global_nx
    !> Full WRF domain north-south dimension.
    integer, intent(in) :: global_ny
    !> Full WRF domain top-bottom dimension.
    integer, intent(in) :: global_nz
    !> Filename index.
    integer, intent(inout) :: output_index

    ! Local variables
    INTEGER                                    :: ialarm
    INTEGER                                    :: rc
    TYPE(WRFU_Time) :: currTime, startTime
    type(aero_state_t) :: large_aero_state
    character(len=100) :: prefix
    character(len=PMC_UUID_LEN) :: uuid
    integer :: i_repeat
    integer :: write_rank, write_n_proc
    real(kind=dp) :: time, del_t
    logical :: record_removals, record_optical, do_particle_output
    integer :: i,j,k

    do_particle_output = grid%do_gridded_output
    record_removals = grid%record_removals
    record_optical = grid%do_optical
    write_rank = 1
    write_n_proc = 1
    prefix = grid%partmc_prefix_out
    time = 60.0d0*(grid%xtime)
    del_t = grid%dt
    uuid =  grid%uuid

    ! Use the WRF_clock function and test for if it is time for output
    call WRFU_ClockGet(grid%domain_clock, CurrTime=currTime, &
         StartTime=startTime )

    if ( WRFU_AlarmIsRinging( grid%alarms( AUXHIST2_ALARM ), rc=rc ) ) then
       call partmc_process(grid, config_flags, env_states, aero_data, &
            aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, &
            pmc_je, pmc_ks, pmc_ke, global_nx, global_ny, global_nz)
    end if
    if (WRFU_AlarmIsRinging( grid%alarms( HISTORY_ALARM ),rc=rc)) then
       call partmc_output(do_particle_output, prefix, aero_data, aero_states, &
            gas_data, gas_states, env_states, &
            pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
            global_nx, global_ny, global_nz, &
            output_index, time, del_t, &
            i_repeat, record_removals, record_optical, write_rank, &
            write_n_proc, uuid)
       ! Increment the counter everytime we actually do output
       output_index = output_index + 1
    end if

  end subroutine partmc_before_solve_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Will do nothing for now unless we need to do something additional after
  ! the chemistry. In WRF, med_after_solve_io only does a time series calc
  subroutine partmc_after_solve_io()

  end subroutine partmc_after_solve_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do final time output regardless if it is the correct time interval.
  subroutine partmc_last_solve_io(grid, config_flags, env_states, aero_data, &
       aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
       pmc_ks, pmc_ke, global_nx, global_ny, global_nz, output_index)

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
    !> Model configuration settings.
    type(grid_config_rec_type), intent(in) :: config_flags
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: env_states
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data 
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: aero_states
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: gas_states
    !> Full WRF domain east-west dimension.
    integer, intent(in) :: global_nx
    !> Full WRF domain north-south dimension.
    integer, intent(in) :: global_ny
    !> Full WRF domain top-bottom dimension.
    integer, intent(in) :: global_nz
    !> Filename index.
    integer, intent(inout) :: output_index

    ! Local variables
    type(aero_state_t) :: large_aero_state
    character(len=100) :: prefix, u
    character(len=PMC_UUID_LEN) :: uuid
    integer :: i_repeat
    integer :: write_rank, write_n_proc
    real(kind=dp) :: time, del_t
    logical :: record_removals, record_optical, do_particle_output
    integer :: i,j,k

    call partmc_process(grid, config_flags, env_states, aero_data, &
         aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, &
         pmc_je, pmc_ks, pmc_ke, global_nx, global_ny, global_nz)

    do_particle_output = grid%do_gridded_output
    record_removals = grid%record_removals
    record_optical = grid%do_optical
    write_rank = 1
    write_n_proc = 1
    prefix = grid%partmc_prefix_out
    time = 60.0d0*(grid%xtime)
    del_t = grid%dt
    uuid = grid%uuid

    call partmc_output(do_particle_output, prefix, aero_data, aero_states, &
         gas_data, gas_states, env_states, &
         pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
         global_nx, global_ny, global_nz, &
         output_index, time, del_t, &
         i_repeat, record_removals, record_optical, write_rank, &
         write_n_proc, uuid)

  end subroutine partmc_last_solve_io

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs entire vertical column.
  subroutine partmc_output(do_particle_output, prefix, aero_data, aero_states, &
       gas_data, gas_states, env_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
       pmc_ks, pmc_ke, global_nx, global_ny, global_nz, &
       output_index, time, del_t, i_repeat, record_removals, &
       record_optical, write_rank, write_n_proc, uuid)

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
    !> Whether to output full particle state.
    logical, intent(in) :: do_particle_output
    !> Filename prefix.
    character(len=100), intent(in) :: prefix
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data 
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: aero_states
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: gas_states
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: env_states
    !> Full WRF domain east-west dimension.
    integer, intent(in) :: global_nx
    !> Full WRF domain north-south dimension.
    integer, intent(in) :: global_ny
    !> Full WRF domain top-bottom dimension.
    integer, intent(in) :: global_nz
    !> Filename index.
    integer, intent(inout) :: output_index
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current repeat number.
    integer, intent(in) :: i_repeat
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol optical properties.
    logical, intent(in) :: record_optical
    !> Rank to write into file.
    integer, intent(in) :: write_rank
    !> Number of processes to write into file.
    integer, intent(in) ::write_n_proc
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid

    ! Local variables
    integer :: i,j,k
    real(kind=dp) :: t1, t2
    character(len=len(prefix)+100) :: filename
    integer :: ncid
    logical :: first

    t1 = MPI_Wtime()

    if (do_particle_output) then
       do j = max(pmc_js,1), min(pmc_je, global_ny)
       do i = max(pmc_is,1), min(pmc_ie, global_nx)
          ! Old method
          !call output_column_to_file(prefix, aero_data, &
          !     aero_states(i,pmc_ks:pmc_ke,j), gas_data, &
          !     gas_states(i,pmc_ks:pmc_ke,j), &
          !     env_states(i,pmc_ks:pmc_ke,j), global_nz, &
          !     output_index,time,del_t,i_repeat, &
          !     record_removals, record_optical, uuid)
          ! New method
          call output_column_to_file_new(prefix, aero_data, &
               aero_states(i,pmc_ks:pmc_ke,j), gas_data, &
               gas_states(i,pmc_ks:pmc_ke,j), &
               env_states(i,pmc_ks:pmc_ke,j), global_nz, &
               output_index,time,del_t,i_repeat, &
               record_removals, record_optical, uuid)
       end do
       end do
    end if

    t2 = MPI_Wtime()
    print*, "particle output timing:", t2-t1

  end subroutine partmc_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Maps and converts WRF variables to PartMC-MOSAIC env_state variables.
  subroutine wrf_to_partmc(grid, env_states, old_env_states, pmc_is, pmc_ie, &
       pmc_js, pmc_je, pmc_ks, pmc_ke)

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
    !> Environmental states to be updated.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Previous environmental state.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: old_env_states

    ! Local variables
    integer :: i,j,k
    real(kind=dp) :: pressure
    real(kind=dp) :: potential_temp,temperature
    real(kind=dp) :: rh

#ifdef PMC_DEBUG
    call wrf_message('wrf_to_partmc: starting update of pressure, &
         & temperature and rh')
#endif

    do j = pmc_js, pmc_je
    do k = pmc_ks, pmc_ke
    do i = pmc_is, pmc_ie
       ! Preserve the values of the previous gas state
       old_env_states(i,k,j) = env_states(i,k,j)
       ! Calculate the pressure
       pressure = grid%p(i,k,j)+grid%pb(i,k,j)
       ! Calculate the temperature
       potential_temp = grid%t_2(i,k,j)+t0
       temperature  = potential_temp*(pressure/p0)**(rcp)
       ! Set the value in the env_state
       env_states(i,k,j)%pressure = pressure
       env_states(i,k,j)%temp = temperature
       env_states(i,k,j)%rel_humid = rel_hum(real( &
            grid%moist(i,k,j,p_qv),dp),pressure,temperature)
       ! Height and altitude
       env_states(i,k,j)%height = grid%z_at_w(i,k+1,j) &
            - grid%z_at_w(i,k,j)
       env_states(i,k,j)%altitude = grid%z(i,k,j)
       env_states(i,k,j)%z_min = grid%z_at_w(i,k,j)
       env_states(i,k,j)%z_max = grid%z_at_w(i,k+1,j)
       ! Inverse density
       env_states(i,k,j)%rrho = grid%alt(i,k,j)
       env_states(i,k,j)%cell_volume = get_grid_cell_volume(grid, i, k, j)
    end do
    end do
    end do

#ifdef PMC_DEBUG
    call wrf_message('wrf_to_partmc: finished update of pressure,  &
         & temperature and rh')
#endif

  end subroutine wrf_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rh given qv, pressure and temperature.
  real(kind=dp) function rel_hum(qv, pressure, temp)

    !> qv
    real(kind=dp) :: qv
    !> Pressure (Pa).
    real(kind=dp) :: pressure
    !> Temperature (K).
    real(kind=dp) :: temp

    real(kind=dp) :: rh

    ! This calculation was taken from module_chem_utilities (chem_prep)
    ! to be most consistent with WRF-CHEM
    rh = qv/(3.8d0*exp(17.27*(temp-273.0d0)/(temp-36.0d0))/(.01*pressure))

    rel_hum = max(.1d0,min(rh,.95d0))

  end function rel_hum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Process grid cell.
  subroutine partmc_process(grid, config_flags, env_states, &
       aero_data, aero_states, gas_data, gas_states, pmc_is, pmc_ie, pmc_js, &
       pmc_je, pmc_ks, pmc_ke, global_nx, global_ny, global_nz)

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
    !> Model configuration settings.
    type(grid_config_rec_type), intent(in) :: config_flags
    !> Full domain of environmental states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: env_states
    !> Full domain of aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Full domain of aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in):: aero_states
    !> Full domain of gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Full domain of gas_states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(in) :: gas_states
    !> Full WRF domain east-west dimension.
    integer :: global_nx
    !> Full WRF domain north-south dimension.
    integer :: global_ny
    !> Full WRF domain top-bottom dimension.
    integer :: global_nz

    !
    real(kind=dp), dimension(:), allocatable :: dry_diameters, num_concs, &
         masses, bc_masses, num_dist, mass_dist, num_concs_average, &
         num_concs_deaverage, diameters, dry_diameters_average, &
         dry_diameters_deaverage
    real(kind=dp), dimension(:), allocatable :: so4_masses, nh4_masses, &
       oc_masses, no3_masses, soa_masses, ss_masses, oin_masses, h2o_masses
    real(kind=dp), dimension(:), allocatable ::  crit_rel_humids, dry_mass_conc
    
    integer :: i, j, k, species, i_bin
    real(kind=dp) :: tot_num_conc, tot_mass_conc
    type(bin_grid_t) :: diam_grid
    real(kind=dp) :: d_gamma, d_alpha, chi
    real(kind=dp), dimension(:), allocatable :: d_gamma_array, d_alpha_array, &
         chi_array
    real(kind=dp), dimension(:), allocatable :: num_conc_by_source
    integer, dimension(:), allocatable :: n_components
    ! Environmental saturation thresholds for CCN
    integer, parameter :: n_sats = 4
    real(kind=dp), dimension(n_sats), parameter :: env_sat = [1.001,1.003,1.006,1.01]
    logical, dimension(:), allocatable :: is_ccn_active, contains_bc, &
       coagulated, is_size_range
    real(kind=dp), dimension(:), allocatable :: num_concs_active, num_concs_bc
    integer :: i_env_sat
    !
    type(aero_state_t) :: aero_state_average, aero_state_deaverage, &
         aero_state_submicron
    type(bin_grid_t) :: grid_model
    real(kind=dp), dimension(11), parameter :: mosaic_bin_edges = &
      [.1,39.0,78.0,156.0,312.0,625.0,1250.0,2500.0,5000.0,10000.0, &
         1000000.0]
    logical, parameter :: mosaic_grid = .false.

    real(kind=dp), dimension(:), allocatable :: b_scat, b_abs

    real(kind=dp) :: del_t, time
    integer :: n_bin, i_part, i_mode
    integer :: i_repeat, output_index, write_n_proc, write_rank
    character(len=PMC_UUID_LEN) :: uuid
    logical :: record_removals, record_optical
    character(len=100) :: prefix
    integer :: counter

    ! Grid cells to output
    integer, dimension(2), parameter :: loc_i = [92, 101]
    integer, dimension(2), parameter :: loc_j = [96, 100]
    integer :: p

    ! optical
    integer, parameter :: i_550nm = 3

    character(len=AERO_NAME_LEN), allocatable :: external_groups(:,:)

!    allocate(external_groups(5, 8))
    allocate(external_groups(4,8))
    external_groups(1,:) = ["OC    ", "BC    ", "      ", "      ", &
         "      ", "      ", "      ", "      "]
    external_groups(2,:) = ["API1  ", "API2  ", "LIM1  ", "LIM2  ", &
         "ARO1  ", "ARO2  ", "ALK1  ", "OLE1  "]
!    external_groups(3,:) = ["SO4   ", "NO3   ", "NH4   ", "      ", &
!         "      ", "      ", "      ", "      "]
!    external_groups(4,:) = ["Na    ", "Cl    ", "      ", "      ", &
!         "      ", "      ", "      ", "      "]
!    external_groups(5,:) = ["Ca    ", "CO3   ", "OIN   ", "      ", &
!         "      ", "      ", "      ", "      "]
     external_groups(3,:) = ["SO4   ", "NO3   ", "NH4   ", "Na    ", &
         "Cl    ", "      ", "      ", "      "]
     external_groups(4,:) = ["Ca    ", "CO3   ", "OIN   ", "      ", &
         "      ", "      ", "      ", "      "]


    ! Grid for the size distribution
    call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 100, 1d-9, 1d-3)
    grid%bin_centers = diam_grid%centers
    grid%bin_edges = diam_grid%edges

    ! make the mosaic grid
    if (mosaic_grid) then
       n_bin = size(mosaic_bin_edges)-1
       allocate(grid_model%centers(n_bin))
       allocate(grid_model%edges(n_bin+1))
       allocate(grid_model%widths(n_bin))
       grid_model%type = BIN_GRID_TYPE_LOG
       grid_model%edges = .5d0 * mosaic_bin_edges * 1d-9
       do i = 1,n_bin
          grid_model%widths(i) = (log(grid_model%edges(i+1)) - &
               log(grid_model%edges(i)))
          grid_model%centers(i) = sqrt(grid_model%edges(i+1) &
             * grid_model%edges(i))
       end do
    else
       ! Note that WRF-PartMC is going to do modes in order of size
       !  WRF PartMC a1 = MAM3 a2
       !  WRF PartMC a2 = MAM3 a1
       !  WRF PartMC a3 = MAM3 a3
       n_bin = 3
       allocate(grid_model%centers(n_bin))
       allocate(grid_model%edges(n_bin+1))
       allocate(grid_model%widths(n_bin))
       grid_model%type = BIN_GRID_TYPE_LOG
       grid_model%edges(1) = 1d-9 / 2
       grid_model%edges(2) = 1d-7 / 2
       grid_model%edges(3) = 1d-6 / 2
       grid_model%edges(4) = 1d-3 / 2
       do i = 1,n_bin
          grid_model%widths(i) = (log(grid_model%edges(i+1)) - &
               log(grid_model%edges(i)))
          grid_model%centers(i) = sqrt(grid_model%edges(i+1) &
             * grid_model%edges(i))
       end do
    end if

    ! Copy gases to the 4D chem array 
    call partmc_to_wrf(grid, aero_states, aero_data, gas_states, &
         gas_data, pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke)

    do j = pmc_js, pmc_je
    do k = pmc_ks, pmc_ke
    do i = pmc_is, pmc_ie
    grid%cell_vol(i,k,j) = get_grid_cell_volume(grid, i,k,j)

    grid%temperature(i,k,j) = env_states(i,k,j)%temp
    grid%rel_humid(i,k,j) = env_states(i,k,j)%rel_humid

    grid%density_dry_air(i,k,j) =  env_state_air_den(env_states(i,k,j))

    if (aero_state_total_particles(aero_states(i,k,j)) > 0) then
       dry_diameters = aero_state_dry_diameters(aero_states(i,k,j), &
            aero_data)
       num_concs = aero_state_num_concs(aero_states(i,k,j), aero_data)

       tot_num_conc = sum(num_concs)
       grid%tot_num_conc(i,k,j) = tot_num_conc

       grid%tot_wet_num_conc(i,k,j) = aero_state_total_num_conc_wet( &
            aero_states(i,k,j), aero_data)

       masses = aero_state_masses(aero_states(i,k,j), aero_data)
       tot_mass_conc = sum(masses * num_concs)
       grid%tot_mass_conc(i,k,j) = tot_mass_conc
       num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)

       grid%num_dist(i,k,j,PARAM_FIRST_SCALAR:PARAM_NUM_num_dist) = num_dist

       mass_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, &
            num_concs * masses)
       grid%mass_dist(i,k,j,PARAM_FIRST_SCALAR:PARAM_NUM_mass_dist) = mass_dist

       num_conc_by_source = aero_state_num_concs_by_source(aero_states(i,k,j), &
            aero_data)
       grid%num_conc_source(i,k,j,PARAM_FIRST_SCALAR:PARAM_FIRST_SCALAR + & 
            aero_data_n_source(aero_data) - 1) = num_conc_by_source
       ! FIXME: We should test for each species with aero_data_spec_by_name
       if (aero_data_n_spec(aero_data) > 1) then
       bc_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"BC"/))
       grid%aero_mass(i,k,j,p_pmc_bc) = sum(bc_masses * num_concs)

       oc_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"OC"/))
       grid%aero_mass(i,k,j,p_pmc_oc) = sum(oc_masses * num_concs)

       oin_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"OIN"/))
       grid%aero_mass(i,k,j,p_pmc_oin) = sum(oin_masses * num_concs)

       so4_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"SO4"/))
       grid%aero_mass(i,k,j,p_pmc_so4) = sum(so4_masses * num_concs)

       nh4_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"NH4"/))
       grid%aero_mass(i,k,j,p_pmc_nh4) = sum(nh4_masses * num_concs)

       no3_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"NO3"/))
       grid%aero_mass(i,k,j,p_pmc_no3) = sum(no3_masses * num_concs)

       ss_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"Na"/))
       grid%aero_mass(i,k,j,p_pmc_na) = sum(ss_masses * num_concs)

       ss_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"Cl"/))
       grid%aero_mass(i,k,j,p_pmc_cl) = sum(ss_masses * num_concs)

       h2o_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"H2O"/))
       grid%aero_mass(i,k,j,p_pmc_h2o) = sum(h2o_masses * num_concs) 

       ! SOA species
       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"ARO1"/))
       grid%aero_mass(i,k,j,p_pmc_aro1) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"ARO2"/))
       grid%aero_mass(i,k,j,p_pmc_aro2) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"ALK1"/))
       grid%aero_mass(i,k,j,p_pmc_alk1) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"OLE1"/))
       grid%aero_mass(i,k,j,p_pmc_ole1) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"API1"/))
       grid%aero_mass(i,k,j,p_pmc_api1) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"API2"/))
       grid%aero_mass(i,k,j,p_pmc_api2) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"LIM1"/))
       grid%aero_mass(i,k,j,p_pmc_lim1) = sum(soa_masses * num_concs)

       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"LIM2"/))
       grid%aero_mass(i,k,j,p_pmc_lim2) = sum(soa_masses * num_concs)

       ! 
       soa_masses = aero_state_masses(aero_states(i,k,j), aero_data, &
            include=["ARO1", "ARO2","ALK1", "OLE1", "API1", "API2", "LIM1", "LIM2"])
       end if

       if (config_flags%do_advanced_process) then
       ! Mixing state information
       call aero_state_mixing_state_metrics(aero_states(i,k,j), &
            aero_data, d_alpha, d_gamma, chi, exclude=["H2O"])

       grid%chi(i,k,j) = chi
       grid%d_alpha(i,k,j) = d_alpha
       grid%d_gamma(i,k,j) = d_gamma

       call aero_state_mixing_state_metrics(aero_states(i,k,j), &
            aero_data, d_alpha, d_gamma, chi, &
            exclude=["H2O"], group=["BC ", "OC "])
  
       grid%chi_ccn(i,k,j) = chi
       grid%d_alpha_ccn(i,k,j) = d_alpha
       grid%d_gamma_ccn(i,k,j) = d_gamma

       call aero_state_mixing_state_metrics(aero_states(i,k,j), &
            aero_data, d_alpha, d_gamma, chi, &
            exclude=["H2O"], group=["BC "])

       grid%chi_opt(i,k,j) = chi
       grid%d_alpha_opt(i,k,j) = d_alpha
       grid%d_gamma_opt(i,k,j) = d_gamma

       call aero_state_mixing_state_metrics_by_size(aero_states(i,k,j), &
            aero_data, grid_model, d_alpha_array, d_gamma_array, &
            chi_array, exclude=["H2O"])

       i_bin = 1
       grid%chi_species_a1(i,k,j) = chi_array(i_bin)
       grid%d_alpha_species_a1(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_species_a1(i,k,j) = d_gamma_array(i_bin)
       i_bin = 2
       grid%chi_species_a2(i,k,j) = chi_array(i_bin)
       grid%d_alpha_species_a2(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_species_a2(i,k,j) = d_gamma_array(i_bin)
       i_bin = 3
       grid%chi_species_a3(i,k,j) = chi_array(i_bin)
       grid%d_alpha_species_a3(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_species_a3(i,k,j) = d_gamma_array(i_bin)

       call aero_state_mixing_state_metrics_by_size(aero_states(i,k,j), &
            aero_data, grid_model, d_alpha_array, d_gamma_array, &
            chi_array, exclude=["H2O"], group=["BC ", "OC "])

       i_bin = 1
       grid%chi_ccn_a1(i,k,j) = chi_array(i_bin)
       grid%d_alpha_ccn_a1(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_ccn_a1(i,k,j) = d_gamma_array(i_bin)
       i_bin = 2
       grid%chi_ccn_a2(i,k,j) = chi_array(i_bin)
       grid%d_alpha_ccn_a2(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_ccn_a2(i,k,j) = d_gamma_array(i_bin)
       i_bin = 3
       grid%chi_ccn_a3(i,k,j) = chi_array(i_bin)
       grid%d_alpha_ccn_a3(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_ccn_a3(i,k,j) = d_gamma_array(i_bin)

       call aero_state_mixing_state_metrics_by_size(aero_states(i,k,j), &
            aero_data, grid_model, d_alpha_array, d_gamma_array, &
            chi_array, exclude=["H2O"], group=["BC "])

       i_bin = 1
       grid%chi_opt_a1(i,k,j) = chi_array(i_bin)
       grid%d_alpha_opt_a1(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_opt_a1(i,k,j) = d_gamma_array(i_bin)
       i_bin = 2
       grid%chi_opt_a2(i,k,j) = chi_array(i_bin)
       grid%d_alpha_opt_a2(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_opt_a2(i,k,j) = d_gamma_array(i_bin)
       i_bin = 3
       grid%chi_opt_a3(i,k,j) = chi_array(i_bin)
       grid%d_alpha_opt_a3(i,k,j) = d_alpha_array(i_bin)
       grid%d_gamma_opt_a3(i,k,j) = d_gamma_array(i_bin)

       ! Submicron mixing state
       aero_state_submicron = aero_states(i,k,j)
       do i_part = aero_state_n_part(aero_state_submicron),1,-1
          if (dry_diameters(i_part) > 1d-6) then
             call aero_state_remove_particle_no_info(aero_state_submicron, &
                  i_part)
          end if
       end do
       call aero_state_mixing_state_metrics(aero_state_submicron, &
            aero_data, d_alpha, d_gamma, chi, exclude=["H2O"])
       grid%chi_submicron(i,k,j) = chi
       grid%d_alpha_submicron(i,k,j) = d_alpha
       grid%d_gamma_submicron(i,k,j) = d_gamma

       call aero_state_mixing_state_metrics(aero_state_submicron, &
            aero_data, d_alpha, d_gamma, chi, &
            exclude=["H2O"], group=["BC ", "OC "])
       grid%chi_ccn_submicron(i,k,j) = chi
       grid%d_alpha_ccn_submicron(i,k,j) = d_alpha
       grid%d_gamma_ccn_submicron(i,k,j) = d_gamma

       call aero_state_mixing_state_metrics(aero_state_submicron, &
            aero_data, d_alpha, d_gamma, chi, &
            exclude=["H2O"], group=["BC "])
       grid%chi_opt_submicron(i,k,j) = chi
       grid%d_alpha_opt_submicron(i,k,j) = d_alpha
       grid%d_gamma_opt_submicron(i,k,j) = d_gamma

       ! CCN information
       ! Ordering: inner loop is modes, outer loop is supersaturation
       crit_rel_humids = aero_state_crit_rel_humids(aero_states(i,k,j), &
            aero_data, env_states(i,k,j))
       counter = 1 
       do i_env_sat = 1,n_sats
          is_ccn_active = crit_rel_humids < env_sat(i_env_sat)
          num_concs_active = pack(num_concs, is_ccn_active)
          grid%pmc_ccn_conc(i,k,j,PARAM_FIRST_SCALAR+i_env_sat-1) = &
             sum(num_concs_active)
          if (.not. mosaic_grid) then
             do i_mode = 1,bin_grid_size(grid_model)
                is_size_range = 2.0d0 * grid_model%edges(i_mode) &
                     < dry_diameters .and. dry_diameters &
                     <= grid_model%edges(i_mode+1) * 2.0d0
                num_concs_active = pack(num_concs, is_ccn_active &
                     .and. is_size_range)
                grid%pmc_ccn_conc_modes_pr(i,k,j,PARAM_FIRST_SCALAR+counter-1) = &
                     sum(num_concs_active)
                counter = counter + 1
             end do
          end if
       end do

       ! Get average
       aero_state_average = aero_states(i,k,j)
       call aero_state_make_dry(aero_state_average, aero_data) 
       call aero_state_bin_average_comp(aero_state_average, grid_model, &
            aero_data)
       dry_diameters_average = aero_state_dry_diameters(aero_state_average, &
            aero_data)
       crit_rel_humids = aero_state_crit_rel_humids(aero_state_average, &
            aero_data, env_states(i,k,j))
       num_concs_average = aero_state_num_concs(aero_state_average, aero_data)
       counter = 1
       do i_env_sat = 1,n_sats
          is_ccn_active = crit_rel_humids < env_sat(i_env_sat)
          num_concs_active = pack(num_concs_average, is_ccn_active)
          grid%pmc_ccn_conc_internal(i,k,j,PARAM_FIRST_SCALAR+i_env_sat-1) = &
             sum(num_concs_active)
          if (.not. mosaic_grid) then
             do i_mode = 1,bin_grid_size(grid_model)
                is_size_range = 2.0d0 * grid_model%edges(i_mode) &
                     < dry_diameters_average .and. dry_diameters_average &
                     <= grid_model%edges(i_mode+1) * 2.0d0
                num_concs_active = pack(num_concs_average, is_ccn_active &
                     .and. is_size_range)
                grid%pmc_ccn_conc_modes_internal(i,k,j,PARAM_FIRST_SCALAR+counter-1) = &
                     sum(num_concs_active)
                counter = counter + 1
             end do
          end if
       end do

       ! Get deaverage
!       aero_state_deaverage = aero_states(i,k,j)
!       call aero_state_make_dry(aero_state_deaverage, aero_data)
!       call aero_state_bin_deaverage_comp(aero_state_deaverage, grid_model, aero_data, &
!            external_groups)
!       dry_diameters_deaverage = aero_state_dry_diameters(aero_state_deaverage, &
!            aero_data)
!       crit_rel_humids = aero_state_crit_rel_humids(aero_state_deaverage, &
!            aero_data, env_states(i,k,j))
!       num_concs_deaverage = aero_state_num_concs(aero_state_deaverage, aero_data)
!       counter = 1
!       do i_env_sat = 1,n_sats
!          is_ccn_active = crit_rel_humids < env_sat(i_env_sat)
!          num_concs_active = pack(num_concs_deaverage, is_ccn_active)
!          grid%pmc_ccn_conc_external(i,k,j,PARAM_FIRST_SCALAR+i_env_sat-1) = &
!             sum(num_concs_active)
!          ! loop over modes
!          if (.not. mosaic_grid) then
!             do i_mode = 1,bin_grid_size(grid_model)
!                is_size_range = 2.0d0 * grid_model%edges(i_mode) &
!                     < dry_diameters_deaverage .and. dry_diameters_deaverage &
!                     <= 2.0d0 * grid_model%edges(i_mode+1)
!                num_concs_active = pack(num_concs_deaverage, is_ccn_active &
!                     .and. is_size_range)
!                grid%pmc_ccn_conc_modes_external(i,k,j,PARAM_FIRST_SCALAR+counter-1) = &
!                     sum(num_concs_active)
!                counter = counter + 1
!             end do
!          end if
!       end do

       ! Total BC number concentration for determining fraction
       contains_bc = bc_masses > 0.0d0
       num_concs_bc = pack(num_concs, contains_bc)
       grid%tot_bc_num_conc(i,k,j) = sum(num_concs_bc)
       is_ccn_active = crit_rel_humids < env_sat(2)
       grid%tot_bc_num_conc_aged(i,k,j) = sum(pack(num_concs, &
            (contains_bc .and. is_ccn_active)))

       ! total number of computational particles
       grid%n_parts(i,k,j) = aero_state_n_part(aero_states(i,k,j))

       n_components = aero_state_num_components(aero_states(i,k,j))
       coagulated = n_components > 1
       grid%tot_coagulation_num_conc(i,k,j) = sum(pack(num_concs, coagulated))

       ! Hydrophobic mass concentration
       grid%tot_hydrophobic_mass_conc(i,k,j) = sum( &
            aero_state_masses(aero_states(i,k,j), aero_data, &
            include=(/"BC ","OC ", "OIN"/)) * num_concs)
       grid%tot_hydrophylic_mass_conc(i,k,j) = sum( &
            aero_state_masses(aero_states(i,k,j), aero_data, &
            exclude=(/"H2O"/)) * num_concs) - &
            grid%tot_hydrophobic_mass_conc(i,k,j)

       ! PM mass concentrations
       dry_mass_conc = num_concs * aero_state_masses(aero_states(i,k,j), &
            aero_data, exclude=(/"H2O"/))
       grid%PM1_mass_conc(i,k,j) = sum(pack(dry_mass_conc, dry_diameters < 1d-6))
       grid%PM25_mass_conc(i,k,j) = sum(pack(dry_mass_conc, dry_diameters < 2.5d-6))
       grid%PM10_mass_conc(i,k,j) = sum(pack(dry_mass_conc, dry_diameters < 10d-6))

       ! Number concentration of each mode
       if (.not. mosaic_grid) then
          i_mode = 1
          is_size_range = 2 * grid_model%edges(i_mode) < dry_diameters &
               .and. dry_diameters <= 2 * grid_model%edges(i_mode+1)
          grid%num_conc_a1 = sum(pack(num_concs, is_size_range))
          grid%mass_conc_a1 = sum(pack(dry_mass_conc, is_size_range))

          i_mode = 2
          is_size_range = 2 * grid_model%edges(i_mode) < dry_diameters &
               .and. dry_diameters <= 2 * grid_model%edges(i_mode+1)
          grid%num_conc_a2 = sum(pack(num_concs, is_size_range))
          grid%mass_conc_a2 = sum(pack(dry_mass_conc, is_size_range))

          i_mode = 3
          is_size_range = 2 * grid_model%edges(i_mode) < dry_diameters &
               .and. dry_diameters <= 2 * grid_model%edges(i_mode+1)
          grid%num_conc_a3 = sum(pack(num_concs, is_size_range))
          grid%mass_conc_a3 = sum(pack(dry_mass_conc, is_size_range))
       end if

       ! Optical properties
       if (grid%do_optical) then
          grid%scat_aer_550(i,k,j) = aero_state_scattering(aero_states(i,k,j), &
               aero_data, i_550nm)
          grid%ext_aer_550(i,k,j) = aero_state_absorption(aero_states(i,k,j), &
               aero_data, i_550nm) + grid%scat_aer_550(i,k,j)
          b_scat = aero_state_get_bin_scat(aero_states(i,k,j), aero_data, &
               grid_model, i_550nm)
          b_abs = aero_state_get_bin_abs(aero_states(i,k,j), aero_data, &
               grid_model, i_550nm)

          grid%scat_aer_550_pr_a1(i,k,j) = b_scat(1)
          grid%scat_aer_550_pr_a2(i,k,j) = b_scat(2)
          grid%scat_aer_550_pr_a3(i,k,j) = b_scat(3)

          grid%ext_aer_550_pr_a1(i,k,j) = b_scat(1) + b_abs(1)
          grid%ext_aer_550_pr_a2(i,k,j) = b_scat(2) + b_abs(2)
          grid%ext_aer_550_pr_a3(i,k,j) = b_scat(3) + b_abs(3)

          ! Internal mixed optical properties
          call condense_equilib_particles(env_states(i,k,j), aero_data, aero_state_average)

!          call mosaic_compute_single_aero_optical(env_states(i,k,j), aero_data, &
!               aero_state_average, gas_data, gas_states(i,k,j))
          grid%scat_aer_550_internal(i,k,j) = aero_state_scattering(aero_state_average, &
               aero_data, i_550nm)
          grid%ext_aer_550_internal(i,k,j) = aero_state_absorption(aero_state_average, &
               aero_data, i_550nm) + grid%scat_aer_550_internal(i,k,j)

          b_scat = aero_state_get_bin_scat(aero_state_average, aero_data, &
               grid_model, i_550nm)
          b_abs = aero_state_get_bin_abs(aero_state_average, aero_data, &
               grid_model, i_550nm)

          grid%scat_aer_550_internal_a1(i,k,j) = b_scat(1)
          grid%scat_aer_550_internal_a2(i,k,j) = b_scat(2)
          grid%scat_aer_550_internal_a3(i,k,j) = b_scat(3)

          grid%ext_aer_550_internal_a1(i,k,j) = b_scat(1) + b_abs(1)
          grid%ext_aer_550_internal_a2(i,k,j) = b_scat(2) + b_abs(2)
          grid%ext_aer_550_internal_a3(i,k,j) = b_scat(3) + b_abs(3)

          ! External mixed optical properties
          call condense_equilib_particles(env_states(i,k,j), aero_data, &
               aero_state_deaverage)

!          call mosaic_compute_single_aero_optical(env_states(i,k,j), aero_data, &
!               aero_state_deaverage, gas_data, gas_states(i,k,j))
          grid%scat_aer_550_external(i,k,j) = aero_state_scattering( &
               aero_state_deaverage, aero_data, i_550nm)
          grid%ext_aer_550_external(i,k,j) = aero_state_absorption( &
               aero_state_deaverage, aero_data, i_550nm) &
               + grid%scat_aer_550_external(i,k,j)
          b_scat = aero_state_get_bin_scat(aero_state_deaverage, &
               aero_data, grid_model, i_550nm)
          b_abs = aero_state_get_bin_abs(aero_state_deaverage, &
               aero_data, grid_model, i_550nm)

          grid%scat_aer_550_external_a1(i,k,j) = b_scat(1)
          grid%scat_aer_550_external_a2(i,k,j) = b_scat(2)
          grid%scat_aer_550_external_a3(i,k,j) = b_scat(3)

          grid%ext_aer_550_external_a1(i,k,j) = b_scat(1) + b_abs(1)
          grid%ext_aer_550_external_a2(i,k,j) = b_scat(2) + b_abs(2)
          grid%ext_aer_550_external_a3(i,k,j) = b_scat(3) + b_abs(3)
       end if
    end if ! End loop for aerosol species
    end if ! End loop for existence of paricles
    end do
    end do
    end do

#ifndef PMC_CONSTANT_VEL
    ! FIXME: Here we will select grid cells to output more frequently
    ! We will need to come up with a way to change this
    do i = pmc_is, pmc_ie
    do j = pmc_js, pmc_je
    do p = 1, size(loc_i)
    if (i == loc_i(p) .and. j == loc_j(p)) then
       record_removals = grid%record_removals
       record_optical = grid%do_optical
       write_rank = 1
       write_n_proc = 1
       prefix = "out2/grid_cell"
       time = 60.0d0*(grid%xtime)
       del_t = grid%dt
       uuid =  grid%uuid
       ! FIXME: This isn't flexible. Needs to be the output frequency in minutes
       output_index = grid%xtime / 10
       call output_column_to_file_new(prefix, aero_data, &
               aero_states(i,1,j), gas_data, gas_states(i,1,j), &
               env_states(i,1,j), 1, output_index, time, del_t, i_repeat, &
               record_removals, record_optical, uuid)
    end if
    end do
    end do
    end do
#endif

    call pmc_mpi_barrier()

  end subroutine partmc_process

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_column_to_file_new(prefix, aero_data, aero_state, gas_data, &
       gas_state, env_state, nz, index, time, del_t, i_repeat, &
       record_removals, record_optical, uuid, write_rank, write_n_proc)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), dimension(nz), intent(in) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), dimension(nz), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), dimension(nz), intent(in) :: env_state
    !>
    integer, intent(in) :: nz
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Current repeat number.
    integer, intent(in) :: i_repeat
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    !> Whether to output aerosol optical properties.
    logical, intent(in) :: record_optical
    !> UUID of the simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid
    !> Rank to write into file.
    integer, intent(in), optional :: write_rank
    !> Number of processes to write into file.
    integer, intent(in), optional :: write_n_proc

    character(len=len(prefix)+100) :: filename
    integer :: ncid, status
    integer :: k, i_part, i_part_g
    integer :: total_particles, total_components
    real(kind=dp), allocatable :: gas_mixing_ratios(:,:), &
         aero_particle_mass(:,:)
    integer :: dimid_gas_species, dimid_num_levels, dimid_aero_particle
    integer :: dimid_aero_components, dimid_aero_species
    integer :: dimid_aero_weight_group, dimid_aero_weight_class
    integer :: dimid_optical_wavelengths
    integer :: n_group, n_class, i_group, i_class
    real(kind=dp), allocatable :: weight_exponent(:,:,:), &
         weight_magnitude(:,:,:)
    integer, allocatable :: weight_type(:,:,:)
    ! All the particle things
    integer :: n_parts(nz)
    integer :: part_start_index(nz)
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
    real(kind=dp), allocatable :: aero_component_create_time(:)
    integer, allocatable :: aero_component_len(:)
    integer, allocatable :: aero_component_start_ind(:)
    integer :: array_position, i_comp
    integer :: n_swbands
   
#ifdef PMC_USE_WRF
    write(filename,'(a,a,i3.3,a,i3.3,a,i8.8,a)') trim(prefix), &
         '_', env_state(1)%ix, '_', env_state(1)%iy, '_', index, '.nc'
    call pmc_nc_open_write(filename, ncid)
    call pmc_nc_write_info(ncid, uuid, &
         "WRF-PartMC version " // trim(PARTMC_VERSION), write_rank, &
          write_n_proc)
    call write_time(ncid, time, del_t, index)

    call gas_data_output_netcdf(gas_data, ncid)
    call aero_data_output_netcdf(aero_data, ncid)

    total_particles = 0
    total_components = 0
    do k = 1,nz
       n_parts(k) = aero_state_n_part(aero_state(k))
       part_start_index(k) = total_particles + 1
       total_particles = total_particles + aero_state_n_part(aero_state(k))
       total_components = total_components + &
            aero_state_n_components(aero_state(k))
    end do

    allocate(gas_mixing_ratios(nz,gas_data_n_spec(gas_data)))  
    allocate(aero_particle_mass(total_particles, aero_data_n_spec(aero_data)))
    allocate(aero_num_conc(total_particles))
    allocate(aero_component_len(total_particles))
    allocate(aero_component_start_ind(total_particles))
    allocate(aero_component_particle_num(total_components))
    allocate(aero_component_source_num(total_components))
    allocate(aero_component_create_time(total_components))
    allocate(aero_particle_weight_group(total_particles))
    allocate(aero_particle_weight_class(total_particles))
    allocate(aero_water_hyst_leg(total_particles))
    allocate(aero_id(total_particles))
    allocate(aero_least_create_time(total_particles))
    allocate(aero_greatest_create_time(total_particles))
    n_swbands = 5 ! size(aero_state(1)%apa%particle(1)%absorb_cross_sect)
    allocate(aero_absorb_cross_sect(total_particles,n_swbands))
    allocate(aero_scatter_cross_sect(total_particles,n_swbands))
    allocate(aero_asymmetry(total_particles,n_swbands))
    allocate(aero_refract_shell_real(total_particles,n_swbands))
    allocate(aero_refract_shell_imag(total_particles,n_swbands))
    allocate(aero_refract_core_real(total_particles,n_swbands))
    allocate(aero_refract_core_imag(total_particles,n_swbands))
    allocate(aero_core_vol(total_particles))

    i_part_g = 1
    do k = 1, nz
       gas_mixing_ratios(k,:) = gas_state(k)%mix_rat
       do i_part = 1,aero_state_n_part(aero_state(k))
          aero_particle_mass(i_part_g, :) &
               = aero_state(k)%apa%particle(i_part)%vol * aero_data%density
          aero_component_len(i_part_g) = aero_particle_n_components( &
              aero_state(k)%apa%particle(i_part))
          if (i_part_g == 1) then
             aero_component_start_ind(i_part_g) = 1
          else
             aero_component_start_ind(i_part_g) = &
                  aero_component_start_ind(i_part_g - 1) &
                  + aero_component_len(i_part_g - 1)
          end if
          do i_comp = 1, aero_component_len(i_part_g)
             array_position = aero_component_start_ind(i_part_g) + i_comp - 1
             aero_component_particle_num(array_position) = i_part
             aero_component_source_num(array_position) =  &
                  aero_state(k)%apa%particle(i_part)%component(i_comp)%source_id
             aero_component_create_time(array_position) = &
                  aero_state(k)%apa%particle(i_part)%component(i_comp)%create_time
          end do
          aero_particle_weight_group(i_part_g) &
               = aero_state(k)%apa%particle(i_part)%weight_group
          aero_particle_weight_class(i_part_g) &
               = aero_state(k)%apa%particle(i_part)%weight_class
          aero_water_hyst_leg(i_part_g) &
               = aero_state(k)%apa%particle(i_part)%water_hyst_leg
          aero_num_conc(i_part_g) &
               = aero_state_particle_num_conc(aero_state(k), &
               aero_state(k)%apa%particle(i_part), aero_data)
          aero_id(i_part_g) = aero_state(k)%apa%particle(i_part)%id
          aero_least_create_time(i_part_g) &
               = aero_state(k)%apa%particle(i_part)%least_create_time
          aero_greatest_create_time(i_part_g) &
               = aero_state(k)%apa%particle(i_part)%greatest_create_time
          if (record_optical) then
             aero_absorb_cross_sect(i_part_g,:) &
                  = aero_state(k)%apa%particle(i_part)%absorb_cross_sect
             aero_scatter_cross_sect(i_part_g,:) &
                  = aero_state(k)%apa%particle(i_part)%scatter_cross_sect
             aero_asymmetry(i_part_g,:) = aero_state(k)%apa%particle(i_part)%asymmetry
             aero_refract_shell_real(i_part_g,:) &
                  = real(aero_state(k)%apa%particle(i_part)%refract_shell)
             aero_refract_shell_imag(i_part_g,:) &
                  = aimag(aero_state(k)%apa%particle(i_part)%refract_shell)
             aero_refract_core_real(i_part_g,:) &
                  = real(aero_state(k)%apa%particle(i_part)%refract_core)
             aero_refract_core_imag(i_part_g,:) &
                  = aimag(aero_state(k)%apa%particle(i_part)%refract_core)
             aero_core_vol(i_part_g) = aero_state(k)%apa%particle(i_part)%core_vol
          end if
       i_part_g = i_part_g + 1
       end do
    end do

    ! try to get the dimension ID
    call aero_data_netcdf_dim_aero_species(aero_data, ncid, &
         dimid_aero_species)
    status = nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle)
    if (status /= NF90_NOERR) then
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "aero_particle", &
         total_particles, dimid_aero_particle))
    call pmc_nc_check(nf90_enddef(ncid))
    end if

    status = nf90_inq_dimid(ncid, "number_levels", dimid_num_levels)
    if (status /= NF90_NOERR) then 
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "number_levels", nz, &
         dimid_num_levels))
    call pmc_nc_check(nf90_enddef(ncid))
    end if

    call gas_data_netcdf_dim_gas_species(gas_data, ncid, &
         dimid_gas_species)

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_components", dimid_aero_components)
    if (status /= NF90_NOERR) then 
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_dim(ncid, "aero_components", &
         total_components, dimid_aero_components))
    call pmc_nc_check(nf90_enddef(ncid))
    end if

    if (record_optical) then
       call aero_state_netcdf_dim_optical_wavelengths(aero_state(1), ncid, &
            dimid_optical_wavelengths)
    end if

    call pmc_nc_write_real_2d(ncid, gas_mixing_ratios, &
         "gas_mixing_ratio", (/ dimid_num_levels, dimid_gas_species /), unit="ppb", &
         long_name="mixing ratios of gas species")

    call pmc_nc_write_integer_1d(ncid, part_start_index, &
         "part_start_index", (/ dimid_num_levels /), unit="(1)", &
         long_name="first particle index for each grid cell")
    call pmc_nc_write_integer_1d(ncid, n_parts, "n_parts", &
         (/ dimid_num_levels /), unit="(1)", &
         long_name="number of particles for each grid cell")

    call pmc_nc_write_real_2d(ncid, aero_particle_mass, &
         "aero_particle_mass", (/ dimid_aero_particle, &
         dimid_aero_species /), unit="kg", &
         long_name="constituent masses of each aerosol particle")

    call pmc_nc_write_integer_1d(ncid, aero_component_len, &
         "aero_component_len", (/ dimid_aero_particle /), &
         long_name="number of aero_components for each aerosol particle")
    call pmc_nc_write_integer_1d(ncid, aero_component_start_ind, &
         "aero_component_start_ind", (/ dimid_aero_particle /), &
         long_name="start index of aero_component for each aerosol particle")
    call pmc_nc_write_integer_1d(ncid, aero_component_particle_num, &
         "aero_component_particle_num", (/ dimid_aero_components /), &
         long_name="associated aerosol particle number for each component")
    call pmc_nc_write_integer_1d(ncid, aero_component_source_num, &
         "aero_component_source_num", (/ dimid_aero_components /), &
         long_name="associated source number for each component")
    call pmc_nc_write_real_1d(ncid, aero_component_create_time, &
         "aero_component_create_time", (/ dimid_aero_components /), &
         long_name="associated creation time for each component")

    call pmc_nc_write_integer_1d(ncid, aero_particle_weight_group, &
         "aero_particle_weight_group", (/ dimid_aero_particle /), &
         long_name="weight group number of each aerosol particle")
    call pmc_nc_write_integer_1d(ncid, aero_particle_weight_class, &
         "aero_particle_weight_class", (/ dimid_aero_particle /), &
         long_name="weight class number of each aerosol particle")
    call pmc_nc_write_integer_1d(ncid, aero_water_hyst_leg, &
         "aero_water_hyst_leg", (/ dimid_aero_particle /), &
         long_name="leg of the water hysteresis curve leg of each "&
         // "aerosol particle")
    call pmc_nc_write_real_1d(ncid, aero_num_conc, &
         "aero_num_conc", (/ dimid_aero_particle /), unit="m^{-3}", &
         long_name="number concentration for each particle")
    call pmc_nc_write_integer64_1d(ncid, aero_id, &
         "aero_id", (/ dimid_aero_particle /), &
         long_name="unique ID number of each aerosol particle")
    call pmc_nc_write_real_1d(ncid, aero_least_create_time, &
         "aero_least_create_time", (/ dimid_aero_particle /), unit="s", &
         long_name="least creation time of each aerosol particle", &
         description="least (earliest) creation time of any original " &
         // "constituent particles that coagulated to form each " &
         // "particle, measured from the start of the simulation")
    call pmc_nc_write_real_1d(ncid, aero_greatest_create_time, &
         "aero_greatest_create_time", (/ dimid_aero_particle /), &
         unit="s", &
         long_name="greatest creation time of each aerosol particle", &
         description="greatest (latest) creation time of any original " &
         // "constituent particles that coagulated to form each " &
         // "particle, measured from the start of the simulation")
    if (record_optical) then
       call pmc_nc_write_real_2d(ncid, aero_absorb_cross_sect, &
            "aero_absorb_cross_sect", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="m^2", &
            long_name="optical absorption cross sections of each " &
            // "aerosol particle")
       call pmc_nc_write_real_2d(ncid, aero_scatter_cross_sect, &
            "aero_scatter_cross_sect", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="m^2", &
            long_name="optical scattering cross sections of each " &
            // "aerosol particle")
       call pmc_nc_write_real_2d(ncid, aero_asymmetry, &
            "aero_asymmetry", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="1", &
            long_name="optical asymmetry parameters of each " &
            // "aerosol particle")
       call pmc_nc_write_real_2d(ncid, aero_refract_shell_real, &
            "aero_refract_shell_real", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="1", &
            long_name="real part of the refractive indices of the " &
            // "shell of each aerosol particle")
       call pmc_nc_write_real_2d(ncid, aero_refract_shell_imag, &
            "aero_refract_shell_imag", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="1", &
            long_name="imaginary part of the refractive indices of " &
            // "the shell of each aerosol particle")
       call pmc_nc_write_real_2d(ncid, aero_refract_core_real, &
            "aero_refract_core_real", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="1", &
            long_name="real part of the refractive indices of the core " &
            // "of each aerosol particle")
       call pmc_nc_write_real_2d(ncid, aero_refract_core_imag, &
            "aero_refract_core_imag", (/ dimid_aero_particle, &
            dimid_optical_wavelengths /), unit="1", &
            long_name="imaginary part of the refractive indices of " &
            // "the core of each aerosol particle")
       call pmc_nc_write_real_1d(ncid, aero_core_vol, &
            "aero_core_vol", (/ dimid_aero_particle /), unit="m^3", &
            long_name="volume of the optical cores of each " &
            // "aerosol particle")
    end if

    ! Weighting information needed for restarts
    call aero_weight_netcdf_dim_aero_weight_group(aero_state(1)%awa, ncid, &
         dimid_aero_weight_group)
    call aero_weight_netcdf_dim_aero_weight_class(aero_state(1)%awa, ncid, &
         dimid_aero_weight_class)

    n_group = aero_weight_array_n_group(aero_state(1)%awa)
    n_class = aero_weight_array_n_class(aero_state(1)%awa)
    allocate(weight_type(nz,n_group,n_class))
    allocate(weight_magnitude(nz,n_group,n_class))
    allocate(weight_exponent(nz,n_group,n_class))
    do k = 1,nz
    do i_group = 1,n_group
    do i_class = 1,n_class
       weight_type(k,i_group,i_class) = &
            aero_state(k)%awa%weight(i_group,i_class)%type
       weight_magnitude(k,i_group,i_class) = &
            aero_state(k)%awa%weight(i_group,i_class)%magnitude
       weight_exponent(k,i_group,i_class) = &
            aero_state(k)%awa%weight(i_group,i_class)%exponent
    end do
    end do
    end do

    call pmc_nc_write_integer_3d(ncid, weight_type, &
         "weight_type", &
         (/ dimid_num_levels, dimid_aero_weight_group, dimid_aero_weight_class /), &
         description="type of each aerosol weighting function: 0 = invalid, " &
         // "1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), " &
         // "3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))")
    call pmc_nc_write_real_3d(ncid, weight_magnitude, &
         "weight_magnitude", &
         (/ dimid_num_levels, dimid_aero_weight_group, dimid_aero_weight_class /), &
         unit="m^{-3}", &
         description="magnitude for each weighting function")
    call pmc_nc_write_real_3d(ncid, weight_exponent, &
         "weight_exponent", &
         (/ dimid_num_levels, dimid_aero_weight_group, dimid_aero_weight_class /), unit="1", &
         description="exponent alpha for the power weight_type, " &
         // "set to -3 for MFA, and zero otherwise")

    call pmc_nc_write_integer64(ncid, next_id, "next_id", &
         unit="1", description="next_id for particle creation")

    call pmc_nc_check(nf90_close(ncid))
#endif

  end subroutine output_column_to_file_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Maps PartMC data structures to WRF data structure.
  subroutine partmc_to_wrf(grid, aero_states, aero_data, gas_states, &
       gas_data, pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke)

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
    !> Full domain of aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in):: aero_states
    type(aero_data_t), intent(in) :: aero_data
    !> Full domain of gas_states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(in) :: gas_states
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data

    integer :: i, j, k, i_class, i_group, n_class, i_spec
    real(kind=dp), allocatable :: num_concs(:)
    real(kind=dp) :: num_conc_class

    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
    do i_spec =1,gas_data_n_spec(gas_data)
       grid%chem(i,k,j,1+i_spec) = real(gas_states(i,k,j)%mix_rat(i_spec)/1000.0d0)
    end do
    end do
    end do
    end do

    n_class = aero_weight_array_n_class(aero_states(pmc_is,pmc_ks,pmc_js)%awa)
    !allocate(num_concs(aero_data_n_source(aero_data)))
    i_group = 1
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       do i_class = 1,n_class
          num_conc_class =  aero_state_total_particles( &
               aero_states(i,k,j), i_group, i_class) &
               * aero_states(i,k,j)%awa%weight(i_group, i_class)%magnitude
          grid%chem(i,k,j,p_NUM_CONC_a01 + i_class -1) = max(real(num_conc_class) &
               * grid%alt(i,k,j),1e-16)
       end do
    end do
    end do
    end do

  end subroutine partmc_to_wrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Maps WRF data to PartMC data structures.
  subroutine partmc_from_wrf(grid, gas_states, gas_data, pmc_is, pmc_ie, &
       pmc_js, pmc_je, pmc_ks, pmc_ke)

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
    !> Full domain of gas_states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: gas_states
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data

    integer :: i, j, k, i_spec

    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
    do i_spec = 1,gas_data_n_spec(gas_data)
       gas_states(i,k,j)%mix_rat(i_spec) = real(grid%chem(i,k,j,1+i_spec),kind=dp) * 1000.0d0
    end do
    end do
    end do
    end do


  end subroutine partmc_from_wrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a sampled sea salt emission using Gong 2003 parameterization.
  subroutine seasalt_emissions(wind_speed, aero_state, aero_data, env_state, &
       delta_t, allow_doubling, allow_halving, sea_salt_opt)

    !> Wind speed (m s^{-1}).
    real(kind=dp), intent(in) :: wind_speed
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environmental state.
    type(env_state_t), intent(in) :: env_state
    !> Timestep (s).
    real(kind=dp), intent(in) :: delta_t
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving
    !> Sea salt parameterization option.
    integer, intent(in) :: sea_salt_opt

    type(aero_dist_t) :: emissions
    real(kind=dp) :: p, emission_rate_scale
    real(kind=dp), allocatable :: sample_radius(:), sample_num_conc(:)
    type(bin_grid_t) :: radius_grid
    real(kind=dp), parameter :: d_min = 1d-9
    real(kind=dp), parameter :: d_max = 1e-5
    real(kind=dp) :: characteristic_factor
    real(kind=dp), parameter :: sample_timescale = 3600.0d0
    real(kind=dp) :: sample_timescale_effective
    integer, parameter :: n_bin = 25 
    integer :: i_bin
    character(len=AERO_SOURCE_NAME_LEN), parameter :: sea_salt = "sea_salt"

    real(kind=dp), parameter, dimension(5) :: sigma = [ 1.37,1.5,1.42,1.53,1.85]
    real(kind=dp), parameter, dimension(5) :: CMD = [.018, .041, .09, .23, .83]
    real(kind=dp), parameter, dimension(5) :: coef = [104.5, .0442, 149.6, 2.96, .51]
    real(kind=dp), parameter, dimension(5) :: offset = [1d5,1d5,1d5,1d5,2d5]
    real(kind=dp), parameter, dimension(5) :: expon = [.556,1.08,.545,.79,.87]
    character(len=AERO_SOURCE_NAME_LEN) :: mode_name
    integer :: i_mode, n_modes
    real(kind=dp) :: Re, F_i, H_s, vu_w, C_d
    real(kind=dp), allocatable :: num_concs(:), ss_masses(:)
    real(kind=dp) :: num_conc_before, ss_before
    integer :: i_ss

    ! Gong
    integer, parameter :: Gong_SS = 1
    ! Ovadnevaite 2014 doi:10.5194/acp-14-1837-2014
    integer, parameter :: Ovadnevaite_SS = 2
    ! 
    logical, parameter :: match_mam3_modes = .true.
    real(kind=dp), parameter, dimension(2) :: mam3_num_conc_flux = [73.7095947,10.1706810]
    real(kind=dp), parameter, dimension(2) :: mam3_diams = [2.03347408e-07, &
         1.09840811e-06]
    real(kind=dp), parameter, dimension(2) :: mam3_std = [2.0,2.2]

    real(kind=dp), parameter :: mw_na = 22.990d0
    real(kind=dp), parameter :: mw_cl = 35.450d0

    if (sea_salt_opt == Gong_SS) then
       if (match_mam3_modes) then
          n_modes = 2
          allocate(emissions%mode(n_modes))
          do i_mode =1,n_modes
             emissions%mode(i_mode)%type = AERO_MODE_TYPE_LOG_NORMAL
             allocate(emissions%mode(i_mode)%vol_frac(aero_data_n_spec(aero_data)))
             allocate(emissions%mode(i_mode)%vol_frac_std(aero_data_n_spec(aero_data)))
             emissions%mode(i_mode)%vol_frac = 0d0
             emissions%mode(i_mode)%vol_frac_std = 0d0
             emissions%mode(i_mode)%vol_frac(aero_data_spec_by_name(aero_data, "Na")) = &
                  mw_na / (mw_na + mw_cl)
             emissions%mode(i_mode)%vol_frac(aero_data_spec_by_name(aero_data, "Cl")) = &
                  1.0 - (mw_na / (mw_na + mw_cl))
             write(mode_name,'(a,i2.2)') 'sea_salt_', i_mode
             emissions%mode(i_mode)%name = mode_name
             emissions%mode(i_mode)%source = string_array_find(aero_data%source_name, &
                  mode_name)
             emissions%mode(i_mode)%weight_class = string_array_find( &
                  aero_data%weight_class_name, mode_name)
             emissions%mode(i_mode)%char_radius =  mam3_diams(i_mode) / 2
             emissions%mode(i_mode)%log10_std_dev_radius = log10(mam3_std(i_mode))
             emissions%mode(i_mode)%num_conc = mam3_num_conc_flux(i_mode) 
          end do
          emission_rate_scale =  wind_speed**3.41
       else
          i_ss = 1
          call bin_grid_make(radius_grid, BIN_GRID_TYPE_LOG, n_bin, d_min, d_max)
          allocate(emissions%mode(1))
          emissions%mode(1)%type = AERO_MODE_TYPE_SAMPLED
          allocate(emissions%mode(1)%sample_radius(n_bin + 1))
          emissions%mode(1)%sample_radius = radius_grid%edges
          allocate(emissions%mode(1)%sample_num_conc(n_bin))
          do i_bin = 1, n_bin
             call seasalt_emissions_1bin(emissions%mode(1)%sample_radius(i_bin), &
                  emissions%mode(1)%sample_radius(i_bin+1), &
                  emissions%mode(1)%sample_num_conc(i_bin))
          end do
          allocate(emissions%mode(1)%vol_frac(aero_data_n_spec(aero_data)))
          allocate(emissions%mode(1)%vol_frac_std(aero_data_n_spec(aero_data)))
          emissions%mode(1)%vol_frac = 0d0
          emissions%mode(1)%vol_frac_std = 0d0
          emissions%mode(1)%vol_frac(aero_data_spec_by_name(aero_data, "Na")) = &
                  mw_na / (mw_na + mw_cl)
          emissions%mode(1)%vol_frac(aero_data_spec_by_name(aero_data, "Cl")) = &
                  1.0 - (mw_na / (mw_na + mw_cl))
          write(mode_name,'(a,i2.2)') 'sea_salt_', i_ss
          emissions%mode(1)%name = sea_salt
          emissions%mode(1)%source = string_array_find(aero_data%source_name, &
               mode_name) 
          emissions%mode(1)%weight_class = string_array_find( &
               aero_data%weight_class_name, mode_name)
          emission_rate_scale = wind_speed**3.41
       end if
    else if (sea_salt_opt == Ovadnevaite_SS) then
       n_modes = 5
       allocate(emissions%mode(n_modes))
       do i_mode =1,n_modes
          emissions%mode(i_mode)%type = AERO_MODE_TYPE_LOG_NORMAL
          allocate(emissions%mode(i_mode)%vol_frac(aero_data_n_spec(aero_data)))
          allocate(emissions%mode(i_mode)%vol_frac_std(aero_data_n_spec(aero_data)))
          emissions%mode(i_mode)%vol_frac = 0d0
          emissions%mode(i_mode)%vol_frac_std = 0d0
          emissions%mode(i_mode)%vol_frac(aero_data_spec_by_name(aero_data, "Na")) = &
               mw_na / (mw_na + mw_cl)
          emissions%mode(i_mode)%vol_frac(aero_data_spec_by_name(aero_data, "Cl")) = &
                  1.0 - (mw_na / (mw_na + mw_cl))
          write(mode_name,'(a,i2.2)') 'sea_salt_', i_mode
          emissions%mode(i_mode)%name = mode_name
          emissions%mode(i_mode)%source = string_array_find(aero_data%source_name, &
               mode_name)
          emissions%mode(i_mode)%weight_class = string_array_find( &
               aero_data%weight_class_name, mode_name)
          H_s = 1.0d0
          vu_w = 1.0d-6
          C_d = 1.0d-3
          emissions%mode(i_mode)%char_radius =  .5*CMD(i_mode) * 1d-6
          emissions%mode(i_mode)%log10_std_dev_radius = log10(sigma(i_mode))
          Re = C_d**.5 * wind_speed * H_s / vu_w
          Re = 3.1d6
          F_i = coef(i_mode) * (Re - offset(i_mode))**expon(i_mode)
          emissions%mode(i_mode)%num_conc = F_i
       end do
       emission_rate_scale = 1.0d0
    else
       print*, 'invalid sea salt option'
    end if

    p = emission_rate_scale * delta_t / env_state%height

    sample_timescale_effective = max(1.0, min(sample_timescale, &
         env_state%elapsed_time))
    characteristic_factor = sample_timescale_effective / delta_t

    call aero_state_add_aero_dist_sample(aero_state, aero_data, &
         emissions, p, characteristic_factor, env_state%elapsed_time, &
         allow_doubling, allow_halving)

  end subroutine seasalt_emissions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the integrated total number for a given size bin.
  subroutine seasalt_emissions_1bin(dry_rad_lo, dry_rad_hi, numb_flux)

    !> Lower particle radius (m).
    real(kind=dp), intent(in) :: dry_rad_lo
    !> Upper particle radius (m).
    real(kind=dp), intent(in) :: dry_rad_hi
    !> Number flux (# m^{-2} s^{-1}).
    real(kind=dp), intent(out) :: numb_flux

    real(kind=dp), parameter :: c1 = 0.7674
    real(kind=dp), parameter :: c2 = 3.079
    real(kind=dp), parameter :: c3 = 2.573e-11
    real(kind=dp), parameter :: c4 = -1.424
    real(kind=dp), parameter :: onethird = 1.0/3.0
    real(kind=dp), parameter :: pi = 3.14
    real(kind=dp) :: relhum
    integer, parameter :: n_bin = 10
    integer :: i_bin
    real(kind=dp) :: rdryhi, rdrylo, rdryaa, rdrybb, rdry
    real(kind=dp) :: rwet, rwetaa, rwetbb
    real(kind=dp) :: sum_numb
    real(kind=dp) :: dum, duma, dumb
    real(kind=dp) :: rwet_cm, rdry_cm
    real(kind=dp) :: df0drwet, drwet
    real(kind=dp) :: alnrdrylo, dlnrdry

    ! Assumed relative humidity
    relhum = 0.80d0

    sum_numb = 0.0d0

    ! convert from m to microns
    rdrylo = dry_rad_lo * 1d6
    rdryhi = dry_rad_hi * 1d6

    alnrdrylo = log( rdrylo )
    dlnrdry = (log( rdryhi ) - alnrdrylo)/n_bin

    rdrybb = exp( alnrdrylo )
    rdry_cm = rdrybb*1.0d-4
    rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
    rwetbb = rwet_cm*1.0d4

    do i_bin = 1,n_bin
       rdryaa = rdrybb
       rwetaa = rwetbb
       dum = alnrdrylo + i_bin*dlnrdry
       rdrybb = exp( dum )
       rdry_cm = rdrybb*1.0d-4
       rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
          ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
       rwetbb = rwet_cm*1.0d4

       rdry = sqrt(rdryaa * rdrybb)
       rwet = sqrt(rwetaa * rwetbb)
       drwet = rwetbb - rwetaa

       dumb = (0.433 - log10(rwet)) / 0.433 !( 0.380 - log10(rwet) ) / 0.650
       duma = 4.7*(1.0 + 30.0* rwet)**(-0.017*rwet**(-1.44))
       df0drwet = 1.373 * (rwet**(-duma)) * &
            (1.0 + 0.057*(rwet**3.45)) * &
            (10.0**(1.607*exp( -dumb*dumb)))
       sum_numb = sum_numb + drwet*df0drwet
    end do

    numb_flux = sum_numb

  end subroutine seasalt_emissions_1bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the total grid cell volume.
  real(kind=dp) function get_grid_cell_volume(grid, i, k, j)
    !> WRF grid.
    type(domain) :: grid
    !> Grid cell index west-east.
    integer :: i
    !> Grid cell index bottom-top.
    integer :: k
    !> Grid cell index south-north.
    integer :: j

    real(kind=dp) :: dx, dy, dz

    dz = real(grid%z_at_w(i,k+1,j) - grid%z_at_w(i,k,j), kind=dp)
    dx = real((1.0d0/grid%rdx)/grid%msftx(i,j), kind=dp)
    dy = real((1.0d0/grid%rdy)/grid%msfty(i,j), kind=dp)

    get_grid_cell_volume = dx*dy*dz

  end function get_grid_cell_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the bulk aerosol optical properties that WRF needs.
  subroutine compute_bulk_optical_props(grid, aero_data, aero_state, &
    env_state)

    !> WRF domain.
    type(domain), intent(inout) :: grid
    !> Full domain of aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Environmental state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), allocatable :: g_top(:), omega(:), tau(:), g_bottom(:)
    real(kind=dp), allocatable :: num_concs(:)
    real(kind=dp) :: B_scat, B_abs, B_ext
    integer :: i_part, i_lambda, i, j, k, n_swbands, idx

    n_swbands = 5 !size(aero_state%apa%particle(1)%absorb_cross_sect)

    allocate(g_top(n_swbands))
    allocate(omega(n_swbands))
    allocate(tau(n_swbands))
    allocate(g_bottom(n_swbands))
    g_top = 0.0d0
    g_bottom = 0.0d0
    omega = 0.0d0
    tau = 0.0d0

    num_concs = aero_state_num_concs(aero_state, aero_data)
    do i_lambda = 1,n_swbands
       B_scat = 0.0d0
       B_abs = 0.0d0
       B_ext = 0.0d0
       do i_part = 1,aero_state_n_part(aero_state)
          g_top(i_lambda) = g_top(i_lambda) + num_concs(i_part) &
               * aero_state%apa%particle(i_part)%asymmetry(i_lambda) &
               * aero_state%apa%particle(i_part)%scatter_cross_sect(i_lambda)
          g_bottom(i_lambda) = g_bottom(i_lambda) + num_concs(i_part) &
               *  aero_state%apa%particle(i_part)%scatter_cross_sect(i_lambda)
          B_scat = B_scat + num_concs(i_part) &
               * aero_state%apa%particle(i_part)%scatter_cross_sect(i_lambda)
          B_abs = B_abs +  num_concs(i_part) &
               * aero_state%apa%particle(i_part)%absorb_cross_sect(i_lambda)
       end do
       omega(i_lambda) = B_scat / (B_scat + B_abs)
       tau(i_lambda) = (B_scat + B_abs) * env_state%height
    end do

    i = env_state%ix
    j = env_state%iy
    k = env_state%iz

    idx = 1
    grid%tauaer1(i,k,j) = tau(idx)
    grid%gaer1(i,k,j) = g_top(idx) / g_bottom(idx)
    grid%waer1(i,k,j) = omega(idx)

    idx = 2
    grid%tauaer2(i,k,j) = tau(idx)
    grid%gaer2(i,k,j) = g_top(idx) / g_bottom(idx)
    grid%waer2(i,k,j) = omega(idx)

    idx = 4 
    grid%tauaer3(i,k,j) = tau(idx)
    grid%gaer3(i,k,j) = g_top(idx) / g_bottom(idx)
    grid%waer3(i,k,j) = omega(idx)

    idx = 5
    grid%tauaer4(i,k,j) = tau(idx)
    grid%gaer4(i,k,j) = g_top(idx) / g_bottom(idx)
    grid%waer4(i,k,j) = omega(idx)

  end subroutine compute_bulk_optical_props

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module wrf_pmc_driver
