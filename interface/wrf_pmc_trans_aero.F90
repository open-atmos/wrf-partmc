! Copyright (C) 2011-2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
! The partmc_trans_aero module.

! Subroutines associated with aerosol transport
module wrf_pmc_trans_aero

  ! PartMC modules
  use pmc_aero_state
  use pmc_aero_data
  use pmc_scenario
  use pmc_env_state
  use pmc_util
  use pmc_rand
  use pmc_mpi
  use pmc_integer_varray

  use module_wrf_error
  use module_domain_type
  use module_domain
  use module_model_constants
  use module_configure, only: p_qv
  USE module_state_description
  USE module_domain, ONLY : domain, domain_clock_get
  USE module_configure, ONLY : grid_config_rec_type
  ! Get the WRF MPI grid
  USE module_dm, ONLY: local_communicator_periodic, local_communicator
#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines all the particle sets to move due to transport.
  subroutine trans_aero(grid, config_flags, aero_data, scenario, env_states, &
       aero_states, pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, &
       global_nx, global_ny, global_nz)

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
    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> WRF namelist configuration.
    type(grid_config_rec_type), intent(in) :: config_flags
    !> Full domain of aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Scenario data.
    type(scenario_t),dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenario
    !> Full domain of environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Full domain of aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: aero_states
     !> Number of grid cells for entire domain in west-east direction.
    integer, intent(in) :: global_nx
    !> Number of grid cells for entire domain in north-south direction.
    integer, intent(in) :: global_ny
    !> Number of grid cells for entire domain in vertical direction.
    integer, intent(in) :: global_nz

    ! Local variables
    type(aero_state_t), allocatable :: delta_aero_states(:,:,:)

    integer :: i, j, k
    integer :: is, ie, js, je, ks, ke
    integer :: n_parts_before_add, n_parts_after_add,n_parts_to_be_added
    real(kind=dp) :: t1, t2

    ! MPI related
    integer :: dir, disp
    integer :: i_send_s, i_send_e, i_recv_s, i_recv_e
    integer :: j_send_s, j_send_e, j_recv_s, j_recv_e
    integer :: k_send_s, k_send_e, k_recv_s, k_recv_e
    integer :: source, dest, ierr

    character, dimension(:), allocatable :: send_buffer, recv_buffer
    integer :: max_size
    integer :: max_size_send, max_size_recv
    integer :: request, send_request, request_size, send_request_size
    integer :: tag, tag_buffer_size
    integer :: n_particles, n_particles_delta
    integer :: i_inc, j_inc, k_inc
    integer :: position
    integer :: status(MPI_STATUS_SIZE)
    type(aero_state_t) :: unpack_aero_state

    logical :: ideal
    logical :: throw_away
    integer :: total_particles, particle_send_count
    integer :: n_part_add, i_group, i_class
    real(kind=dp) :: prob, total_num_conc
    real(kind=dp), allocatable, dimension(:,:,:,:,:) :: expected_num_conc
    integer :: max_class, max_group
    integer :: i_start, i_end, j_start, j_end
    real(kind=dp) :: dest_vol, source_vol

    type(aero_dist_t) :: background
    real(kind=dp) :: dilution_rate 
    integer :: p_t, p_tm1
    logical :: change_boundary_conditions

    character, dimension(:), allocatable :: recv_buffer_cell_volume, &
       send_buffer_cell_volume, recv_buffer_rrho, send_buffer_rrho
    integer :: i_pos, j_pos, k_pos
    real(kind=dp), allocatable, dimension(:,:,:) :: send_array_cell_volume, &
       send_array_rrho, recv_cell_volume, recv_rrho
    integer :: tag_cv, tag_rrho
    integer :: cell_volume_send_request, rrho_send_request, &
         cell_volume_recv_request, rrho_recv_request
    type(env_state_t), allocatable, dimension(:,:,:) :: env_states_temp

    if (config_flags%periodic_x .and. config_flags%periodic_y) then
       ideal = .true.
    else
       ideal = .false.
    end if

    ! Determine the extents
    call get_array_extents(pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
         global_nx, global_ny, global_nz, is, ie, js, je, ks, ke, config_flags)

    allocate(env_states_temp(is:ie,ks:ke,js:je))
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       env_states_temp(i,k,j) = env_states(i,k,j)
    end do
    end do
    end do

    tag_cv = 11
    tag_rrho = 12
    do dir = 1,2
    do disp = -1,1,2

       if (ideal) then
          call MPI_cart_shift(local_communicator_periodic, (dir - 1), disp, &
               source, dest, ierr)
       else
          call MPI_cart_shift(local_communicator, (dir - 1), disp, &
               source, dest, ierr)
       end if

       call get_array_index_for_send_recv_to_ghost(disp, dir, is, ie, js, je, &
            ks, ke, i_send_s, i_send_e, j_send_s, j_send_e, k_send_s, &
            k_send_e, i_recv_s, i_recv_e, j_recv_s, j_recv_e, k_recv_s, &
            k_recv_e)
!       allocate(send_array_cell_volume(i_send_s:i_send_e,k_send_s:k_send_e, &
!            j_send_s:j_send_e))
!       allocate(send_array_rrho(i_send_s:i_send_e,k_send_s:k_send_e, &
!            j_send_s:j_send_e))
       allocate(send_array_cell_volume(i_send_e - i_send_s + 1, &
            k_send_e - k_send_s + 1, j_send_e - j_send_s + 1))
       allocate(send_array_rrho(i_send_e - i_send_s + 1, &
            k_send_e - k_send_s + 1, j_send_e - j_send_s + 1))

       send_array_rrho = 0.0d0
       send_array_cell_volume = 0.0d0

       max_size = pmc_mpi_pack_size_real_array_3d(send_array_cell_volume)
       if (dest /= MPI_PROC_NULL) then
       do j_inc = j_send_s, j_send_e
       do k_inc = k_send_s, k_send_e
       do i_inc = i_send_s, i_send_e
          i_pos = i_inc - i_send_s + 1
          j_pos = j_inc - j_send_s + 1
          k_pos = k_inc - k_send_s + 1
          send_array_cell_volume(i_pos,k_pos,j_pos) = &
               env_states(i_inc,k_inc,j_inc)%cell_volume
          send_array_rrho(i_pos,k_pos,j_pos) = &
               env_states(i_inc,k_inc,j_inc)%rrho
       end do
       end do
       end do
       end if

       ! Send and receive
       if (dest /= MPI_PROC_NULL) then
          ! Pack the buffer
          allocate(send_buffer_cell_volume(max_size))
          allocate(send_buffer_rrho(max_size))
          position = 0
          call pmc_mpi_pack_real_array_3d(send_buffer_cell_volume, position, &
               send_array_cell_volume)
          position = 0
          call pmc_mpi_pack_real_array_3d(send_buffer_rrho, position, &
               send_array_rrho)

          if (ideal) then
             call MPI_isend(send_buffer_cell_volume, max_size, MPI_PACKED, dest, &
                  tag_cv, local_communicator_periodic, cell_volume_send_request, ierr)
             call MPI_isend(send_buffer_rrho, max_size, MPI_PACKED, dest, &
                  tag_rrho, local_communicator_periodic, rrho_send_request, ierr)
          else
             call MPI_isend(send_buffer_cell_volume, max_size, MPI_PACKED, dest, &
                  tag_cv, local_communicator, cell_volume_send_request, ierr)
             call MPI_isend(send_buffer_rrho, max_size, MPI_PACKED, dest, &
                  tag_rrho, local_communicator, rrho_send_request, ierr)
          end if

       end if

       if (source /= MPI_PROC_NULL) then
          allocate(recv_buffer_cell_volume(max_size))
          allocate(recv_buffer_rrho(max_size))
          ! Non-blocking receiving of buffer from source
          if (ideal) then
             call MPI_Irecv(recv_buffer_cell_volume, max_size, MPI_PACKED, source, &
                  tag_cv, local_communicator_periodic, cell_volume_recv_request, ierr)
             call MPI_Irecv(recv_buffer_rrho, max_size, MPI_PACKED, source, &
                  tag_rrho, local_communicator_periodic, rrho_recv_request, ierr)
          else
             call MPI_Irecv(recv_buffer_cell_volume, max_size, MPI_PACKED, source, &
                  tag_cv, local_communicator, cell_volume_recv_request, ierr)
             call MPI_Irecv(recv_buffer_rrho, max_size, MPI_PACKED, source, &
                  tag_rrho, local_communicator, rrho_recv_request, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          call MPI_wait(cell_volume_send_request, status, ierr)
          call MPI_wait(rrho_send_request, status, ierr)
       end if

       if (source /= MPI_PROC_NULL) then
          ! Wait
          call MPI_wait(cell_volume_recv_request, status, ierr)
          call MPI_wait(rrho_recv_request, status, ierr)
          ! Unpack
          position = 0
          call pmc_mpi_unpack_real_array_3d(recv_buffer_cell_volume, position, &
               recv_cell_volume)
          position = 0
          call pmc_mpi_unpack_real_array_3d(recv_buffer_rrho, position, &
               recv_rrho)
          do j_inc = j_recv_s, j_recv_e
          do k_inc = k_recv_s, k_recv_e
          do i_inc = i_recv_s, i_recv_e
             i_pos = i_inc - i_recv_s + 1
             j_pos = j_inc - j_recv_s + 1
             k_pos = k_inc - k_recv_s + 1
             env_states_temp(i_inc,k_inc,j_inc)%cell_volume = recv_cell_volume( &
                  i_pos,k_pos,j_pos)
             env_states_temp(i_inc,k_inc,j_inc)%rrho = recv_rrho(i_pos,k_pos,j_pos)
          end do
          end do
          end do
       end if
       if (allocated(recv_buffer_cell_volume)) deallocate(recv_buffer_cell_volume)
       if (allocated(recv_buffer_rrho)) deallocate(recv_buffer_rrho)
       if (allocated(send_buffer_cell_volume)) deallocate(send_buffer_cell_volume)
       if (allocated(send_buffer_rrho)) deallocate(send_buffer_rrho)
       if (allocated(send_array_cell_volume)) deallocate(send_array_cell_volume)
       if (allocated(send_array_rrho)) deallocate(send_array_rrho)
    end do
    end do

#ifdef PMC_DEBUG
    total_particles = 0
    total_num_conc = 0.0d0
    ! Determine total particles
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       total_particles = total_particles + &
            aero_state_total_particles(aero_states(i,k,j))
       total_num_conc = total_num_conc + &
            aero_state_total_num_conc(aero_states(i,k,j), aero_data)
    end do
    end do
    end do
    print "(A50,I10,E16.7)", 'Total particles before preweighting = ', &
         total_particles, total_num_conc
#endif

    ! We want to precompute transport to adjust weightings in advance
    call trans_aero_preweight(grid, config_flags, aero_data, scenario, env_states_temp, &
         aero_states, is, ie, js, je, ks, ke, &
         pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, global_nx, &
         global_ny, global_nz, expected_num_conc)
    call pmc_mpi_barrier()
#ifdef PMC_DEBUG
    total_particles = 0
    total_num_conc = 0.0d0
    ! Determine total particles
    do j = pmc_js,pmc_je !0, ny + 1
    do k = pmc_ks,pmc_ke !1, nz
    do i = pmc_is,pmc_ie !0, nx + 1
       total_particles = total_particles + &
            aero_state_total_particles(aero_states(i,k,j))
       total_num_conc = total_num_conc + &
            aero_state_total_num_conc(aero_states(i,k,j), aero_data)
    end do
    end do
    end do
    print "(A50,I10,E16.7)", 'Total particles after preweighting = ', &
        total_particles, total_num_conc
#endif

#ifdef PMC_DEBUG
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
    do i_class = 1, min(3,aero_weight_array_n_class(aero_states(i,k,j)%awa))
      print*, 'Status:', i,k,j, i_class, &
           aero_state_total_particles(aero_states(i,k,j),1,i_class), &
           'expected', expected_num_conc(i,k,j,1,i_class), &
           'actual', aero_state_group_class_num_conc(aero_states(i,k,j), &
           aero_data, 1, i_class)
    end do
    end do
    end do
    end do
#endif

    ! Overall workflow strategy
    ! 1.) Allocate these delta_aero_states to hold the transition particles
    ! 2.) Allocate the computational volumes of the delta_aero_states
    !      and communicate to get the correct ghost point comp_vols
    ! 3.) Loop over all points and sample into the correct delta_states
    !      using a multinomial sampling approach.
    ! 3.) Send the delta_states that are in the ghost points to where they
    !      should be on the neighboring domain
    ! 4.) Add the particles to the right spots 
    ! 5.) Deallocate all those delta_aero_states

    ! Allocate the delta_aero_states (with ghost points)
   
!    ! Determine the extents
!    call get_array_extents(pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
!         global_nx, global_ny, global_nz, is, ie, js, je, ks, ke, config_flags)

#ifdef PMC_DEBUG
    print*, 'physical particles are in cells:', pmc_is,pmc_ie, pmc_js, pmc_je, &
         pmc_ks, pmc_ke
    print*, 'This processor has transport cells:', is,ie,js,je,ks,ke
#endif

    allocate(delta_aero_states(is:ie,ks:ke,js:je))
    do j = js,je
    do k = ks,ke
    do i = is,ie
       ! FIXME: Do we really need aero_data from elsewhere. We shouldn't.
       ! aero_data should be constant. it is needed for the number of sources.
       call aero_state_zero(delta_aero_states(i,k,j))
       call aero_state_set_weight(delta_aero_states(i,k,j), &
            aero_data, AERO_STATE_WEIGHT_FLAT_SPECIFIED)
       call aero_state_set_n_part_ideal(delta_aero_states(i,k,j), &
            real(grid%num_particles,kind=dp))
    end do
    end do
    end do

    ! Copy the local computational volume information
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
        delta_aero_states(i,k,j)%awa = aero_states(i,k,j)%awa
    end do
    end do
    end do

    ! Get ghost point comp_vols
    tag = 0
    do dir = 1,2
    do disp = -1,1,2
       ! Set the loop variables
       call get_array_index_for_send_recv_to_ghost(disp, dir, is, ie, js, je, &
            ks, ke, i_send_s, i_send_e, j_send_s, j_send_e, k_send_s, &
            k_send_e, i_recv_s, i_recv_e, j_recv_s, j_recv_e, k_recv_s, &
            k_recv_e)
       if (dir == 1) then ! Moving in the y direction
          max_size = 0
          do i_inc = i_recv_s, i_recv_e
          do k_inc = k_recv_s, k_recv_e
               max_size = max_size + pmc_mpi_pack_size_aero_weight_array( &
               aero_states(i_inc,k_inc,j_send_s)%awa)
          end do
          end do
       elseif (dir == 2) then ! Moving in the x direction
          max_size = 0
          do j_inc = j_recv_s, j_recv_e
          do k_inc = k_recv_s, k_recv_e
             max_size = max_size + pmc_mpi_pack_size_aero_weight_array( & 
               aero_states(i_send_s,k_inc,j_inc)%awa)
          end do
          end do
       end if

       if (ideal) then
          call MPI_cart_shift(local_communicator_periodic, (dir - 1), disp, &
               source, dest, ierr)
       else
          call MPI_cart_shift(local_communicator, (dir - 1), disp, &
               source, dest, ierr)
       end if

#ifdef PMC_DEBUG
       print*, 'Rank: ',pmc_mpi_rank(), 'source: ',source, 'dest', dest, &
           'dir', dir-1, 'disp',disp
#endif
       ! Allocate a send buffer and a receive buffer
       allocate(send_buffer(max_size))
       allocate(recv_buffer(max_size))

       if (source /= MPI_PROC_NULL) then
          ! Non-blocking receiving of buffer from source
          if (ideal) then
             call MPI_Irecv(recv_buffer, max_size, MPI_PACKED, source, tag, &
                  local_communicator_periodic, request, ierr)
          else
             call MPI_Irecv(recv_buffer, max_size, MPI_PACKED, source, tag, &
                  local_communicator, request, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          ! Pack the buffer with aero weights
          position = 0
          do j_inc = j_send_s, j_send_e
          do k_inc = k_send_s, k_send_e
          do i_inc = i_send_s, i_send_e
              call pmc_mpi_pack_aero_weight_array(send_buffer, position, &
                   aero_states(i_inc, k_inc, j_inc)%awa)
          end do
          end do
          end do
          ! Send the buffer to dest
          if (ideal) then
             call MPI_send(send_buffer, max_size, MPI_PACKED, dest, tag, &
                  local_communicator_periodic, ierr)
          else
             call MPI_send(send_buffer, max_size, MPI_PACKED, dest, tag, &
                  local_communicator, ierr)
          end if
       end if

       if (source /= MPI_PROC_NULL) then
          ! Wait
          call MPI_wait(request, status, ierr)
          ! Unpack the aero_weights
          position = 0
          do j_inc = j_recv_s, j_recv_e
          do k_inc = k_recv_s, k_recv_e
          do i_inc = i_recv_s, i_recv_e
              call pmc_mpi_unpack_aero_weight_array(recv_buffer, position, &
                   delta_aero_states(i_inc,k_inc,j_inc)%awa) 
          end do
          end do
          end do
       end if

       if (allocated(send_buffer)) deallocate(send_buffer)
       if (allocated(recv_buffer)) deallocate(recv_buffer)
    end do
    end do

    call pmc_mpi_barrier()

!    total_particles = 0
!    ! Determine total particles
!    do j = pmc_js,pmc_je
!    do k = pmc_ks,pmc_ke
!    do i = pmc_is,pmc_ie
!       total_particles = total_particles + &
!            aero_state_total_particles(aero_states(i,k,j))
!    end do
!    end do
!    end do
!    print*, 'Number of particles before sampling', total_particles

!    call aero_states_debug(aero_states, pmc_is, pmc_ie, pmc_ks, pmc_ke, &
!         pmc_js, pmc_je)

#ifdef PMC_DEBUG
    print*, 'start sampling particles'
    t1 =  MPI_Wtime()
#endif

    if (config_flags%do_uniform) then
       i_start = pmc_is
       i_end = pmc_ie
       j_start = pmc_js
       j_end = pmc_je
    else if (config_flags%do_scm) then
       i_start = pmc_is
       i_end = pmc_ie !global_nx
       j_start = pmc_js
       j_end = pmc_je !global_ny
    else if (config_flags%periodic_x .and. config_flags%periodic_y) then
       i_start = pmc_is
       i_end = pmc_ie !global_nx
       j_start = pmc_js
       j_end = pmc_je !global_ny
    else
       i_start = max(pmc_is,2)
       i_end = min(pmc_ie,global_nx-1)
       j_start = max(pmc_js,2)
       j_end = min(pmc_je,global_ny-1)
    end if

    do j = j_start,j_end
    do k = pmc_ks,pmc_ke
    do i = i_start,i_end
       call aero_state_multisample(env_states_temp, aero_data, &
            aero_states(i,k,j), &
            delta_aero_states,is,ie,js,je,ks,ke, &
            pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, &
            i,j,k)
    end do
    end do
    end do

#ifdef PMC_DEBUG
    t2 = MPI_Wtime()
    print*, 'done sampling particles. time spent = ', t2-t1
#endif
    call pmc_mpi_barrier()

#ifdef PMC_DEBUG
    print*, 'start sending particles'
#endif

    t1 = MPI_Wtime()
    ! Send the delta_states that are in ghost points
    tag = 1
    tag_buffer_size = 2
    do dir = 1,2
    do disp = -1,1,2
       call get_array_index_for_send_recv_from_ghost(disp, dir, is, ie, js, je, &
            ks, ke, i_send_s, i_send_e, j_send_s, j_send_e, k_send_s, &
            k_send_e, i_recv_s, i_recv_e, j_recv_s, j_recv_e, k_recv_s, &
            k_recv_e)

       ! Buffer size to send all the particles
       max_size_send = 0
       max_size_recv = 0
       particle_send_count = 0
       do j_inc = j_send_s, j_send_e
       do k_inc = k_send_s, k_send_e
       do i_inc = i_send_s, i_send_e
          max_size_send = max_size_send + pmc_mpi_pack_size_aero_state( &
             delta_aero_states(i_inc,k_inc,j_inc))
          particle_send_count = particle_send_count +  &
               aero_state_total_particles(delta_aero_states(i_inc,k_inc,j_inc))
       end do
       end do
       end do

       if (ideal) then
          call MPI_cart_shift(local_communicator_periodic, (dir - 1), disp, &
               source, dest, ierr)
       else
          call MPI_cart_shift(local_communicator, (dir - 1), disp, &
               source, dest, ierr)
       end if

       ! Allocate a send buffer
       allocate(send_buffer(max_size_send))
       ! Get the size of the receive buffer
       if (source /= MPI_PROC_NULL) then
          if (ideal) then
             call MPI_Irecv(max_size_recv, 1, MPI_INTEGER, source, tag_buffer_size, &
                  local_communicator_periodic, request_size, ierr)
          else
             call MPI_Irecv(max_size_recv, 1, MPI_INTEGER, source, tag_buffer_size, &
                  local_communicator, request_size, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          if (ideal) then
             call MPI_isend(max_size_send, 1, MPI_INTEGER, dest, tag_buffer_size, &
                  local_communicator_periodic, send_request_size, ierr)
          else
             call MPI_isend(max_size_send, 1, MPI_INTEGER, dest, tag_buffer_size, &
                  local_communicator, send_request_size, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          call MPI_wait(send_request_size, status, ierr)
       end if

       if (source /= MPI_PROC_NULL) then
          call MPI_wait(request_size, status, ierr)
          allocate(recv_buffer(max_size_recv))
       end if
#ifdef PMC_DEBUG
       print*, 'sending', particle_send_count, 'particles to ', dest, max_size_send, &
                max_size_recv
#endif

       ! Advance to get the actual particles
       if (source /= MPI_PROC_NULL) then
          ! Non-blocking receiving of buffer from source
          if (ideal) then
             call MPI_Irecv(recv_buffer, max_size_recv, MPI_PACKED, source, tag, &
                  local_communicator_periodic, request, ierr)
          else
             call MPI_Irecv(recv_buffer, max_size_recv, MPI_PACKED, source, tag, &
                  local_communicator, request, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          ! Pack the buffer
          position = 0
          do j_inc = j_send_s, j_send_e
          do k_inc = k_send_s, k_send_e
          do i_inc = i_send_s, i_send_e
             call pmc_mpi_pack_aero_state(send_buffer, position, &
                  delta_aero_states(i_inc, k_inc, j_inc))
          end do
          end do
          end do
          do j_inc = j_send_s, j_send_e
          do k_inc = k_send_s, k_send_e
          do i_inc = i_send_s, i_send_e
             call aero_state_zero(delta_aero_states(i_inc, k_inc, j_inc))
          end do
          end do
          end do
          ! Send the buffer to dest
          if (ideal) then
             call MPI_isend(send_buffer, max_size_send, MPI_PACKED, dest, tag, &
                  local_communicator_periodic, send_request, ierr)
          else
             call MPI_isend(send_buffer, max_size_send, MPI_PACKED, dest, tag, &
                  local_communicator, send_request, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          call MPI_wait(send_request, status, ierr)
       end if

       if (source /= MPI_PROC_NULL) then
          ! Wait
          call MPI_wait(request, status, ierr)
          position = 0
          do j_inc = j_recv_s,j_recv_e
          do k_inc = k_recv_s,k_recv_e
          do i_inc = i_recv_s,i_recv_e
             call pmc_mpi_unpack_aero_state(recv_buffer, position, &
                  unpack_aero_state)
             call aero_state_add_particles( &
                  delta_aero_states(i_inc,k_inc,j_inc), unpack_aero_state, &
                  aero_data)
             call aero_state_zero(unpack_aero_state)
          end do
          end do
          end do
       end if
       if (allocated(send_buffer)) deallocate(send_buffer)
       if (allocated(recv_buffer)) deallocate(recv_buffer)
    end do
    end do

#ifdef PMC_DEBUG
    t2 = MPI_Wtime()
    print*, 'done sending particles. time spent = ', t2-t1
#endif

    call pmc_mpi_barrier()

#ifdef PMC_DEBUG
    print*, 'start boundary conditions'
#endif
    t1 = MPI_Wtime()
    !print*, pmc_is, pmc_ie, pmc_js, pmc_je, global_nx, global_ny
    if (.not. config_flags%periodic_x) then
       ! FIXME: Fix the background to change over time
       if (pmc_is == 1) then
          i = 1
          do j = max(pmc_js,2), min(pmc_je,global_ny-1)
             do k = pmc_ks, pmc_ke
             dest_vol = env_states(i+1,k,j)%cell_volume
             source_vol = env_states(i,k,j)%cell_volume
             do i_class = 1,aero_weight_array_n_class(aero_states(i,k,j)%awa)
                 prob = min(env_states(i,k,j)%prob_advection(1,0,0,i_class) + &
                      env_states(i,k,j)%prob_diffusion(1,0,0,i_class), 1.0d0)
                 call aero_state_sample_particles_by_class(aero_states(i,k,j), &
                      delta_aero_states(i+1,k,j), aero_data, prob, &
                      AERO_INFO_NONE, env_states(i,k,j)%rrho, &
                      env_states(i+1,k,j)%rrho, source_vol, dest_vol, i_class)
             end do
          end do
          end do
       end if
       if (pmc_ie == global_nx) then
          i = global_nx
          do j = max(pmc_js,2), min(pmc_je,global_ny-1)
             do k = pmc_ks, pmc_ke
             dest_vol = env_states(i-1,k,j)%cell_volume
             source_vol = env_states(i,k,j)%cell_volume
             do i_class = 1,aero_weight_array_n_class(aero_states(i,k,j)%awa)
                prob = min(env_states(i,k,j)%prob_advection(-1,0,0,i_class) + &
                     env_states(i,k,j)%prob_diffusion(-1,0,0,i_class), 1.0d0)
                call aero_state_sample_particles_by_class(aero_states(i,k,j), &
                     delta_aero_states(i-1,k,j), aero_data, prob, &
                     AERO_INFO_NONE, env_states(i,k,j)%rrho, &
                     env_states(i-1,k,j)%rrho, &
                     source_vol, dest_vol, i_class)
             end do
             end do
          end do
       end if
    end if
    if (.not. config_flags%periodic_y) then
       ! South boundary
       if (pmc_js == 1) then
          j = 1
          do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
             do k = pmc_ks, pmc_ke
             dest_vol = env_states(i,k,j+1)%cell_volume
             source_vol = env_states(i,k,j)%cell_volume
             do i_class = 1,aero_weight_array_n_class(aero_states(i,k,j)%awa)
                prob = min(env_states(i,k,j)%prob_advection(0,0,1,i_class) + &
                     env_states(i,k,j)%prob_diffusion(0,0,1,i_class), 1.0d0)
                call aero_state_sample_particles_by_class(aero_states(i,k,j), &
                     delta_aero_states(i,k,j+1), aero_data, prob, &
                     AERO_INFO_NONE, env_states(i,k,j)%rrho, env_states(i,k,j+1)%rrho, &
                     source_vol, dest_vol, i_class)
                end do
             end do
          end do
       end if
       ! North boundary
       if (pmc_je == global_ny) then
          j = global_ny
          do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
             do k = pmc_ks, pmc_ke
             dest_vol = env_states(i,k,j-1)%cell_volume
             source_vol = env_states(i,k,j)%cell_volume
             do i_class = 1,aero_weight_array_n_class(aero_states(i,k,j)%awa)
                prob = min(env_states(i,k,j)%prob_advection(0,0,-1,i_class) + &
                     env_states(i,k,j)%prob_diffusion(0,0,-1,i_class), 1.0d0)
                call aero_state_sample_particles_by_class(aero_states(i,k,j), &
                     delta_aero_states(i,k,j-1), aero_data, prob, &
                     AERO_INFO_NONE, env_states(i,k,j)%rrho, env_states(i,k,j-1)%rrho, &
                     source_vol, dest_vol, i_class)
             end do
             end do
          end do
       end if
    end if
#ifdef PMC_DEBUG
    t2 = MPI_Wtime()
    print*, 'end boundary conditions. time spent = ', t2-t1
    print*, 'adding particles to aero_states'
#endif

!    print*, 'investigate the delta aero_states'
!    call aero_states_debug(delta_aero_states, pmc_is, pmc_ie, pmc_ks, pmc_ke, &
!         pmc_js, pmc_je)

    do j = j_start, j_end
    do k = pmc_ks,pmc_ke
    do i = i_start, i_end
       n_parts_before_add = aero_state_total_particles(aero_states(i,k,j))
       call aero_state_add_particles(aero_states(i,k,j), &
            delta_aero_states(i,k,j), aero_data)
       n_parts_after_add = aero_state_total_particles(aero_states(i,k,j))
       ! Warning if too many particles are being added
#ifdef PMC_DEBUG
       if ((n_parts_after_add - n_parts_before_add) > 2.0 * &
            real(grid%num_particles,kind=dp)) then
          write(*,*) "Adding too many particles to ",i,j,k, &
                n_parts_after_add - n_parts_before_add, 'particles were added'
          ! Debug
          call aero_state_debug(aero_states(i,k,j), max_group, max_class)
          call aero_state_debug(delta_aero_states(i,k,j), max_group, max_class)
          write(*,*) 'expected', expected_num_conc(i,k,j,max_group,max_class), &
            'actual', aero_state_group_class_num_conc(aero_states(i,k,j), &
             aero_data, max_group, max_class)
       end if
#endif
       call aero_state_rebalance(aero_states(i,k,j), aero_data, &
            grid%allow_doubling, grid%allow_halving, &
            initial_state_warning=.false.)
       call aero_state_sort(aero_states(i,k,j), aero_data)
    end do
    end do
    end do

!    print*, 'done adding particles and rebalanced'

!    call aero_states_debug(aero_states, pmc_is, pmc_ie, pmc_ks, pmc_ke, &
!         pmc_js, pmc_je)

    if (.not. config_flags%periodic_x) then
       if (pmc_is == 1) then
          i = 1
          j = pmc_js
          k = pmc_ks
          change_boundary_conditions = .false.
          p_tm1 = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time - grid%dt)
          p_t = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time)
          if (p_t > p_tm1 .and. p_t > 1) then
             change_boundary_conditions = .true.
          end if

          do j = max(pmc_js,2), min(pmc_je,global_ny-1)
             do k = pmc_ks, pmc_ke
                if (grid%u_2(i,k,j) < 0.0d0) then
                   call aero_state_dup_all_particles(aero_states(i+1,k,j), &
                        aero_states(i,k,j))
                else
                   throw_away = .false.
                   call aero_dist_interp_1d(scenario(i,k,j)%aero_background, &
                        scenario(i,k,j)%aero_dilution_time, &
                        scenario(i,k,j)%aero_dilution_rate, &
                        env_states(i,k,j)%elapsed_time, background, &
                        dilution_rate)
                   if (grid%u_2(i,k,j) * grid%u_1(i,k,j) < 0.0d0) throw_away = .true.
                   if (change_boundary_conditions) throw_away = .true. 
                   call aero_state_resample(aero_states(i,k,j), aero_data, &
                        background, (1.0d0 / env_states(i,k,j)%rrho), 1.0d0, &
                        0.0d0, .true., .true., throw_away, n_part_add)
                end if
             end do
          end do
       end if

       if (pmc_ie == global_nx) then
          i = global_nx
          j = pmc_js
          k = pmc_ks
          change_boundary_conditions = .false.
          p_tm1 = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time - grid%dt)
          p_t = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time)
          if (p_t > p_tm1 .and. p_t > 1) then
             change_boundary_conditions = .true.
          end if

          do j = max(pmc_js,2), min(pmc_je,global_ny-1)
             do k = pmc_ks, pmc_ke
                if (grid%u_2(i+1,k,j) > 0.0d0) then
                   call aero_state_dup_all_particles(aero_states(i-1,k,j), &
                        aero_states(i,k,j))
                else
                  throw_away = .false.
                   call aero_dist_interp_1d(scenario(i,k,j)%aero_background, &
                        scenario(i,k,j)%aero_dilution_time, scenario(i,k,j)%aero_dilution_rate, &
                        env_states(i,k,j)%elapsed_time, background, dilution_rate)
                   if (grid%u_2(i+1,k,j) * grid%u_1(i+1,k,j) < 0.0d0) throw_away = .true.
                   if (change_boundary_conditions) throw_away = .true.
                   call aero_state_resample(aero_states(i,k,j), aero_data, &
                        background, (1.0d0 / env_states(i,k,j)%rrho), 1.0d0, &
                        0.0d0, .true., .true., throw_away, n_part_add) 
                end if 
             end do
          end do
       end if
    end if     

    if (.not. config_flags%periodic_y) then
       if (pmc_js == 1) then
          j = 1
          i = pmc_is
          k = pmc_ks
          change_boundary_conditions = .false.
          p_tm1 = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time - grid%dt)
          p_t = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time)
          if (p_t > p_tm1 .and. p_t > 1) then
             change_boundary_conditions = .true.
          end if

          do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
             do k = pmc_ks, pmc_ke
                if (grid%v_2(i,k,j) < 0.0d0) then
                   call aero_state_dup_all_particles(aero_states(i,k,j+1), &
                        aero_states(i,k,j))
                else
                   throw_away = .false.
                   call aero_dist_interp_1d(scenario(i,k,j)%aero_background, &
                        scenario(i,k,j)%aero_dilution_time, scenario(i,k,j)%aero_dilution_rate, &
                        env_states(i,k,j)%elapsed_time, background, dilution_rate)
                   if (grid%v_2(i,k,j) * grid%v_1(i,k,j) < 0.0d0) throw_away = .true.
                   if (change_boundary_conditions) throw_away = .true.
                   call aero_state_resample(aero_states(i,k,j), aero_data, &
                        background, (1.0d0 / env_states(i,k,j)%rrho), 1.0d0, &
                        0.0d0, .true., .true., throw_away, n_part_add) 
                end if
             end do
          end do
       end if

       if (pmc_je == global_ny) then
          j = global_ny
          i = pmc_is
          k = pmc_ks
          change_boundary_conditions = .false.
          p_tm1 = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time - grid%dt)
          p_t = find_1d(size(scenario(i,k,j)%aero_background), &
               scenario(i,k,j)%aero_dilution_time, &
               env_states(i,k,j)%elapsed_time)
          if (p_t > p_tm1 .and. p_t > 1) then
             change_boundary_conditions = .true.
          end if

          do i = max(pmc_is,2), min(pmc_ie,global_nx-1)
             do k = pmc_ks, pmc_ke
                if (grid%v_2(i,k,j+1) > 0.0d0) then
                   call aero_state_dup_all_particles(aero_states(i,k,j-1), &
                        aero_states(i,k,j))
                else
                   throw_away = .false.
                   call aero_dist_interp_1d(scenario(i,k,j)%aero_background, &
                        scenario(i,k,j)%aero_dilution_time, scenario(i,k,j)%aero_dilution_rate, &
                        env_states(i,k,j)%elapsed_time, background, dilution_rate)
                   if (grid%v_2(i,k,j+1) * grid%v_1(i,k,j+1) < 0.0d0) throw_away = .true.
                   if (change_boundary_conditions) throw_away = .true.
                   call aero_state_resample(aero_states(i,k,j), aero_data, &
                        background, (1.0d0 / env_states(i,k,j)%rrho), 1.0d0, &
                        0.0d0, .true., .true., throw_away, n_part_add)
                end if
             end do
          end do
       end if
    end if

    deallocate(delta_aero_states)

#ifdef PMC_DEBUG
    total_particles = 0
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       total_particles = total_particles + aero_state_total_particles(aero_states(i,k,j))
    end do
    end do
    end do
    print "(A50,I10)", 'Total particles after transport = ', total_particles
#endif

  end subroutine trans_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resamples a distribution to replenish removed particles.
  subroutine aero_state_resample(aero_state, aero_data, &
       aero_dist, density, sample_prop, create_time, allow_doubling, allow_halving, &
       throw_away_previous, n_part_add)

    !> Aero state to add to.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Distribution to sample.
    type(aero_dist_t), intent(in) :: aero_dist
    !>
    real(kind=dp), intent(in) :: density 
    !> Fraction to sample (1).
    real(kind=dp), intent(in) :: sample_prop
    !> Creation time for new particles (s).
    real(kind=dp), intent(in) :: create_time
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving
    !>
    logical, intent(in) :: throw_away_previous
    !> Number of particles added.
    integer, intent(out), optional :: n_part_add

    real(kind=dp) :: n_samp_avg, radius, total_vol
    real(kind=dp) :: vols(aero_data_n_spec(aero_data))
    integer :: n_samp, i_mode, i_samp, i_group, i_class, n_group, n_class
    type(aero_particle_t) :: aero_particle

    ! Local variables
    integer, allocatable :: counter(:)
    integer :: i_source(1)
    integer :: n_diff, n_current, n_part
    integer :: n_part_before
    integer :: i_part, i_entry
    integer :: n_spec, n_source, n_samp_mode
    real(kind=dp) :: mix_rat_before
    real(kind=dp) :: total_num_conc, sample_proportion, n_parts_new
    integer :: n_parts_to_add
    real(kind=dp), parameter :: sample_timescale = 1.0d0
    real(kind=dp) :: size_factor
    real(kind=dp) :: num_conc_threshold
    real(kind=dp), parameter :: factor = .1d0

    n_group = size(aero_state%awa%weight, 1)
    n_class = size(aero_state%awa%weight, 2)

    if (present(n_part_add)) then
       n_part_add = 0
    end if

    if (throw_away_previous) then
       call aero_state_zero(aero_state)
       call aero_state_add_aero_dist_sample(aero_state, aero_data, &
            aero_dist, density, sample_timescale, create_time, allow_doubling, allow_halving, &
            n_part_add)
    else
    do i_group = 1,n_group
    do i_class = 1,n_class
       n_part_before = aero_state_total_particles(aero_state, i_group, &
            i_class)
       total_num_conc = 0.0d0
       do i_mode = 1, aero_dist_n_mode(aero_dist)
          if(aero_dist%mode(i_mode)%weight_class == i_class) then
             total_num_conc = total_num_conc + &
                  aero_mode_total_num_conc(aero_dist%mode(i_mode)) * density
          end if
       end do
       n_parts_new = total_num_conc &
            / aero_state%awa%weight(i_group, i_class)%magnitude
       n_samp = rand_poisson(n_parts_new)
       n_parts_to_add = n_samp - n_part_before 
       if (n_parts_to_add > 0) then
          do i_mode = 1,aero_dist_n_mode(aero_dist)
             if(aero_dist%mode(i_mode)%weight_class == i_class) then
             sample_proportion = (density * (aero_dist%mode(i_mode)%num_conc)) &
                  / total_num_conc
             n_samp_mode = rand_poisson(sample_proportion * n_parts_to_add)
             n_samp_avg = sample_proportion * n_parts_to_add
             num_conc_threshold = (density * aero_mode_total_num_conc(aero_dist%mode(i_mode)) &
                  / n_samp_avg) * factor
             size_factor = min((1.0d0 / (sample_timescale*n_samp_avg)), 0.32d0)
             if (n_samp_mode > 0) then
             do i_samp = 1,n_samp_mode
                call aero_particle_zero(aero_particle, aero_data)
                call aero_mode_sample_radius(aero_dist%mode(i_mode), &
                     aero_data, aero_state%awa%weight(i_group, i_class), &
                     radius, size_factor)
                total_vol = aero_data_rad2vol(aero_data, radius)
                call aero_mode_sample_vols(aero_dist%mode(i_mode), total_vol, &
                     vols)
                call aero_particle_set_vols(aero_particle, vols)
                call aero_particle_new_id(aero_particle)
                call aero_particle_set_weight(aero_particle, i_group, i_class)
                call aero_particle_set_component(aero_particle, &
                     aero_dist%mode(i_mode)%source, create_time)
                if (aero_weight_array_num_conc(aero_state%awa, aero_particle, &
                     aero_data) > num_conc_threshold) then
                   call aero_state_add_particle(aero_state, aero_particle, &
                        aero_data)
                end if
             end do
             end if
             end if
          end do
       else
          n_part = integer_varray_n_entry( &
              aero_state%aero_sorted%group_class%inverse(i_group, i_class))
          do i_samp = 1,min(-n_parts_to_add,n_part)
             i_entry = pmc_rand_int(integer_varray_n_entry( &
                  aero_state%aero_sorted%group_class%inverse(i_group, &
                  i_class)))
             i_part = aero_state%aero_sorted%group_class%inverse(i_group, &
                  i_class)%entry(i_entry)
             call aero_state_remove_particle_no_info(aero_state, i_part)
          end do
       end if
       end do
       end do
       end if

  end subroutine aero_state_resample

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adjust the weight magnitudes before sampling based on expected values.
  subroutine trans_aero_preweight(grid, config_flags, aero_data, scenario, &
       env_states, aero_states, is, ie, js, je, ks, ke, &
       pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, &
       global_nx, global_ny, global_nz, expected_num_conc)

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
    ! FIXME:
    integer, intent(in) :: is,ie,js,je,ks,ke
    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> WRF namelist configuration.
    type(grid_config_rec_type), intent(in) :: config_flags
    !> Full domain of aerosol data.
    type(aero_data_t), intent(in) :: aero_data 
    !> Scenario data.
    type(scenario_t),dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenario
    !> Full domain of environment states.
    type(env_state_t), dimension(is:ie,ks:ke,js:je), &
         intent(in) :: env_states
    !> Full domain of aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: aero_states
!    integer, intent(in) :: is,ie,js,je,ks,ke
!    !> Starting PartMC index in east-west direction.
!    integer, intent(in) :: pmc_is
!    !> Ending PartMC index in east-west direction.
!    integer, intent(in) :: pmc_ie
!    !> Starting PartMC index in vertical direction.
!    integer, intent(in) :: pmc_ks
!    !> Ending PartMC index in vertical direction.
!    integer, intent(in) :: pmc_ke
!    !> Starting PartMC index in north-south direction.
!    integer, intent(in) :: pmc_js
!    !> Ending PartMC index in north-south direction.
!    integer, intent(in) :: pmc_je
    !> Number of mass points in east-west direction.
    integer, intent(in) :: global_nx
    !> Number of mass points in the north-south direction.
    integer, intent(in) :: global_ny
    !> Number of mass points in the vertical direction.
    integer, intent(in) :: global_nz
    real(kind=dp), allocatable, intent(inout), dimension(:,:,:,:,:) :: expected_num_conc

    ! Local variables
    integer :: i,j,k
    integer :: i_group, i_class
    integer :: n_group, n_class
    integer :: i_part
    real(kind=dp) :: weight_ratio, num_conc, val
    real(kind=dp), allocatable, dimension(:,:,:,:,:) :: &
         num_conc_old, num_conc_new
    real(kind=dp), allocatable, dimension(:,:,:,:,:) :: temp_5d_array, &
         send_array
    logical :: ideal
    character, dimension(:), allocatable :: send_buffer, recv_buffer
    integer :: n_parts_group_class
    real(kind=dp) :: n_parts_new
    integer :: dest, source, ierr, disp, dir
    integer :: i_send_s, i_send_e, i_recv_s, i_recv_e
    integer :: j_send_s, j_send_e, j_recv_s, j_recv_e
    integer :: k_send_s, k_send_e, k_recv_s, k_recv_e
    integer :: i_inc, j_inc, k_inc
    integer :: tag, position, max_size_send, max_size_recv, request
    integer :: send_request_size, send_request
    integer :: status(MPI_STATUS_SIZE)
    integer :: i_entry, n_part

    ideal = .false.
    if (config_flags%periodic_x .and. config_flags%periodic_y) ideal = .true. 

    ! For each grid cell
    !  - Compute your number concentrations per source.
    !  - Compute the fluxes using losses
    !  - Send your fluxes per source to neighbors.
    !  - Compute future number concentration.
    !     - If this increases a lot, we want to adjust the weighting
    !       magnitude such that we don't have extreme values of gamma.

    n_group = aero_weight_array_n_group(aero_states(pmc_is,pmc_ks,pmc_js)%awa)
    n_class = aero_weight_array_n_class(aero_states(pmc_is,pmc_ks,pmc_js)%awa)

    allocate(num_conc_old(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je, n_group, &
         n_class))

!    call get_array_extents(pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
!         global_nx, global_ny, global_nz, is, ie, js, je, ks, ke, config_flags)

    allocate(num_conc_new(is:ie,ks:ke,js:je,n_group,n_class))
    allocate(expected_num_conc(is:ie,ks:ke,js:je,n_group,n_class))
    num_conc_old = 0.0d0
    num_conc_new = 0.0d0

    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
       call aero_state_sort(aero_states(i,k,j), aero_data)
    do i_group = 1,n_group
    do i_class = 1,n_class
       n_parts_group_class = aero_state_total_particles(aero_states(i,k,j), &
            i_group, i_class)
       num_conc_old(i,k,j,i_group,i_class) = n_parts_group_class * &
            aero_states(i,k,j)%awa%weight(i_group, i_class)%magnitude !&
!            * env_states(i,k,j)%rrho
       num_conc_new(i,k,j,i_group,i_class) = num_conc_old(i,k,j,i_group,i_class)
    end do
    end do    
    end do
    end do
    end do

#ifdef PMC_DEBUG
    print*, 'Total of the array before', sum(num_conc_old)
#endif

    ! Fluxes per source for grid cells
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
    do i_group = 1,n_group
    do i_class = 1,n_class
       call compute_fluxes(grid, num_conc_old(i,k,j,i_group,i_class), &
            env_states(i,k,j), num_conc_new(:,:,:,i_group,i_class), &
            pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, is, ie, ks, ke, &
            js, je, global_nx, global_ny, global_nz, env_states, ideal, &
            i_group, i_class)
    end do
    end do
    end do
    end do
    end do
#ifdef PMC_DEBUG
    print*, "Total of the array after:", sum(num_conc_new)
#endif
#ifdef PMC_DEBUG
    write(*,*) 'extents', pmc_mpi_rank(),is,ie,js,je,pmc_is,pmc_ie,pmc_js, pmc_je
#endif
    ! Num_conc_new is updated with all losses and local gains
    ! Send/receive any additional fluxes to add



    tag = 0
    do dir = 1,2
    do disp = -1,1,2

       if (ideal) then
          call MPI_cart_shift(local_communicator_periodic, (dir - 1), disp, &
               source, dest, ierr)
       else
          call MPI_cart_shift(local_communicator, (dir - 1), disp, &
               source, dest, ierr)
       end if

       ! Get the index information for arrays
       call get_array_index_for_send_recv_from_ghost(disp, dir, is, ie, js, je, &
            ks, ke, i_send_s, i_send_e, j_send_s, j_send_e, k_send_s, &
            k_send_e, i_recv_s, i_recv_e, j_recv_s, j_recv_e, k_recv_s, &
            k_recv_e)

       if (source /= MPI_PROC_NULL) then
          if (ideal) then
             call MPI_Irecv(max_size_recv, 1, MPI_INTEGER, source, tag, &
                  local_communicator_periodic, request, ierr)
          else
             call MPI_Irecv(max_size_recv, 1, MPI_INTEGER, source, tag, &
                  local_communicator, request, ierr)
          end if
       end if
#ifdef PMC_DEBUG
       print*, 'array sending away', &
            sum(num_conc_new(i_send_s:i_send_e,k_send_s:k_send_e, &
            j_send_s:j_send_e,:,:))
#endif
       if (dest /= MPI_PROC_NULL) then
          send_array = num_conc_new(i_send_s:i_send_e,k_send_s:k_send_e, &
               j_send_s:j_send_e,:,:)
          max_size_send = pmc_mpi_pack_size_real_array_5d(send_array)
          if (ideal) then
             call MPI_isend(max_size_send, 1, MPI_INTEGER, dest, tag, &
                  local_communicator_periodic, send_request_size, ierr)
          else
             call MPI_isend(max_size_send, 1, MPI_INTEGER, dest, tag, &
                  local_communicator, send_request_size, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          call MPI_wait(send_request_size, status, ierr)
       end if

       ! Send and receive
       if (dest /= MPI_PROC_NULL) then
          ! Pack the buffer
          position = 0
          allocate(send_buffer(max_size_send))
          call pmc_mpi_pack_real_array_5d(send_buffer, position, send_array)

          if (ideal) then
             call MPI_isend(send_buffer, max_size_send, MPI_PACKED, dest, tag, &
                  local_communicator_periodic, send_request, ierr)
          else
             call MPI_isend(send_buffer, max_size_send, MPI_PACKED, dest, tag, &
                  local_communicator, send_request, ierr)
          end if
       end if

       if (source /= MPI_PROC_NULL) then
          call MPI_wait(request, status, ierr)
          ! Non-blocking receiving of buffer from source
          allocate(recv_buffer(max_size_recv))
          if (ideal) then
             call MPI_Irecv(recv_buffer, max_size_recv, MPI_PACKED, source, tag, &
                  local_communicator_periodic, request, ierr)
          else
             call MPI_Irecv(recv_buffer, max_size_recv, MPI_PACKED, source, tag, &
                  local_communicator, request, ierr)
          end if
       end if

       if (dest /= MPI_PROC_NULL) then
          call MPI_wait(send_request, status, ierr)
       end if

       if (source /= MPI_PROC_NULL) then
          ! Wait
          call MPI_wait(request, status, ierr)
          ! Unpack
          position = 0
          call pmc_mpi_unpack_real_array_5d(recv_buffer, position, &
               temp_5d_array)
          do j_inc = j_recv_s, j_recv_e
          do k_inc = k_recv_s, k_recv_e
          do i_inc = i_recv_s, i_recv_e
             ! temp_5d_array is indexed starting at 1
             num_conc_new(i_inc,k_inc,j_inc,:,:) = &
                num_conc_new(i_inc,k_inc,j_inc,:,:) + &
                temp_5d_array(i_inc - i_recv_s + 1, k_inc - k_recv_s + 1, &
                     j_inc - j_recv_s + 1,:,:)
          end do
          end do
          end do
       end if

       if (allocated(send_buffer)) deallocate(send_buffer)
       if (allocated(recv_buffer)) deallocate(recv_buffer)
       if (allocated(temp_5d_array)) deallocate(temp_5d_array)
    end do
    end do

#ifdef PMC_DEBUG
    print*, 'array after all movements (real points only):', &
       sum(num_conc_new(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je, &
          :, :))
#endif

    ! Compare new and old number concentrations per class and group
!    if (pmc_mpi_rank() == 981) then
!    print*, 'starting the rebalance'
!    end if
    do j = pmc_js,pmc_je
    do k = pmc_ks,pmc_ke
    do i = pmc_is,pmc_ie
    do i_group = 1,n_group
    do i_class = 1,n_class
      n_parts_new = num_conc_new(i,k,j,i_group,i_class) / &
           aero_states(i,k,j)%awa%weight(i_group, i_class)%magnitude !&
!      if (pmc_mpi_rank() == 981) then
!         print*, i,k,j,i_class,i_group, num_conc_new(i,k,j,i_group,i_class), &
!              num_conc_old(i,k,j,i_group,i_class), n_parts_new, aero_states(i,k,j)%awa%weight(i_group, i_class)%magnitude
!      end if
!           / env_states(i,k,j)%rrho
      if (n_parts_new > 2.0d0 * &
         aero_states(i,k,j)%n_part_ideal(i_group, i_class)) then
         weight_ratio = n_parts_new / &
              aero_states(i,k,j)%n_part_ideal(i_group, i_class)
         call aero_state_scale_weight(aero_states(i,k,j), &
              aero_data, i_group, &
              i_class, weight_ratio, .true., .true.)
!         if (pmc_mpi_rank() == 981) then
!         print*, 'adjusting', i,k,j, i_group, i_class, weight_ratio, &
!              aero_state_total_particles(aero_states(i,k,j), i_group, i_class)
!         end if
      end if
    end do
    end do
    end do
    end do
    end do

!    if (pmc_mpi_rank() == 981) then
!    print*, 'looping over cells after balancing'
!    do j = pmc_js,pmc_je
!    do k = pmc_ks,pmc_ke
!    do i = pmc_is,pmc_ie
!    do i_group = 1,n_group
!    do i_class = 1,n_class
!        print*, i,k,j, aero_states(i,k,j)%awa%weight(i_group, i_class)%magnitude, &
!            aero_state_total_particles(aero_states(i,k,j), i_group, i_class), &
!            num_conc_new(i,k,j,i_group,i_class)
!    end do
!    end do
!    end do
!    end do
!    end do
!    end if

    expected_num_conc = num_conc_new

  end subroutine trans_aero_preweight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the fluxes of a grid cell to determine new value.
  subroutine compute_fluxes(grid, num_conc_old, env_state, num_conc_new, &
      pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, is, ie, ks, ke, js, je, &
      global_nx, global_ny, global_nz, env_states, ideal, i_group, i_class)

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
    !> FIXME:
    integer, intent(in) :: is,ie,js,je,ks,ke
    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> Number concentration at time t.
    real(kind=dp), intent(in) :: num_conc_old
    !> Environmental state.
    type(env_state_t), intent(in) :: env_state
    !> Number concentrations at time t+dt.
    real(kind=dp), intent(inout), dimension(is:ie,ks:ke,js:je) :: num_conc_new
    !> 
    integer, intent(in) :: global_nx
    !>
    integer, intent(in) :: global_ny
    !>
    integer, intent(in) :: global_nz
    !> Environmental states
    type(env_state_t), dimension(is:ie,ks:ke,js:je), &
         intent(in) :: env_states
    !> Logical flag for periodic boundary conditions
    logical, intent(in) :: ideal
    !>
    integer, intent(in) :: i_group
    !>
    integer, intent(in) :: i_class


    integer :: i_index, j_index, k_index
    integer :: i,j,k, i_part
    real(kind=dp) :: flux, prob, dest_vol, source_vol

    i_index = env_state%ix
    j_index = env_state%iy
    k_index = env_state%iz

    source_vol = env_states(i_index,k_index,j_index)%cell_volume
    if (ideal) then
       k = 0
       i = 0
       do j = -1,1,2
          dest_vol = env_states(i_index,k_index,j_index+j)%cell_volume
          prob = env_state%prob_advection(i,k,j,i_class) + &
               env_state%prob_diffusion(i,k,j,i_class)
          flux = prob * num_conc_old
          num_conc_new(i_index,k_index,j_index+j) = &
               num_conc_new(i_index,k_index,j_index+j) + flux &
               * (dest_vol / source_vol)
          num_conc_new(i_index,k_index,j_index) =  &
               num_conc_new(i_index,k_index,j_index) - flux
       end do

       j = 0
       do i = -1,1,2
          dest_vol = env_states(i_index+i,k_index,j_index)%cell_volume
          prob = env_state%prob_advection(i,k,j,i_class) + &
               env_state%prob_diffusion(i,k,j,i_class)
          flux = prob * num_conc_old
          num_conc_new(i_index+i,k_index,j_index) = &
              num_conc_new(i_index+i,k_index,j_index) + flux &
              * (dest_vol / source_vol)
          num_conc_new(i_index,k_index,j_index) = &
               num_conc_new(i_index,k_index,j_index) - flux
       end do

       do k = pmc_ks,pmc_ke
          if (k .ne. k_index) then
             prob = env_state%prob_vert_diffusion(k,i_class)
             dest_vol = env_states(i_index,k,j_index)%cell_volume
             flux = prob * num_conc_old
             num_conc_new(i_index,k,j_index) = &
                  num_conc_new(i_index,k,j_index) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) = &
                  num_conc_new(i_index,k_index,j_index) - flux
          end if
       end do
    else ! Real case
       if (i_index == 1) then
          if (j_index .ne. 1 .and. j_index .ne. global_ny) then
             k = 0
             j = 0
             i = 1
             dest_vol = env_states(i_index+i,k_index,j_index)%cell_volume
             prob = env_state%prob_advection(i,k,j,i_class) + &
                  env_state%prob_diffusion(i,k,j,i_class)
             flux = prob * num_conc_old
             num_conc_new(i_index+i,k_index,j_index) = &
                  num_conc_new(i_index+i,k_index,j_index) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) = &
                  num_conc_new(i_index,k_index,j_index) - flux
          end if
       else if (i_index == global_nx) then
          if (j_index .ne. 1 .and. j_index .ne. global_ny) then
             k = 0
             j = 0
             i = -1
             dest_vol = env_states(i_index+i,k_index,j_index)%cell_volume
             prob = env_state%prob_advection(i,k,j,i_class) + &
                  env_state%prob_diffusion(i,k,j,i_class)
             flux = prob * num_conc_old
             num_conc_new(i_index+i,k_index,j_index) = &
                  num_conc_new(i_index+i,k_index,j_index) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) = &
                 num_conc_new(i_index,k_index,j_index) - flux
          end if
       else if (j_index == 1) then
          if (i_index .ne. 1 .and. i_index .ne. global_nx) then
             k = 0
             i = 0
             j = 1
             dest_vol = env_states(i_index,k_index,j_index+j)%cell_volume
             prob = env_state%prob_advection(i,k,j,i_class) + &
                  env_state%prob_diffusion(i,k,j,i_class)
             flux = prob * num_conc_old
             num_conc_new(i_index,k_index,j_index+j) = &
                  num_conc_new(i_index,k_index,j_index+j) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) =  &
                  num_conc_new(i_index,k_index,j_index) - flux
          end if
       else if (j_index == global_ny) then
          if (i_index .ne. 1 .and. i_index .ne. global_nx) then
             k = 0
             i = 0
             j = -1
             dest_vol = env_states(i_index,k_index,j_index+j)%cell_volume
             prob = env_state%prob_advection(i,k,j,i_class) + &
                  env_state%prob_diffusion(i,k,j,i_class)
             flux = prob * num_conc_old
             num_conc_new(i_index,k_index,j_index+j) = &
                  num_conc_new(i_index,k_index,j_index+j) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) =  &
                  num_conc_new(i_index,k_index,j_index) - flux
          end if
       else
          k = 0
          i = 0
          do j = -1,1,2
             dest_vol = env_states(i_index,k_index,j_index+j)%cell_volume
             prob = env_state%prob_advection(i,k,j,i_class) + &
                  env_state%prob_diffusion(i,k,j,i_class)
             flux = prob * num_conc_old
             num_conc_new(i_index,k_index,j_index+j) = &
                  num_conc_new(i_index,k_index,j_index+j) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) =  &
                  num_conc_new(i_index,k_index,j_index) - flux
          end do

          j = 0
          do i = -1,1,2
             dest_vol = env_states(i_index+i,k_index,j_index)%cell_volume
             prob = env_state%prob_advection(i,k,j,i_class) + &
                  env_state%prob_diffusion(i,k,j,i_class)
             flux = prob * num_conc_old
             num_conc_new(i_index+i,k_index,j_index) = &
                  num_conc_new(i_index+i,k_index,j_index) + flux * &
                  (dest_vol / source_vol)
             num_conc_new(i_index,k_index,j_index) = &
                  num_conc_new(i_index,k_index,j_index) - flux
          end do

          do k = pmc_ks,pmc_ke
             if (k .ne. k_index) then
                prob = env_state%prob_vert_diffusion(k,i_class)
                dest_vol = env_states(i_index,k,j_index)%cell_volume
                flux = prob * num_conc_old
                num_conc_new(i_index,k,j_index) = &
                     num_conc_new(i_index,k,j_index) + flux  * &
                     (dest_vol / source_vol)
                num_conc_new(i_index,k_index,j_index) = &
                     num_conc_new(i_index,k_index,j_index) - flux
             end if
          end do
       end if
    end if

  end subroutine compute_fluxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the local domain extent.
  subroutine get_array_extents(pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, &
      pmc_ke, global_nx, global_ny, global_nz, is, ie, js, je, ks, ke, &
      config_flags)

    !>
    integer, intent(in) :: pmc_is
    !>
    integer, intent(in) :: pmc_ie
    !>
    integer, intent(in) :: pmc_js
    !>
    integer, intent(in) :: pmc_je
    !>
    integer, intent(in) :: pmc_ks
    !>
    integer, intent(in) :: pmc_ke
    !>
    integer, intent(in) :: global_nx
    !>
    integer, intent(in) :: global_ny
    !>
    integer, intent(in) :: global_nz
    !>
    integer, intent(out) :: is
    !>
    integer, intent(out) :: ie
    !>
    integer, intent(out) :: js
    !>
    integer, intent(out) :: je
    !>
    integer, intent(out) :: ks
    !>
    integer, intent(out) :: ke
    !> WRF namelist configuration.
    type(grid_config_rec_type), intent(in) :: config_flags

    if (config_flags%periodic_x .and. config_flags%periodic_y) then
       is = pmc_is - 1
       ie = pmc_ie + 1
       js = pmc_js - 1
       je = pmc_je + 1
    else
       if (pmc_is == 1) then
          if (pmc_ie /= global_nx) then
             is = pmc_is
             ie = pmc_ie + 1
          else
             is = pmc_is
             ie = pmc_ie
          end if
       ! East boundary
       else if (pmc_ie == global_nx) then
          is = pmc_is - 1
          ie = pmc_ie
       else
          is = pmc_is - 1
          ie = pmc_ie + 1
       end if
       ! South boundary
       if (pmc_js == 1) then
          if (pmc_je /= global_ny) then
             js = pmc_js
             je = pmc_je + 1
          else
             js = pmc_js
             je = pmc_je
          end if
       ! North boundary
       else if (pmc_je == global_ny) then
          js = pmc_js - 1
          je = pmc_je
       else
          js = pmc_js - 1
          je = pmc_je + 1
       end if
    end if

    ks = pmc_ks
    ke = pmc_ke

  end subroutine get_array_extents

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Array index for sending to neighbor ghosts and receiving in local ghosts.
  subroutine get_array_index_for_send_recv_to_ghost(disp, dir, is, ie, js, je, &
       ks, ke, i_send_s, i_send_e, j_send_s, j_send_e, k_send_s, &
       k_send_e, i_recv_s, i_recv_e, j_recv_s, j_recv_e, k_recv_s, &
       k_recv_e)

    !>
    integer, intent(in) :: disp
    !>
    integer, intent(in) :: dir
    !>
    integer, intent(in) :: is
    !>
    integer, intent(in) :: ie
    !>
    integer, intent(in) :: js
    !>
    integer, intent(in) :: je
    !>
    integer, intent(in) :: ks
    !>
    integer, intent(in) :: ke
    !>
    integer, intent(out) :: i_send_s
    !>
    integer, intent(out) :: i_send_e
    !>
    integer, intent(out) :: j_send_s
    !>
    integer, intent(out) :: j_send_e
    !>
    integer, intent(out) :: k_send_s
    !>
    integer, intent(out) :: k_send_e
    !>
    integer, intent(out) :: i_recv_s
    !>
    integer, intent(out) :: i_recv_e
    !>
    integer, intent(out) :: j_recv_s
    !>
    integer, intent(out) :: j_recv_e
    !>
    integer, intent(out) :: k_recv_s
    !>
    integer, intent(out) :: k_recv_e

    if (dir == 1) then ! Moving in the y direction
       if (disp == -1) then ! Moving in the negative direction
          ! what we send
          j_send_s = js + 1 !1
          j_send_e = js + 1 !1
          ! where it goes in the array
          j_recv_s = je !ny+1
          j_recv_e = je !ny+1
       else if (disp == 1) then ! Moving in the positive direction
          ! what we send
          j_send_s = je - 1 !ny
          j_send_e = je - 1 !ny
          ! where it goes in the array
          j_recv_s = js !0
          j_recv_e = js !0
       end if
       !
       i_send_s = is + 1 !1
       i_send_e = ie - 1 !nx
       i_recv_s = is + 1 !1
       i_recv_e = ie - 1 !nx
       !
       k_send_s = ks !1
       k_send_e = ke !nz
       k_recv_s = ks !1
       k_recv_e = ke !nz
    elseif (dir == 2) then ! Moving in the x direction
       if (disp == -1) then ! Moving in the negative direction
          ! what we send
          i_send_s = is + 1 !1
          i_send_e = is + 1 !1
          ! where it goes in the array
          i_recv_s = ie ! nx+1
          i_recv_e = ie ! nx+1
       else if (disp == 1) then ! Moving in the positive direction
          ! what we send
          i_send_s = ie - 1 !nx
          i_send_e = ie - 1 !nx
          ! where it goes in the array
          i_recv_s = is !0
          i_recv_e = is !0
       end if
       !
       j_send_s = js + 1 !1
       j_send_e = je - 1 !ny
       j_recv_s = js + 1 !1
       j_recv_e = je - 1 !ny
       !
       k_send_s = ks !1
       k_send_e = ke !nz
       k_recv_s = ks !1
       k_recv_e = ke !nz
    end if

#ifdef PMC_DEBUG
    print*,'to ghost', pmc_mpi_rank(),' sending', i_send_s, i_send_e, &
       j_send_s, j_send_e, 'recv ', i_recv_s, i_recv_e, j_recv_s, j_recv_e, &
       'dir', dir, 'disp', disp
#endif

  end subroutine get_array_index_for_send_recv_to_ghost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Array index to send values to neighbor and to receive.
  subroutine get_array_index_for_send_recv_from_ghost(disp, dir, is, ie, js, je, &
       ks, ke, i_send_s, i_send_e, j_send_s, j_send_e, k_send_s, &
       k_send_e, i_recv_s, i_recv_e, j_recv_s, j_recv_e, k_recv_s, &
       k_recv_e)

    !> Displacment.
    integer, intent(in) :: disp
    !> Direction.
    integer, intent(in) :: dir
    !> Start index of subdmain in west-east direction.
    integer, intent(in) :: is
    !> End index of subdomain in west-east direction.
    integer, intent(in) :: ie
    !> Start index of subdoamin in north-south direction.
    integer, intent(in) :: js
    !> End index of subdomain in north-south direction.
    integer, intent(in) :: je
    !> Start index in vertical direction.
    integer, intent(in) :: ks
    !> End index in vertical direction.
    integer, intent(in) :: ke
    !> Start index to send of west-east direction.
    integer, intent(out) :: i_send_s
    !> End index to send of west-east direction.
    integer, intent(out) :: i_send_e
    !> Start index to send of north-south direction.
    integer, intent(out) :: j_send_s
    !> End index to send of north-south direction.
    integer, intent(out) :: j_send_e
    !> Start index to send of vertical.
    integer, intent(out) :: k_send_s
    !> End index to send of vertical.
    integer, intent(out) :: k_send_e
    !> Start index of receive in west-east direction.
    integer, intent(out) :: i_recv_s
    !> End index of receive in west-east direction.
    integer, intent(out) :: i_recv_e
    !> Start index of receive in north-south direction.
    integer, intent(out) :: j_recv_s
    !> End index of receive in north-south direction.
    integer, intent(out) :: j_recv_e
    !> Start index of receive in vertical direction.
    integer, intent(out) :: k_recv_s
    !> End index of receive in vertical direction.
    integer, intent(out) :: k_recv_e

    if (dir == 1) then ! Moving in the y direction
       if (disp == -1) then ! Moving in the negative direction
          ! what we send
          j_send_s = js
          j_send_e = js
          ! where it goes in the array
          j_recv_s = je - 1
          j_recv_e = je - 1
       else if (disp == 1) then ! Moving in the positive direction
          ! what we send
          j_send_s = je
          j_send_e = je
          ! where it goes in the array
          j_recv_s = js + 1
          j_recv_e = js + 1
       end if
       !
       i_send_s = is + 1 !1
       i_send_e = ie - 1 !nx
       i_recv_s = is + 1 !1
       i_recv_e = ie - 1 !nx
       !
       k_send_s = ks !1
       k_send_e = ke !nz
       k_recv_s = ks !1
       k_recv_e = ke !nz

    elseif (dir == 2) then ! Moving in the x direction
       if (disp == -1) then ! Moving in the negative direction
          ! what we send
          i_send_s = is
          i_send_e = is
          ! where it goes in the array
          i_recv_s = ie - 1
          i_recv_e = ie - 1
       else if (disp == 1) then ! Moving in the positive direction
          ! what we send
          i_send_s = ie
          i_send_e = ie
          ! where it goes in the array
          i_recv_s = is + 1
          i_recv_e = is + 1
       end if
       !
       j_send_s = js + 1 !1
       j_send_e = je - 1 !ny
       j_recv_s = js + 1 !1
       j_recv_e = je - 1 !ny
       !
       k_send_s = ks !1
       k_send_e = ke !nz
       k_recv_s = ks !1
       k_recv_e = ke !nz
    end if

#ifdef PMC_DEBUG
    print*,'from ghost', pmc_mpi_rank(),' sending', i_send_s, i_send_e, &
       j_send_s, j_send_e, 'recv ', i_recv_s, i_recv_e, j_recv_s, j_recv_e, &
       'dir', dir, 'disp', disp
#endif

  end subroutine get_array_index_for_send_recv_from_ghost

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine aero_state_debug(aero_state, max_group, max_class)

    type(aero_state_t), intent(in) :: aero_state

    integer :: i_group, i_class
    integer :: n_group, n_class
    integer :: max_val, max_group, max_class
    integer :: n_part

    n_group = aero_weight_array_n_group(aero_state%awa)
    n_class = aero_weight_array_n_class(aero_state%awa)
    max_val = 0
    max_group = 1 
    max_class = 1 
    do i_class = 1, n_class
    do i_group = 1, n_group
      n_part = aero_state_total_particles(aero_state, i_group, i_class)
      if (n_part > max_val) then
         max_val = n_part
         max_group = i_group
         max_class = i_class
      end if
    end do
    end do

    print*, 'biggest group', max_val, max_group, max_class, &
       aero_state%awa%weight(max_group,max_class)%magnitude

  end subroutine aero_state_debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine aero_states_debug(aero_states, pmc_is, pmc_ie, pmc_ks, pmc_ke, &
       pmc_js, pmc_je)

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
    !> Full domain of aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: aero_states

    integer :: i_group, i_class
    integer :: n_group, n_class
    integer :: max_val, max_group, max_class
    integer :: n_part
    integer :: grid_cell_max(3)
    real(kind=dp) :: n_part_max, n_part_total
    integer :: i,k,j 

    n_part_max = 0.0d0 
    n_part_total = 0.0d0
    do i = pmc_is,pmc_ie
    do k = pmc_ks,pmc_ke
    do j = pmc_js,pmc_je
       n_part = aero_state_total_particles(aero_states(i,k,j))
       if (n_part > n_part_max) then
          n_part_max = n_part
          grid_cell_max(1) = i
          grid_cell_max(2) = k
          grid_cell_max(3) = j
       end if
       n_part_total = n_part_total + n_part
    end do
    end do
    end do

    print*, 'total particles', n_part_total, 'worst cell', n_part_max

    i = grid_cell_max(1)
    k = grid_cell_max(2)
    j = grid_cell_max(3)

    n_group = aero_weight_array_n_group(aero_states(i,k,j)%awa)
    n_class = aero_weight_array_n_class(aero_states(i,k,j)%awa)
    max_val = 0
    max_group = 1
    max_class = 1
    do i_class = 1, n_class
    do i_group = 1, n_group
      n_part = aero_state_total_particles(aero_states(i,k,j), i_group, i_class)
      if (n_part > max_val) then
         max_val = n_part
         max_group = i_group
         max_class = i_class
      end if
    end do
    end do

    print*, 'Largest group', max_val, max_group, max_class, &
       aero_states(i,k,j)%awa%weight(max_group,max_class)%magnitude

  end subroutine aero_states_debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine aero_state_sample_particles_vertical(aero_state_from, &
       aero_state_to, aero_data, sample_prob, removal_action, env_state_from, &
       env_state_to, i_class)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Probability of sampling each particle.
    real(kind=dp), intent(in) :: sample_prob
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    type(env_state_t), intent(in) :: env_state_from
    type(env_state_t), intent(in) :: env_state_to
    !>
    integer, intent(in) :: i_class

    integer :: n_transfer, i_transfer
    integer :: i_group
    logical :: do_add, do_remove
    real(kind=dp) :: num_conc_from, num_conc_to
    type(aero_info_t) :: aero_info
    integer :: n_parts_before
    type(aero_state_t) :: aero_state_temp
    real(kind=dp) :: factor, ratio
    real(kind=dp) :: rho_to, rho_from, mix_rat_to, mix_rat_from
    integer :: i_entry, i_part

    call assert(721006964, (sample_prob >= 0d0) .and. (sample_prob <= 1d0))
#ifndef PMC_USE_WRF
    call aero_state_zero(aero_state_to)
    call aero_state_copy_weight(aero_state_from, aero_state_to)
#else
    call aero_state_zero(aero_state_temp)
    call aero_state_set_weight(aero_state_temp, aero_data, &
         AERO_STATE_WEIGHT_FLAT_SPECIFIED)
    call aero_state_copy_weight(aero_state_to, aero_state_temp)
    factor = env_state_from%height / env_state_to%height
    call aero_weight_array_scale(aero_state_temp%awa, 1.0d0/factor)
#endif

#ifdef PMC_DEBUG
    n_parts_before = aero_state_total_particles(aero_state_to)
#endif

    ! Do we have a valid sort on the from state?
    if (.not. aero_state_from%valid_sort) then
       call aero_state_sort(aero_state_from, aero_data)
    end if

    rho_to = 1.0d0 !/ env_state_to%rrho
    rho_from = 1.0d0 !/ env_state_from%rrho 

    ! Loop over groups
    do i_group = 1,aero_weight_array_n_group(aero_state_from%awa)
       ratio = (aero_state_temp%awa%weight(i_group, i_class)%magnitude / rho_to) / &
            (aero_state_from%awa%weight(i_group, i_class)%magnitude  / rho_from)

       if (ratio >= 1.0d0) then
          n_transfer = rand_binomial(aero_state_total_particles(aero_state_from, &
               i_group, i_class), sample_prob)
       else
          n_transfer = rand_poisson(aero_state_total_particles( &
               aero_state_from, i_group, i_class) * (sample_prob / ratio))
       end if

       i_transfer = 0
       do while (i_transfer < n_transfer)
          if (aero_state_total_particles( &
               aero_state_from,i_group,i_class) <= 0) exit
          i_entry = pmc_rand_int(integer_varray_n_entry( &
               aero_state_from%aero_sorted%group_class%inverse(i_group, &
               i_class)))
          i_part = aero_state_from%aero_sorted%group_class%inverse(i_group, &
               i_class)%entry(i_entry)
          num_conc_from = aero_weight_array_num_conc(aero_state_from%awa, &
               aero_state_from%apa%particle(i_part), aero_data)
          num_conc_to = aero_weight_array_num_conc(aero_state_temp%awa, &
               aero_state_from%apa%particle(i_part), aero_data)
          mix_rat_from = num_conc_from / rho_from
          mix_rat_to = num_conc_to / rho_to
          if (mix_rat_to == mix_rat_from) then ! add and remove
             do_add = .true.
             do_remove = .true.
             i_transfer = i_transfer + 1
          elseif (mix_rat_to < mix_rat_from) then ! always add, maybe remove
             do_add = .true.
             do_remove = .false.
             if (pmc_random() < mix_rat_to / mix_rat_from) then
                do_remove = .true.
             end if
             ! Count adds
             i_transfer = i_transfer + 1
          else ! always remove, maybe add
             do_add = .false.
             if (pmc_random() < mix_rat_from / mix_rat_to) then
                do_add = .true.
             end if
             do_remove = .true.
             ! Count removes
             i_transfer = i_transfer + 1
          end if

          if (do_add) then
             call aero_state_add_particle(aero_state_temp, &
                  aero_state_from%apa%particle(i_part), aero_data)
             ! If we don't remove, we need to change the ID
             if (.not. do_remove) then
                call aero_particle_new_id(aero_state_temp%apa% &
                    particle(aero_state_temp%apa%n_part))
             end if
          end if
          if (do_remove) then
             if (removal_action /= AERO_INFO_NONE) then
                aero_info%id = aero_state_from%apa%particle(i_part)%id
                aero_info%action = removal_action
                call aero_state_remove_particle_with_info(aero_state_from, &
                     i_part, aero_info)
             else
                call aero_state_remove_particle_no_info(aero_state_from, &
                     i_part)
             end if
          end if
       end do
    end do

    call aero_state_add_particles(aero_state_to, aero_state_temp, aero_data)

  end subroutine aero_state_sample_particles_vertical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine aero_state_sample_particles_by_class(aero_state_from, &
       aero_state_to, aero_data, sample_prob, removal_action, &
       rrho_source, rrho_dest, dz_source, dz_dest, i_class)

    !> Original state.
    type(aero_state_t), intent(inout) :: aero_state_from
    !> Destination state.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Probability of sampling each particle.
    real(kind=dp), intent(in) :: sample_prob
    !> Action for removal (see pmc_aero_info module for action
    !> parameters). Set to AERO_INFO_NONE to not log removal.
    integer, intent(in) :: removal_action
    !> Inverse density of source aero_state.
    real(kind=dp), intent(in) :: rrho_source
    !> Inverse density of destination aero_state.
    real(kind=dp), intent(in) :: rrho_dest
    real(kind=dp), intent(in) :: dz_source
    real(kind=dp), intent(in) :: dz_dest
    !> Weight class to sample.
    integer, intent(in) :: i_class

    type(aero_state_t) :: aero_state_temp
    integer :: n_transfer, i_transfer
    integer :: i_group
    logical :: do_add, do_remove
    real(kind=dp) :: num_conc_from, num_conc_to
    type(aero_info_t) :: aero_info
    integer :: n_parts_before
    real(kind=dp) :: factor
    integer :: i_entry, i_part
    real(kind=dp) :: ratio, rho_to, rho_from, mix_rat_to, mix_rat_from
    integer :: i,j,k
    integer :: n_added, n_removed, n_part_from
    logical :: first_problem

    call assert(721006963, (sample_prob >= 0d0) .and. (sample_prob <= 1d0))
#ifndef PMC_USE_WRF
    call aero_state_zero(aero_state_to)
    call aero_state_copy_weight(aero_state_from, aero_state_to)
#else
    call aero_state_zero(aero_state_temp)
    call aero_state_set_weight(aero_state_temp, aero_data, &
         AERO_STATE_WEIGHT_FLAT_SPECIFIED)
    call aero_state_copy_weight(aero_state_to, aero_state_temp)
    factor = dz_source / dz_dest ! env_state_from%height / env_state_to%height
    call aero_weight_array_scale(aero_state_temp%awa, 1.0d0/factor)
#endif

#ifdef PMC_DEBUG
    n_parts_before = aero_state_total_particles(aero_state_to)
#endif

    ! Do we have a valid sort on the from state?
    if (.not. aero_state_from%valid_sort) then
       call aero_state_sort(aero_state_from, aero_data)
    end if

    ! Loop over groups
    i_group = 1
    rho_to = 1.0d0 !/ rrho_dest
    rho_from = 1.0d0 !/ rrho_source
    ratio = (aero_state_temp%awa%weight(i_group, i_class)%magnitude / rho_to) / &
         (aero_state_from%awa%weight(i_group, i_class)%magnitude  / rho_from)

    n_part_from = aero_state_total_particles(aero_state_from, i_group, i_class)

    if (ratio >= 1.0d0) then
       n_transfer = rand_binomial(n_part_from, sample_prob)
    else
       n_transfer = rand_poisson(n_part_from * (sample_prob / ratio))
    end if

    i_transfer = 0
    n_removed = 0 
    n_added = 0
    first_problem = .true.
    do while (i_transfer < n_transfer)
       if (aero_state_total_particles( &
            aero_state_from,i_group,i_class) <= 0) exit
       i_entry = pmc_rand_int(integer_varray_n_entry( &
            aero_state_from%aero_sorted%group_class%inverse(i_group, &
            i_class)))
       i_part = aero_state_from%aero_sorted%group_class%inverse(i_group, &
            i_class)%entry(i_entry)
       num_conc_from = aero_weight_array_num_conc(aero_state_from%awa, &
            aero_state_from%apa%particle(i_part), aero_data)
       num_conc_to = aero_weight_array_num_conc(aero_state_temp%awa, &
            aero_state_from%apa%particle(i_part), aero_data)
       mix_rat_from = num_conc_from / rho_from
       mix_rat_to = num_conc_to / rho_to 
       if (mix_rat_to == mix_rat_from) then ! add and remove
          do_add = .true.
          do_remove = .true.
          i_transfer = i_transfer + 1
       ! ratio < 1
       elseif (mix_rat_to < mix_rat_from) then ! always add, maybe remove
          do_add = .true.
          do_remove = .false.
          if (pmc_random() < mix_rat_to / mix_rat_from) then
             do_remove = .true.
          end if
          ! Count adds
          i_transfer = i_transfer + 1
       else ! always remove, maybe add
          do_add = .false.
          if (pmc_random() < mix_rat_from / mix_rat_to) then
             do_add = .true.
          end if
          do_remove = .true.
          ! Count removes
          i_transfer = i_transfer + 1
       end if
       if (do_add) then
          call aero_state_add_particle(aero_state_temp, &
               aero_state_from%apa%particle(i_part), aero_data)
          ! If we don't remove, we need to change the ID
          if (.not. do_remove) then
             call aero_particle_new_id(aero_state_temp%apa% &
                 particle(aero_state_temp%apa%n_part))
          end if
          n_added = n_added + 1
       end if
       if (do_remove) then
          if (removal_action /= AERO_INFO_NONE) then
             aero_info%id = aero_state_from%apa%particle(i_part)%id
             aero_info%action = removal_action
             call aero_state_remove_particle_with_info(aero_state_from, &
                  i_part, aero_info)
          else
             call aero_state_remove_particle_no_info(aero_state_from, &
                  i_part)
          end if
          n_removed = n_removed + 1
       end if
    end do

    call aero_state_add_particles(aero_state_to, aero_state_temp, aero_data)

  end subroutine aero_state_sample_particles_by_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the scaling factor for a grid cell.
  real(kind=dp) function get_scaling_factor(grid, i, k, j)

    !> WRF grid.
    type(domain) :: grid
    !> Grid cell index in x direction.
    integer :: i
    !> Grid cell index in y direction.
    integer :: k
    !> Grid cell index in z direction.
    integer :: j

    real(kind=dp) :: dx, dy, dz

    dz = real(grid%z_at_w(i,k+1,j) - grid%z_at_w(i,k,j), kind=dp)
    dx = real((1.0d0/grid%rdx)/grid%msftx(i,j), kind=dp)
    dy = real((1.0d0/grid%rdy)/grid%msfty(i,j), kind=dp)

    get_scaling_factor = dx*dy*dz

  end function get_scaling_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the multinomial sample using sorted probabilities.
  subroutine aero_state_multisample(env_states, aero_data, aero_state, &
       delta_aero_states, is, ie, js, je, ks, ke, pmc_is, pmc_ie, &
       pmc_ks, pmc_ke, pmc_js, pmc_je, source_i, source_j, source_k)

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
    !> Start index in west-east direction.
    integer, intent(in) :: is
    !> End index in west-east direction.
    integer, intent(in) :: ie
    !> Start index in north-south direction.
    integer, intent(in) :: js
    !> End index in north-south direction.
    integer, intent(in) :: je
    !> Start index in vertical direction.
    integer, intent(in) :: ks
    !> End index in vertical direction.
    integer, intent(in) :: ke
    !> 3D array of environmental states.
    type(env_state_t), dimension(is:ie,ks:ke,js:je), &
         intent(in) :: env_states
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state to sample.
    type(aero_state_t), intent(inout) :: aero_state
    !> Delta aerosol states to move particles to.
    type(aero_state_t), dimension(is:ie,ks:ke,js:je), intent(inout) :: &
         delta_aero_states
    !> Index of source grid cell in west-east direction.
    integer, intent(in) :: source_i
    !> Index of source grid cell in north-south direction.
    integer, intent(in) :: source_j
    !> Index of source grid cell in vertical direction.
    integer, intent(in) :: source_k

    ! Local variables
    real(kind=dp) :: remainder_prob
    real(kind=dp) :: adjusted_prob
    real(kind=dp) :: prob
    real(kind=dp), parameter :: remainder_eps = 1.0d-8
    type(env_state_t) :: env_state
    real(kind=dp) :: dest_vol, source_vol
    integer :: i, j, k, i_class, prob_count
    integer :: dest_i, dest_j, dest_k
    real(kind=dp) :: probs(5+pmc_ke-pmc_ks)
    integer, allocatable :: perm(:)
    integer :: i_offset(5+pmc_ke-pmc_ks), j_offset(5+pmc_ke-pmc_ks), &
         k_offset(5+pmc_ke-pmc_ks)

    env_state = env_states(source_i,source_k,source_j)

    do i_class = 1, aero_weight_array_n_class(aero_state%awa)
       remainder_prob = 1.0d0
       prob_count = 0 
       probs = 0.0d0
       k = 0
       j = 0
       do i = -1,1,2
          prob = env_state%prob_advection(i,k,j,i_class) &
               + env_state%prob_diffusion(i,k,j,i_class)
          if (prob > 0.0d0) then
             prob_count = prob_count + 1
             probs(prob_count) = prob
             i_offset(prob_count) = i
             j_offset(prob_count) = j
             k_offset(prob_count) = k
          end if
       end do
       k = 0
       i = 0
       do j = -1,1,2
          prob = env_state%prob_advection(i,k,j,i_class) &
               + env_state%prob_diffusion(i,k,j,i_class)
          if (prob > 0.0d0) then
             prob_count = prob_count + 1
             probs(prob_count) = prob
             i_offset(prob_count) = i
             j_offset(prob_count) = j
             k_offset(prob_count) = k
          end if    
       end do
       i = 0
       j = 0 
       do k = ks,ke
          if (k .ne. source_k) then
             prob = env_state%prob_vert_diffusion(k,i_class)
             if (prob > 0.0d0) then
                prob_count = prob_count + 1
                probs(prob_count) = prob
                i_offset(prob_count) = i
                j_offset(prob_count) = j
                k_offset(prob_count) = k - source_k
             end if 
          end if
       end do

       ! Sort the array
       if (allocated(perm)) deallocate(perm)
       allocate(perm(prob_count))
       call real_sort(probs(1:prob_count), perm)
       ! Loop over array
       do i = 1,prob_count
          prob = probs(i)
          adjusted_prob = min(prob / remainder_prob, 1.0d0)
          if (remainder_prob >= remainder_eps) then
             dest_i = source_i + i_offset(perm(i))
             dest_j = source_j + j_offset(perm(i))
             dest_k = source_k + k_offset(perm(i))
             source_vol = env_states(source_i,source_k,source_j)%cell_volume
             dest_vol = env_states(dest_i,dest_k,dest_j)%cell_volume
             call aero_state_sample_particles_by_class(aero_state, &
                  delta_aero_states(dest_i,dest_k,dest_j), aero_data, &
                  adjusted_prob, 0, env_states(source_i,source_k,source_j)%rrho, &
                  env_states(dest_i,dest_k,dest_j)%rrho, source_vol, dest_vol, &
                  i_class)
             remainder_prob = remainder_prob - prob
          end if   
       end do
    end do

  end subroutine aero_state_multisample

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module wrf_pmc_trans_aero
