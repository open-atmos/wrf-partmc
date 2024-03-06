! Copyright (C) 2011-2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc_transport module.

!> Functions that link both aerosol and gas species transport.
module wrf_pmc_trans

  use pmc_aero_state
  use pmc_gas_state
  use pmc_gas_data
  use pmc_aero_data
  use pmc_scenario
  use pmc_env_state
  use pmc_nucleate
  use pmc_util
  use pmc_coagulation
  use pmc_output
  use pmc_rand
  use wrf_pmc_trans_aero
  use wrf_pmc_dep

  use module_wrf_error
  use module_domain_type
  use module_domain
  use module_model_constants
  use module_configure
  USE module_state_description

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Top level subroutine for aerosol and gas species transport.
  subroutine wrf_pmc_trans_driver(grid, config_flags, env_states, scenario, &
       aero_data, aero_states, gas_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
       pmc_ks, pmc_ke, global_nx, global_ny, global_nz)

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
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Scenario data
    type(scenario_t),dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: scenario
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: aero_states
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: gas_states
    !> East-west domain dimension size of WRF domain.
    integer, intent(in) :: global_nx
    !> North-south domain dimension size of WRF domain.
    integer, intent(in) :: global_ny
    !> Top-bottom domain dimension size of WRF domain.
    integer, intent(in) :: global_nz

    integer :: numgas
    integer :: ids,ide, jds,jde, kds,kde,              &
               ims,ime, jms,jme, kms,kme,              &
               its,ite, jts,jte, kts,kte,               &
               ips,ipe, jps,jpe, kps,kpe
    integer :: i,j,k, i_class
    integer :: i_start, i_end, j_start, j_end
    integer :: nx,ny,nz
    real(kind=dp) :: t1, t2

    if (grid%do_transport) then
#ifdef PMC_DEBUG
       call wrf_message("PartMC_trans: starting transport")
#endif
#ifdef PMC_DEBUG
       t1 = MPI_Wtime()
#endif
       ! Determine the transport probabilities
       ! FIXME: These are currently only the horizontal ones
       do j = pmc_js, pmc_je
       do k = pmc_ks, pmc_ke
       do i = pmc_is, pmc_ie
          env_states(i,k,j)%prob_advection = 0.0d0
          env_states(i,k,j)%prob_diffusion = 0.0d0
          call compute_advect_probs_wrf(grid, aero_states(i,k,j), aero_data, &
               env_states(i,k,j))
          env_states(i,k,j)%prob_vert_diffusion = 0.0d0
          call compute_diffusion_probs(grid, env_states(i,k,j))
          !call eliminate_boundary_influence(env_states(i,k,j), &
          !   global_nx, global_ny)
       end do
       end do
       end do

       ! Compute vertical diffusion probabilities
       ! Compute B for a small timestep
       ! Then B^N_subcycle
       ! We do this column wise
       if (config_flags%do_uniform) then
          i_start = pmc_is
          i_end = pmc_ie
          j_start = pmc_js
          j_end = pmc_je
       else if (config_flags%do_scm) then
          i_start = pmc_is
          i_end = pmc_ie
          j_start = pmc_js
          j_end = pmc_je
       else if (config_flags%periodic_x .and. config_flags%periodic_y) then
          i_start = pmc_is
          i_end = pmc_ie
          j_start = pmc_js
          j_end = pmc_je
       else
          i_start = max(pmc_is,2)
          i_end = min(pmc_ie,global_nx-1)
          j_start = max(pmc_js,2)
          j_end = min(pmc_je,global_ny-1)
       end if

       do j = j_start,j_end 
       do i = i_start,i_end
          call compute_vertical_probs(grid, &
               env_states(i,pmc_ks:pmc_ke,j), global_nz, &
               aero_weight_array_n_class(aero_states(i,pmc_ks,j)%awa), &
               config_flags%vertmix_onoff)
       end do
       end do

       do j = pmc_js, pmc_je
       do k = pmc_ks, pmc_ke
       do i = pmc_is, pmc_ie
          do i_class = 1,aero_weight_array_n_class(aero_states(i,k,j)%awa)
             call normalize_probs(env_states(i,k,j), i_class)
          end do
       end do
       end do
       end do

#ifdef PMC_DEBUG
       t2 = MPI_Wtime()
       print*, 'Probability time', t2-t1
#endif

#ifdef PMC_DEBUG
       call wrf_message("PartMC_trans: done computing probabilities")
#endif

       call pmc_mpi_barrier()
       t1 = MPI_Wtime()
       ! Transport
       call transport(grid, config_flags, env_states, scenario, aero_data, &
            aero_states, gas_states, pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
            global_nx, global_ny, global_nz)

       t2 = MPI_Wtime()
       call pmc_mpi_barrier()
       print*, 'time spent in transport total = ', t2-t1

#ifdef PMC_DEBUG
       call wrf_message("PartMC_trans: finished transport")
#endif
    end if

    ! Do deposition on lowest layer
    if (pmc_ks == 1) then
       if (grid%do_deposition) then

#ifdef PMC_DEBUG
          call wrf_message("PartMC_trans: starting deposition")
#endif
          call wrf_pmc_dep_dry_driver(grid, config_flags, &
               aero_data, aero_states, gas_states, env_states, &
               pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
               global_nx, global_ny, global_nz)

#ifdef PMC_DEBUG
          call wrf_message("PartMC_trans: finished deposition")
#endif
       end if
    end if

  end subroutine wrf_pmc_trans_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map fluxes from WRF to PartMC
  subroutine compute_advect_probs_wrf(grid, aero_state, aero_data, env_state)

    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environmental state.
    type(env_state_t), intent(inout) :: env_state

    integer :: i, j, k
    real(kind=dp) :: prob_sum, rho_mass, mix_rat
    integer :: i_class, i_group, n_parts_group_class, wrf_class
    real(kind=dp), parameter :: min_mix_rat = 0.0d-9

    env_state%prob_advection = 0.0d0

    i = env_state%cell_ix
    j = env_state%cell_iy
    k = env_state%cell_iz
    i_group = 1
    rho_mass = 1.0d0 / real(grid%alt(i,k,j),kind=dp)
    do i_class = 1,aero_weight_array_n_class(aero_state%awa)
       wrf_class = i_class + 1
       n_parts_group_class = aero_state_total_particles(aero_state, &
            i_group, i_class)
       mix_rat = n_parts_group_class * &
            aero_state%awa%weight(i_group, i_class)%magnitude &
            / rho_mass

       ! All fluxes are in the positive direction
       ! X direction
       if (grid%u_flux(i+1,k,j,wrf_class) > 0.0e0) then
          if (mix_rat > min_mix_rat) then
             env_state%prob_advection(1,0,0,i_class) = &
                  real(grid%u_flux(i+1,k,j,wrf_class),kind=dp) / mix_rat
          end if
       end if
       if (grid%u_flux(i,k,j,wrf_class) < 0.0e0) then
          if (mix_rat > min_mix_rat) then
             env_state%prob_advection(-1,0,0,i_class) = -1.0d0 &
                  * real(grid%u_flux(i,k,j,wrf_class),kind=dp) / mix_rat
          end if
       end if

       ! Y direction
       if (grid%v_flux(i,k,j+1,wrf_class) > 0.0e0) then
          if (mix_rat > min_mix_rat) then
             env_state%prob_advection(0,0,1,i_class) = & 
                  real(grid%v_flux(i,k,j+1,wrf_class),kind=dp) / mix_rat
          end if
       end if
       if (grid%v_flux(i,k,j,wrf_class) < 0.0e0) then
          if (mix_rat > min_mix_rat) then
             env_state%prob_advection(0,0,-1,i_class) = -1.0d0 &
                  * real(grid%v_flux(i,k,j,wrf_class),kind=dp) / mix_rat
          end if
       end if

       ! Z direction
       if (grid%w_flux(i,k+1,j,wrf_class) > 0.0e0) then
          if (mix_rat > min_mix_rat) then
            env_state%prob_advection(0,1,0,i_class) = &
                  min(real(grid%w_flux(i,k+1,j,wrf_class),kind=dp) / mix_rat, .95d0)
          end if
       end if
       if (grid%w_flux(i,k,j,wrf_class) < 0.0e0) then
          if (mix_rat > min_mix_rat) then
             env_state%prob_advection(0,-1,0,i_class) = min(-1.0d0 &
                  * real(grid%w_flux(i,k,j,wrf_class),kind=dp) / mix_rat, .95d0)
          end if
       end if

       prob_sum = sum((env_state%prob_advection(:,:,:,i_class)))
       if (prob_sum > 1.0d0) then
          env_state%prob_advection(:,:,:,i_class) = &
               env_state%prob_advection(:,:,:,i_class) / (prob_sum + 1.0d-8)
       end if
    end do

  end subroutine compute_advect_probs_wrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute probabilities from upwind method fluxes
  subroutine compute_advect_probs(grid, env_state)

    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> Environmental state.
    type(env_state_t), intent(inout) :: env_state

    integer :: i, j, k
    real(kind=dp) :: dx, dy, dt, dz
    real(kind=dp) :: vel, flux

    env_state%prob_advection = 0.0d0

    i = env_state%cell_ix
    j = env_state%cell_iy
    k = env_state%cell_iz

    dx = grid%dx
    dy = grid%dy
    dz = - 1.0d0 * grid%dnw(k)
    dt = grid%dt

    ! X direction
    ! F i+1/2
    vel = grid%u_2(i+1,k,j) 
    flux = (dt * (vel + abs(vel)))/(2.0d0*dx)
    env_state%prob_advection(1,0,0,:) = flux

    ! F i-1/2
    vel = grid%u_2(i,k,j) 
    flux = -(dt * (vel - abs(vel)))/(2.0d0*dx)
    env_state%prob_advection(-1,0,0,:) = flux

    ! Y direction
    ! F j + 1/2
    vel = grid%v_2(i,k,j+1) 
    flux = (dt * (vel + abs(vel)))/(2.0d0*dy)
    env_state%prob_advection(0,0,1,:) = flux

    ! F j - 1/2
    vel = grid%v_2(i,k,j)
    flux = -(dt * (vel - abs(vel)))/(2.0d0*dy)
    env_state%prob_advection(0,0,-1,:) = flux

    ! Vertical advection
    if (k == 1) then
       vel =  -1.0d0 * grid%ww(i,k+1,j) / grid%mut(i,j)
       flux = (dt * (vel + abs(vel)))/(2.0d0*dz)
       env_state%prob_advection(0,1,0,:) = flux
    else if (k == 39) then
       vel =  - 1.0d0 * grid%ww(i,k,j) / grid%mut(i,j)
       flux = -(dt * (vel - abs(vel)))/(2.0d0*dz)
       env_state%prob_advection(0,-1,0,:) = flux
    else
       vel = - 1.0d0 * grid%ww(i,k+1,j) / grid%mut(i,j)
       flux = (dt * (vel + abs(vel)))/(2.0d0*dz)
       env_state%prob_advection(0,1,0,:) = flux
       vel =  - 1.0d0 * grid%ww(i,k,j) / grid%mut(i,j)
       flux = -(dt * (vel - abs(vel)))/(2.0d0*dz)
       env_state%prob_advection(0,-1,0,:) = flux
    end if

  end subroutine compute_advect_probs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute diffusion probabilities using second order.
  subroutine compute_diffusion_probs(grid, env_state)

    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> Environmental state.
    type(env_state_t), intent(inout) :: env_state

    integer :: i, j, k
    real(kind=dp) :: dx, dy, dt
    real(kind=dp) :: diff_coef, flux
    real(kind=dp) :: rho_edge, rho_mass
    ! Set to have a minimum value for diffusion.
    real(kind=dp), parameter :: min_diff = 0.0d0

    env_state%prob_diffusion = 0.0d0

    dx = grid%dx
    dy = grid%dy
    dt = grid%dt 
    i = env_state%cell_ix
    j = env_state%cell_iy
    k = env_state%cell_iz

    rho_mass = (1.0d0/grid%alt(i,k,j))

    ! X direction
    ! F i+1/2
    diff_coef = max(.5d0 * (grid%xkhh(i,k,j) + grid%xkhh(i+1,k,j)), min_diff)
    rho_edge = .5d0*(rho_mass + (1.0d0/grid%alt(i+1,k,j)))
    flux = (dt * diff_coef * rho_edge) / (dx * dx * rho_mass) 
    env_state%prob_diffusion(1,0,0,:) = flux

    ! F i-1/2
    diff_coef = max(.5d0 * (grid%xkhh(i,k,j) + grid%xkhh(i-1,k,j)), min_diff)
    rho_edge = .5d0*(rho_mass + (1.0d0/grid%alt(i-1,k,j)))
    flux = (dt * diff_coef * rho_edge) / (dx * dx * rho_mass)
    env_state%prob_diffusion(-1,0,0,:) = flux

    ! Y direction
    ! F j + 1/2
    diff_coef = max(.5d0 * (grid%xkhh(i,k,j) + grid%xkhh(i,k,j+1)), min_diff)
    rho_edge = .5d0*(rho_mass + (1.0d0/grid%alt(i,k,j+1)))
    flux = (dt * diff_coef * rho_edge) / (dy * dy * rho_mass)
    env_state%prob_diffusion(0,0,1,:) = flux

    ! F j - 1/2
    diff_coef = max(.5d0 * (grid%xkhh(i,k,j) + grid%xkhh(i,k,j-1)), min_diff)
    rho_edge = .5d0*(rho_mass + (1.0d0/grid%alt(i,k,j-1)))
    flux = (dt * diff_coef * rho_edge) / (dy * dy * rho_mass)
    env_state%prob_diffusion(0,0,-1,:) = flux

  end subroutine compute_diffusion_probs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute vertical probabilities.
  subroutine compute_vertical_probs(grid, env_states, nz, n_class, &
       do_vertmix)

    !> WRF grid.
    type(domain), intent(inout) :: grid
    !> Environmental states.
    type(env_state_t), dimension(nz), intent(inout) :: env_states
    !> Number of vertical levels.
    integer, intent(in) :: nz
    !> Number of classes.
    integer, intent(in) :: n_class
    !> Vertical mixing flag due to turbulence.
    integer, intent(in) :: do_vertmix

    integer :: t, k
    integer :: num_steps, i_class
    real(kind=dp) :: small_dt, dt

    real(kind=dp), dimension(nz,nz) :: B
    real(kind=dp), dimension(nz,nz) :: B_tmp, R
    real(kind=dp), dimension(nz,nz) :: A

    ! Temp arrays
    real(kind=dp), dimension(nz) :: z_mass, rho_mass
    real(kind=dp), dimension(nz+1) :: z_w,k_w,rho_w

    real(kind=dp) :: tmp
    integer :: pmc_i, pmc_j, kk

    call stable_timestep(grid, env_states(1)%cell_ix, env_states(1)%cell_iy, &
         small_dt, num_steps)
  
    pmc_i = env_states(1)%cell_ix
    pmc_j = env_states(1)%cell_iy

    ! Geometric height, diffusion coef. and density values at mass point
    do k = 1,nz
       z_mass(k) = real(grid%z(pmc_i,k,pmc_j),kind=dp)
       rho_mass(k) = (1.0d0/real(grid%alt(pmc_i,k,pmc_j),kind=dp))
    end do

    ! Vertical column geometric heights at cell boundaries
    do k =1,nz+1
       z_w(k) = real(grid%z_at_w(pmc_i,k,pmc_j),kind=dp)
    end do

    ! Compute K on the w points (boundaries between cells only)
    k_w = 0.0d0
    do k = 1,nz-1
       k_w(k)=max(1.0d-6,real(grid%exch_h(pmc_i,k+1,pmc_j),kind=dp))
    end do
    ! Compute inverse density on the w points (boundaries between cells)
    do k = 1,nz-1
       tmp = real(z_mass(k+1) - z_mass(k),kind=dp)
       rho_w(k) = ((z_mass(k+1)-z_w(k+1))/tmp)*rho_mass(k) &
            + ((z_w(k+1)-z_mass(k))/tmp)*rho_mass(k+1)
    end do

    ! Construct the B matrix for the column - this is tridiagonal
    B = 0.0d0
    ! Interior grid cells
    do k = 2,nz-1
       ! Down
       B(k,k-1) = (small_dt/((z_w(k+1)-z_w(k))*(z_mass(k) &
          -z_mass(k-1))))*(rho_w(k-1)/rho_mass(k))*k_w(k-1) 
       ! Up
       B(k,k+1) = (small_dt/((z_w(k+1)-z_w(k))*(z_mass(k+1) &
            -z_mass(k))))*(rho_w(k)/rho_mass(k))*k_w(k)
       ! Stays
       B(k,k) = 1.0d0 - B(k,k-1) - B(k,k+1)
    end do

    ! Boundary conditions at surface 
    B(1,2) = (small_dt/((z_w(2)-z_w(1)) &
         *(z_mass(2)-z_mass(1))))*(rho_w(1)/rho_mass(1))*k_w(1)
    B(1,1) = 1.0d0 - B(1,2)
    ! Boundary conditions at top
    B(nz,nz-1) = (small_dt/((z_w(nz+1)-z_w(nz)) &
       *(z_mass(nz)-z_mass(nz-1))))*(rho_w(nz-1)/rho_mass(nz))*k_w(nz-1)
    B(nz,nz) = 1.0d0 - B(nz,nz-1)

    ! Multiply the B matrix to get full time step probabilities as B*
    B_tmp = B
    do t = 2, num_steps 
       R = matmul(B_tmp, B)
       B_tmp = R
       ! Here we could throw away values to reduce the stencil size
       ! For now, we will just normalize
       do k = 1, nz
           B_tmp(k,1:nz) = B_tmp(k,1:nz) / sum(B_tmp(k,1:nz)) 
       end do
    end do

    ! Construct Advection matrix
    do i_class = 1, n_class
       A = 0.0d0
       do k = 2,nz-1
          A(k,k) = 1.0d0 - (env_states(k)%prob_advection(0,1,0,i_class) + &
             env_states(k)%prob_advection(0,-1,0,i_class))
          A(k,k-1) = env_states(k)%prob_advection(0,-1,0,i_class)
          A(k,k+1) = env_states(k)%prob_advection(0,1,0,i_class)
       end do
       A(1,1) = 1.0d0 - env_states(1)%prob_advection(0,1,0,i_class)
       A(1,2) = env_states(1)%prob_advection(0,1,0,i_class)

       A(nz,nz) = 1.0d0 - env_states(nz)%prob_advection(0,-1,0,i_class)
       A(nz,nz-1) = env_states(nz)%prob_advection(0,-1,0,i_class)

       ! Construct final matrix
       if (do_vertmix > 0) then
          R = matmul(B_tmp, A)
       else
          R = A
       end if

       ! Give the correct row of the B* matrix as a vector to an env_state
       ! representing the probabilities of particles moving 
       ! FIXME: This is no longer only diffusion
       do k = 1, nz
          env_states(k)%prob_vert_diffusion(1:nz,i_class) = R(k,1:nz)
          ! Test to investigate if anything is wrong
          if (B_tmp(k,k) > 1.0d0) then
             print*, pmc_i, k, pmc_j
             print*, B_tmp(k,1:nz)
             print*, env_states(k)%prob_vert_diffusion(1:nz,i_class)
             print*, k_w
             print*, small_dt, num_steps
             print*, z_mass
             print*, z_w
             print*, rho_mass
             print*, rho_w
             ! Dump full state
             do kk = 1,nz
                print*, A(kk,1:nz)
             end do
             do kk = 1,nz
                print*, B(kk,1:nz)
             end do
             do kk = 1,nz
                print*, B_tmp(kk,1:nz)
             end do
             do kk = 1,nz
                print*, R(kk,1:nz)
             end do
          end if
       end do
    end do

  end subroutine compute_vertical_probs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calls aerosol transport.
  subroutine transport(grid, config_flags, env_states, scenario, aero_data, &
       aero_states, gas_states, pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke, &
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
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: env_states
    !> Scenario data.
    type(scenario_t),dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: scenario
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: aero_states
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
        intent(inout) :: gas_states
    !> East-west domain dimension size of WRF domain.
    integer, intent(in) :: global_nx
    !> North-south domain dimension size of WRF domain.
    integer, intent(in) :: global_ny
    !> Top-bottom domain dimension size of WRF domain.
    integer, intent(in) :: global_nz

#ifdef PMC_DEBUG
       call wrf_message("PartMC_trans: starting aerosol transport")
#endif

    ! Advection-diffusion of aerosols
    call trans_aero(grid, config_flags, aero_data, scenario, env_states, &
         aero_states, pmc_is, pmc_ie, pmc_ks, pmc_ke, pmc_js, pmc_je, &
         global_nx, global_ny, global_nz)

#ifdef PMC_DEBUG
       call wrf_message("PartMC_trans: completed aerosol transport")
#endif

  end subroutine transport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes a stable timestep for vertical diffusion.
  subroutine stable_timestep(grid, tar_x, tar_y, small_timestep, num_substeps)

    !> WRF grid.
    type(domain), intent(in) :: grid
    !> Grid cell target x.
    integer, intent(in) :: tar_x
    !> Grid cell target y.
    integer, intent(in) :: tar_y
    !> Desired time step for transport.
    real(kind=dp), intent(out) :: small_timestep 
    !> Number of diffusion time steps required.
    integer,intent(out) :: num_substeps

    real(kind=dp) :: factor
    real(kind=dp) :: min_dz,max_k,max_time
    real(kind=dp) :: dt

    dt = grid%dt

    min_dz = grid%z(tar_x,2,tar_y) - grid%z(tar_x,1,tar_y)
    max_k = max(real(maxval(grid%exch_h(tar_x,:,tar_y)),kind=dp),1.0d-6)
    factor = .1d0
    max_time = factor*(min_dz**2.0d0)/max_k
    num_substeps = ceiling(dt/max_time)
    small_timestep = dt/num_substeps

  end subroutine stable_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eliminate_boundary_influence(env_state, global_nx, global_ny)

    type(env_state_t), intent(inout) :: env_state
    integer, intent(in) :: global_nx
    integer, intent(in) :: global_ny

    integer :: i, j

    i = env_state%cell_ix
    j = env_state%cell_iy

    if (i == 1 .or. i == global_nx .or. j == 1 .or. j ==global_ny) then
       env_state%prob_advection = 0.0d0
       env_state%prob_diffusion = 0.0d0
    end if
    ! If we want reflection at the boundary
    if (i == 2) then
       env_state%prob_advection(-1,0,0,:) = 0.0d0
       env_state%prob_diffusion(-1,0,0,:) = 0.0d0
    elseif ( i == global_nx-1) then
       env_state%prob_advection(1,0,0,:) = 0.0d0
       env_state%prob_diffusion(1,0,0,:) = 0.0d0
    end if
    if (j == 2) then
       env_state%prob_advection(0,0,-1,:) = 0.0d0
       env_state%prob_diffusion(0,0,-1,:) = 0.0d0
    else if (j == global_ny-1) then
       env_state%prob_advection(0,0,1,:) = 0.0d0
       env_state%prob_diffusion(0,0,1,:) = 0.0d0
    end if

  end subroutine eliminate_boundary_influence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Normalize the probabilities to be no greater than 1
  subroutine normalize_probs(env_state, i_class)
    type(env_state_t), intent(inout) :: env_state
    integer, intent(in) :: i_class

    real(kind=dp) :: prob_sum
    real(kind=dp) :: prob_epsilon
    integer :: i, j, k
    logical :: valid_probs

    ! Check to see if any values are negative
    valid_probs = all(env_state%prob_advection(:,:,:,i_class) >= 0.0d0)
    valid_probs = valid_probs .and. &
         all(env_state%prob_diffusion(:,:,:,i_class) >= 0.0d0)
    valid_probs = valid_probs .and. &
         all(env_state%prob_vert_diffusion(:,i_class) >= 0.0d0)
    if (.not. valid_probs) then
       print*, env_state%cell_iz, i_class
       print*, env_state%prob_advection(:,:,:,i_class) 
       print*, env_state%prob_diffusion(:,:,:,i_class) 
       print*, env_state%prob_vert_diffusion(:,i_class)
    end if
    call assert(414897882, valid_probs)

    prob_sum = sum((env_state%prob_advection(:,0,:,i_class)))
    prob_sum = prob_sum + sum(env_state%prob_diffusion(:,:,:,i_class))

    do k = 1,size(env_state%prob_vert_diffusion,1)
       if (k .ne. env_state%cell_iz) then
          prob_sum = prob_sum + env_state%prob_vert_diffusion(k,i_class)
       end if
    end do

    if (prob_sum >= 1.0d0) then
       env_state%prob_vert_diffusion(:,i_class) = &
            env_state%prob_vert_diffusion(:,i_class) / (prob_sum)
       env_state%prob_diffusion(:,:,:,i_class) = &
            env_state%prob_diffusion(:,:,:,i_class) / (prob_sum)
       env_state%prob_advection(:,0,:,i_class) = &
            env_state%prob_advection(:,0,:,i_class) / (prob_sum)
    end if


    ! Check if any values are less than zero
    valid_probs = all(env_state%prob_advection(:,:,:,i_class) >= 0.0d0)
    valid_probs = valid_probs .and. &
         all(env_state%prob_diffusion(:,:,:,i_class) >= 0.0d0)
    valid_probs = valid_probs .and. &
         all(env_state%prob_vert_diffusion(:,i_class) >= 0.0d0)
    if (.not. valid_probs) then
       print*, env_state%cell_iz, i_class
       print*, env_state%prob_advection(:,:,:,i_class)
       print*, env_state%prob_diffusion(:,:,:,i_class)
       print*, env_state%prob_vert_diffusion(:,i_class)
    end if
    call assert(184874367, valid_probs)
    ! Check is any values are greater than one
    valid_probs = all(env_state%prob_advection(:,:,:,i_class) <= 1.0d0)
    valid_probs = valid_probs .and. &
         all(env_state%prob_diffusion(:,:,:,i_class) <= 1.0d0)
    valid_probs = valid_probs .and. &
         all(env_state%prob_vert_diffusion(:,i_class) <= 1.0d0)
    if (.not. valid_probs) then
       print*, env_state%cell_iz, i_class
       print*, env_state%prob_advection(:,:,:,i_class)
       print*, env_state%prob_diffusion(:,:,:,i_class)
       print*, env_state%prob_vert_diffusion(:,i_class)
    end if
    call assert(184874369, valid_probs)

    ! Check that sum of probablities is less than 1 + 1e-8
    prob_sum = sum((env_state%prob_advection(:,0,:,i_class)))
    prob_sum = prob_sum + sum(env_state%prob_diffusion(:,:,:,i_class))

    do k = 1,size(env_state%prob_vert_diffusion,1)
       if (k .ne. env_state%cell_iz) then
          prob_sum = prob_sum + env_state%prob_vert_diffusion(k,i_class)
       end if
    end do

    call assert(414897881, prob_sum < 1.0d0 + 1.0d-8)

  end subroutine normalize_probs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module wrf_pmc_trans
