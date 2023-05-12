! Copyright (C) 2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc_dep_aero module.

!> Aerosol particle deposition.
!> FIXME: Kinematic and dynamic visc. are assumed constant.
!> FIXME: g should be made a const variable and used here
!> FIXME: Give A and gamma descriptions
module wrf_pmc_dep_aero

  use pmc_aero_state
  use pmc_aero_data
  use pmc_constants
  use pmc_aero_particle
  use pmc_rand
  use pmc_env_state

  USE module_domain, ONLY : domain

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Driver for aerosol dry deposition.
  subroutine dep_aero_dry_driver(grid, aero_states, aero_data, env_states, ustar, &
       nx, ny, dt, aer_res_a, pmc_is, pmc_ie, pmc_js, pmc_je, pmc_ks, pmc_ke)

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
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: aero_states
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: env_states
    !> Friction velocity (m/s).
    real(kind=dp), dimension(pmc_is:pmc_ie,pmc_js:pmc_je) :: ustar
    !> East-west dimension.
    integer :: nx
    !> North-south dimension.
    integer :: ny
    !> Timestep (s).
    real(kind=dp) :: dt
    !> Aerodynamic resistance (s/m).
    real(kind=dp), dimension(pmc_is:pmc_ie,pmc_js:pmc_je) :: aer_res_a
    integer :: i, j
    ! From WRF-Chem to map landuse
    integer, save  :: isurftype_zhang_from_usgs(25) ! maps usgs surface type to zhang type
    data isurftype_zhang_from_usgs / &  
               ! zhang type         <---  usgs type                    
        15, &  ! urban                    1: urban and built-up land
         7, &  ! crops, mixed farming     2: dryland cropland & pasture
         7, &  !                          3: irrigated crop & pasture
         7, &  !                          4: mix dry/irrig crop & pasture
         7, &  !                          5: cropland/grassland mosaic 
        10, &  ! shrubs & woodland        6: cropland/woodland mosaic 
         6, &  ! grass                    7: grassland             
         6, &  !                          8: shrubland            
         6, &  !                          9: mixed shrubland/grassland  
         6, &  !                         10: savanna                   
         4, &  ! deciduous broadleaf     11: deciduous broadleaf forest
         3, &  ! deciduous needle-lf     12: deciduous needleleaf forest
         2, &  ! evergreen broadleaf     13: evergreen broadleaf forest
         1, &  ! evergreen needle-lf     14: evergreen needleleaf forest
         5, &  ! mixed broad/needle-lf   15: mixed forest          
        14, &  ! ocean                   16: water bodies      
        11, &  ! wetland w plants        17: herbaceous wetland
        10, &  !                         18: wooded wetland   
         8, &  ! desert                  19: barren or sparsely vegetated
         9, &  ! tundra                  20: herbaceous tundra         
        10, &  !                         21: wooded tundra            
         9, &  !                         22: mixed tundra            
         8, &  !                         23: bare ground tundra     
        12, &  ! ice cap & glacier       24: snow or ice           
         8  /  !                         25: no data              

    ! From Zhang 2001
    real(kind=dp) :: gamma(15)     ! exponent of schmidt number
    data gamma/ 0.56, 0.58, 0.56, 0.56, 0.56, &
                   0.54, 0.54, 0.54, 0.54, 0.54, &
                   0.54, 0.54, 0.50, 0.50, 0.56/

    real(kind=dp) :: alpha(15)      ! parameter for impaction
    data alpha/   1.0,   0.6,   1.1,   0.8,   0.8, &
                    1.2,   1.2,  50.0,  50.0,   1.3, &
                    2.0,  50.0, 100.0, 100.0,   1.5/

    real(kind=dp) :: A(5,15) ! radius (m) of surface collectors
    data A / &
         0.002, 0.002, 0.002, 0.002, 0.002, &  ! luc 1
         0.005, 0.005, 0.005, 0.005, 0.005, &  ! luc 2
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 3
         0.005, 0.005, 0.010, 0.010, 0.005, &  ! luc 4
         0.005, 0.005, 0.005, 0.005, 0.005, &  ! luc 5
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 6
         0.002, 0.002, 0.005, 0.005, 0.002, &  ! luc 7
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 8
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 9
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 10
         0.010, 0.010, 0.010, 0.010, 0.010, &  ! luc 11
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 12
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 13
         9.999, 9.999, 9.999, 9.999, 9.999, &  ! luc 14
         0.010, 0.010, 0.010, 0.010, 0.010  /  ! luc 15

    integer :: isurftype, iseason, isurftype_in, modesurftype
    integer :: k

    k = 1
    modesurftype = 2
    do i = max(pmc_is,2),min(pmc_ie, nx-1)
    do j = max(pmc_js,2),min(pmc_je, ny-1)
       isurftype = 7 
       isurftype_in = grid%ivgtyp(i,j)
       if (modesurftype <= 1) then
          if ((isurftype_in >= 1) .and. (isurftype_in <= 15)) &
             isurftype = isurftype_in
       else
          if ((isurftype_in >= 1) .and. (isurftype_in <= 25)) &
             isurftype = isurftype_zhang_from_usgs(isurftype_in)
       end if
       iseason = 1
!       if ((iseason_in >= 1) .and. (iseason_in <= 5)) &
!         iseason = iseason_in
!       end if

       call dry_dep_aero_state(aero_states(i,k,j), aero_data, &
            env_states(i,k,j), aer_res_a(i,j), ustar(i,j), alpha(isurftype), &
            gamma(isurftype), A(iseason,isurftype), dt)
    end do
    end do

  end subroutine dep_aero_dry_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove particles from an aero state by dry deposition.
  subroutine dry_dep_aero_state(aero_state, aero_data, env_state, aer_res_a, &
       ustar, alpha, gamma, A, dt)

    !> Aerosol states.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment states.
    type(env_state_t), intent(in) :: env_state
    !> Aerodynamic resistance (s/m).
    real(kind=dp), intent(in) :: aer_res_a
    !> Friction velocity (m/s).
    real(kind=dp), intent(in) :: ustar
    !>
    real(kind=dp), intent(in) :: alpha
    !> Parameter for Brownian collection efficiency.
    real(kind=dp), intent(in) :: gamma
    !> Characteristic radius of collectors (m).
    real(kind=dp), intent(in) :: A
    !> Timestep (s).
    real(kind=dp), intent(in) :: dt

    ! Deposition velocity array
    real(kind=dp), allocatable :: vd(:)
    ! Particle removal probability array
    real(kind=dp), allocatable :: remove_prob(:)

    allocate(vd(aero_state%apa%n_part))
    allocate(remove_prob(aero_state%apa%n_part))

    vd = 0.0d0
    remove_prob = 0.0d0

    ! Compute deposition velocity array
    call compute_dep_vel(aero_state, aero_data, aer_res_a, ustar, &
         env_state%temp, env_state%rrho, alpha, gamma, A, vd)
    ! Compute the removal probability array
    call compute_dep_prob(vd, dt, env_state%height, remove_prob)
    ! Error checking
    if ((maxval(remove_prob) >= 1.0d0).or.(minval(remove_prob) <= 0.0d0)) then
       print*, "Unrealistic deposition", minval(remove_prob), maxval(remove_prob)
    end if
    ! Remove particles based on probabilities
    call aero_particles_remove_by_dep(aero_state, remove_prob)

    deallocate(vd)
    deallocate(remove_prob)

  end subroutine dry_dep_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute deposition velocities.
  subroutine compute_dep_vel(aero_state, aero_data, aer_res_a, ustar, temp, &
       rrho, alpha, gamma, A, vd)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerodynamic resistance (s/m)
    real(kind=dp), intent(in) :: aer_res_a
    !> Friction velocity (m/s).
    real(kind=dp), intent(in) :: ustar
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Inverse density
    real(kind=dp), intent(in) :: rrho
    !> Dry deposition velocity (m/s).
    real(kind=dp), intent(inout) :: vd(:)
    !>
    real(kind=dp), intent(in) :: alpha
    !> Parameter for Brownian collection efficiency.
    real(kind=dp), intent(in) :: gamma
    !> Characteristic radius of collectors (m).
    real(kind=dp), intent(in) :: A

    integer :: i_part
    real(kind=dp),allocatable :: diam(:)
    real(kind=dp),allocatable :: dens(:)
    real(kind=dp) :: rs, vs

    ! Compute properties of the aerosols
    allocate(diam(aero_state%apa%n_part))
    allocate(dens(aero_state%apa%n_part))

    do i_part = 1, aero_state%apa%n_part
       diam(i_part) = aero_particle_diameter(aero_state%apa%particle(i_part), &
            aero_data)
       dens(i_part) = aero_particle_density(aero_state%apa%particle(i_part), &
            aero_data)
    end do

    do i_part = 1,aero_state%apa%n_part
       vs = calculate_vs(diam(i_part), dens(i_part))
       rs = calculate_rs(diam(i_part), vs, temp, rrho, ustar, alpha, gamma, A)
       vd(i_part) = calculate_vd(aer_res_a, rs, vs)
    end do

    deallocate(diam)
    deallocate(dens)

  end subroutine compute_dep_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the probability of each particle being removed by dry deposition.
  subroutine compute_dep_prob(vd, dt, dz, remove_prob)

    !> Array of dry deposition velocities (m/s).
    real(kind=dp), intent(in) :: vd(:)
    !> Time step (s).
    real(kind=dp), intent(in) :: dt
    !> Height of lowest level box (m).
    real(kind=dp), intent(in) :: dz
    !> Probabilities of each particle being removed.
    real(kind=dp), intent(inout) :: remove_prob(:)

    integer :: i_part

    ! Two options:
    ! 1.) prob = (dt / dz) * vd
    ! 2.) prob = 1d0 - exp(-dt * vd / dz)
    do i_part = 1,size(vd)
       remove_prob(i_part) = (dt / dz) * vd(i_part)
!       remove_prob(i_part) = 1.0d0 - exp(-dt * vd(i_part) / dz)
    end do

  end subroutine compute_dep_prob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests the particle array to see if the particle is to be removed by
  !> dry deposition based on the provided probability.
  subroutine aero_particles_remove_by_dep(aero_state, remove_prob)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Probabilities of each particle being removed.
    real(kind=dp), intent(in) :: remove_prob(:)

    integer :: i_part

    do i_part = aero_state%apa%n_part,1,-1
         if (pmc_random() < remove_prob(i_part)) then
            call aero_state_remove_particle_no_info(aero_state, i_part)
         end if
    end do

  end subroutine aero_particles_remove_by_dep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the particle settling velocity.
  !> Seinfeld and Pandis eq. 19.18.
  real(kind=dp) function calculate_vs(diameter, density)

     !> Particle diameter.
     real(kind=dp) :: diameter
     !> Particle density.
     real(kind=dp) :: density

     real(kind=dp) :: c_c

     c_c = slip_correction_factor(diameter)

     calculate_vs = ((diameter)**2.0d0 * density * const%grav * c_c) &
          / (18.0d0 * const%air_dyn_visc)

  end function calculate_vs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate surface resistance of a given particle.
  real(kind=dp) function calculate_rs(diameter, vs, temperature, rrho, ustar, &
       alpha, gamma, A)

    !> Particle diameter (m).
    real(kind=dp) :: diameter
    !> Settling velocity (m/s).
    real(kind=dp) :: vs
    !> Temperature (K).
    real(kind=dp) :: temperature
    !> Inverse density
    real(kind=dp) :: rrho
    !> Friction velocity (m/s).
    real(kind=dp) :: ustar
    !> Alpha value.
    real(kind=dp) :: alpha
    !> Parameter for Brownian collection efficiency.
    real(kind=dp) :: gamma
    !> Characteristic radius of collectors (m).
    real(kind=dp) :: A

    real(kind=dp) :: eps_0
    real(kind=dp) :: diff
    real(kind=dp) :: nu
    real(kind=dp) :: EB,EIM,EIN,R1
    real(kind=dp) :: beta
    real(kind=dp) :: Sc
    real(kind=dp) :: st
    real(kind=dp) :: C_b, C_IM

    logical, parameter :: updated_version = .true.

    ! Compute Schmidt number
    diff = compute_brownian_diff(diameter, temperature)
    Sc = 1.45d-5 / diff

     if (.not. updated_version) then
        nu = 2.0
        C_b = 1.0
        beta = 2.0
        C_IM = 1.0
    else
       gamma = 2.0/3.0
       C_b = .2
       nu = 0.8
       beta = 1.7
       C_IM =  0.4
    end if

    ! Compute Stokes number
    St = calculate_st(vs, rrho, ustar, A)

    ! EB: Collection efficiency from Brownian diffusion
    ! Equation 6 from Zhang et al. (2001)
    EB = C_b * Sc**(-gamma)

    ! EIM: Collection efficiency from impaction
    ! Equation 7b from Zhang et al. (2001)
    beta = 2.0d0
    EIM = C_IM * (St / (alpha + St))**(beta)

    ! EIN: Collection efficiency from interception
    ! Equation 8 from Zhang et al. (2001)
    EIN = .5 * (diameter / A)**nu

    ! R1: Correction factor for sticking
    ! Equation 9 from Zhang et al. (2001)
    R1 = exp(-St**.5)

    ! Empirical constant
    eps_0 = 3.0d0

    ! Equation 5 from Zhang et al. (2001)
    calculate_rs = 1.0d0 / (eps_0 * ustar * (EB + EIM + EIN) * R1)

  end function calculate_rs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the dry deposition velocity.
  !> Seinfeld and Pandis eq. 19.7.
  real(kind=dp) function calculate_vd(ra, rs, vs)

    !> Aerodynamic resistance (s/m)
    real(kind=dp) :: ra
    !> Surface resistance (s/m).
    real(kind=dp) :: rs
    !> Settling velocity (m/s).
    real(kind=dp) :: vs

    real(kind=dp) :: tmp

    tmp = ra + rs + ra * rs * vs

    calculate_vd = (1.0d0 / tmp) + vs

  end function calculate_vd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the Cunningham slip correction factor for a single particle.
  !> Seinfeld and Pandis eq. 9.34.
  real(kind=dp) function slip_correction_factor(D_p)

     !> Particle diameter (m).
     real(kind=dp) :: D_p

     real(kind=dp) :: lambda

     ! Mean free path (m)
     lambda = .065d-6

     slip_correction_factor = 1.0d0 + ((2.0d0 * lambda)/D_p) * ( 1.257 &
         + .4 * exp(-(1.1 * D_p)/(2.0d0 * lambda)))

  end function slip_correction_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the Stokes number.
  !> Zhang et al (2001) for vegetative surface.
  real(kind=dp) function calculate_st(vs, rrho, ustar, A)

   !> Settling velocity (m/s).
   real(kind=dp) :: vs
   !> Inverse density
   real(kind=dp) :: rrho
   !> Friction velocity (m/s).
   real(kind=dp) :: ustar
   !> Characteristic radius of collectors (mm).
   real(kind=dp) :: A

   if (A > 1.0d0) then
      calculate_st = (vs * ustar**2) / (const%air_dyn_visc * rrho)
   else
      calculate_st = (vs * ustar) / (const%grav * A)
   end if

  end function calculate_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the Brownian diffusivity of a particle.
  !> Seinfeld and Pandis eq. 9.73.
  real(kind=dp) function compute_brownian_diff(diameter, temp)

    !> Particle diameter (m).
    real(kind=dp) :: diameter
    !> Temperature (K).
    real(kind=dp) :: temp

    real(kind=dp) :: C_c

    C_c = slip_correction_factor(diameter)

    compute_brownian_diff = (const%boltzmann * temp * C_c) / &
        (3.0d0 * const%pi * const%air_dyn_visc * diameter)


  end function compute_brownian_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module wrf_pmc_dep_aero
