! Copyright (C) 2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc_dep module.

!> Main driver for partmc deposition in WRF.
module wrf_pmc_dep

  use wrf_pmc_dep_aero
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_state
  use pmc_env_state

  use module_wrf_error
  use module_domain_type
  use module_domain
  use module_model_constants
  use module_configure, only: p_qv
  USE module_state_description
  USE module_domain, ONLY : domain, domain_clock_get
  USE module_configure, ONLY : grid_config_rec_type

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do dry deposition.
  subroutine wrf_pmc_dep_dry_driver(grid, config_flags, aero_data, &
         aero_states, gas_states, env_states, pmc_is, pmc_ie, pmc_js, pmc_je, &
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
    !>
    type(grid_config_rec_type) :: config_flags
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data 
    !> Aerosol states.
    type(aero_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: aero_states
    !> Gas states.
    type(gas_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(inout) :: gas_states
    !> Environment states.
    type(env_state_t), dimension(pmc_is:pmc_ie,pmc_ks:pmc_ke,pmc_js:pmc_je), &
         intent(in) :: env_states
    integer, intent(in) :: global_nx
    integer, intent(in) :: global_ny
    integer, intent(in) :: global_nz

    ! Local variables
    integer :: nx, ny, nz
    real(kind=dp) :: dtstep
    real(kind=dp), dimension(pmc_is:pmc_ie,pmc_js:pmc_je) :: zmid
    real(kind=dp), dimension(pmc_is:pmc_ie,pmc_js:pmc_je) :: aer_res, ust
    integer :: i,j,k
    real(kind=dp) :: zr

    zmid = 0.0d0

    do i = pmc_is, pmc_ie
    do j = pmc_js, pmc_je
       zmid(i,j)= grid%z_at_w(i,2,j) - grid%z_at_w(i,1,j)
       ust(i,j) = real(grid%ust(i,j),dp)
    end do
    end do

    dtstep = real(grid%dt,dp)

    zr = 2.0
    do i = pmc_is, pmc_ie
    do j = pmc_js, pmc_je
       call depvel_partmc(grid%rmol(i,j), zr, real(grid%znt(i,j),kind=dp), &
            ust(i,j), aer_res(i,j))
    end do
    end do

    call dep_aero_dry_driver(grid, aero_states, aero_data, env_states, &
         ust, global_nx, global_ny, dtstep, aer_res, pmc_is, pmc_ie, &
         pmc_js, pmc_je, pmc_ks, pmc_ke)

  end subroutine wrf_pmc_dep_dry_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out deposition velocities for gas species.
  subroutine print_chem(ddvel)

    !> Dry deposition velocities.
    real, dimension(:) :: ddvel

    print*, p_so2, "so2", ddvel(p_so2)
    print*, p_sulf, "sulf", ddvel(p_sulf)
    print*, p_no2, "no2", ddvel(p_no2)
    print*, p_no, "no", ddvel(p_no)
    print*, p_o3, "o3", ddvel(p_o3)
    print*, p_hno3, "hno3", ddvel(p_hno3)
    print*, p_h2o2, "h2o2", ddvel(p_h2o2)
    print*, p_ald, "ald", ddvel(p_ald)
    print*, p_hcho, "hcho", ddvel(p_hcho)
    print*, p_op1, "op1", ddvel(p_op1)
    print*, p_op2, "op2", ddvel(p_op2)
    print*, p_paa, "paa", ddvel(p_paa)
    print*, p_ora1, "ora1", ddvel(p_ora1)
    print*, p_ora2, "ora2", ddvel(p_ora2)
    print*, p_nh3, "nh3", ddvel(p_nh3)
    print*, p_n2o5, "n2o5", ddvel(p_n2o5)
    print*, p_no3, "no3", ddvel(p_no3)
    print*, p_pan, "pan", ddvel(p_pan)
    print*, p_hc3, "hc3", ddvel(p_hc3)
    print*, p_hc5, "hc5", ddvel(p_hc5)
    print*, p_hc8, "hc8", ddvel(p_hc8)
    print*, p_eth, "eth", ddvel(p_eth)
    print*, p_ch4, "ch4", ddvel(p_ch4)
    print*, p_co, "co", ddvel(p_co)
    print*, p_ol2, "ol2", ddvel(p_ol2)
    print*, p_olt, "olt", ddvel(p_olt)
    print*, p_oli, "oli", ddvel(p_oli)
    print*, p_tol, "tol", ddvel(p_tol)
    print*, p_xyl, "xyl", ddvel(p_xyl)
    print*, p_aco3, "aco3", ddvel(p_aco3)
    print*, p_tpan, "tpan", ddvel(p_tpan)
    print*, p_hono, "hono", ddvel(p_hono)
    print*, p_hno4, "hno4", ddvel(p_hno4)
    print*, p_ket, "ket", ddvel(p_ket)
    print*, p_gly, "gly", ddvel(p_gly)
    print*, p_mgly, "mgly", ddvel(p_mgly)
    print*, p_dcb, "dcb", ddvel(p_dcb)
    print*, p_onit, "onit", ddvel(p_onit)
    print*, p_csl, "csl", ddvel(p_csl)
    print*, p_iso, "iso", ddvel(p_iso)
    print*, p_ho, "ho", ddvel(p_ho)
    print*, p_ho2, "ho2", ddvel(p_ho2)
    print*, p_hcl, "hcl", ddvel(p_hcl)
    print*, p_ch3o2, "ch3o2", ddvel(p_ch3o2)
    print*, p_ethp, "ethp", ddvel(p_ethp)
    print*, p_ch3oh, "ch3oh", ddvel(p_ch3oh)
    print*, p_c2h5oh, "c2h5oh", ddvel(p_c2h5oh)
    print*, p_par, "par", ddvel(p_par)
    print*, p_to2, "to2", ddvel(p_to2)
    print*, p_cro, "cro", ddvel(p_cro)
    print*, p_open, "open", ddvel(p_open)
    print*, p_op3, "op3", ddvel(p_op3)
    print*, p_c2o3, "c2o3", ddvel(p_c2o3)
    print*, p_ro2, "ro2", ddvel(p_ro2)
    print*, p_ano2, "ano2", ddvel(p_ano2)
    print*, p_nap, "nap", ddvel(p_nap)
    print*, p_xo2, "xo2", ddvel(p_xo2)
    print*, p_xpar, "xpar", ddvel(p_xpar)
    print*, p_isoprd, "isoprd", ddvel(p_isoprd)
    print*, p_isopp, "isopp", ddvel(p_isopp)
    print*, p_isopn, "isopn", ddvel(p_isopn)
    print*, p_isopo2, "isopo2", ddvel(p_isopo2)

    ! Typically not included species
    ! These are marine gas species that PartMC-MOSAIC includes
    !print*,p_dms,"dms"
    !print*,p_msa,"msa"
    !print*,p_dmso,"dmso"
    !print*,p_dmso2,"dmso2"
    !print*,p_ch3so2h,"ch3so2h"
    !print*,p_ch3sch2oo,"ch3sch2oo"
    !print*,p_ch3so2,"ch3so2"
    !print*,p_ch3so3,"ch3so3"
    !print*,p_ch3so2oo,"ch3so2oo"
    !print*,p_ch3so2ch2oo,"ch3so2ch2oo"
    !print*,p_mtf,"mtf"

  end subroutine print_chem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine depvel_partmc(rmol, zr, z0, ustar, aer_res)

    real(kind=dp), intent(in)    :: ustar, z0
    real(kind=dp) :: zr
    real(kind=dp), intent(out)   :: aer_res
    real, intent(inout) :: rmol

    real(kind=dp) :: ao, ar, polint, vk

    vk = 0.40d0

    if(abs(rmol) < 1.0e-6) rmol = 0.

    IF (rmol < 0.0) THEN
      ar = ((1.0-9.0*zr*rmol)**(0.25)+0.001)**2
      ao = ((1.0-9.0*z0*rmol)**(0.25)+0.001)**2
      polint = 0.74*(dlog((ar-1.0)/(ar+1.0))-dlog((ao-1.0)/(ao+1.0)))
    ELSE IF (rmol == 0.0) THEN
      polint = 0.74*dlog(zr/z0)
    ELSE
      polint = 0.74*dlog(zr/z0) + 4.7*rmol*(zr-z0)
    END IF

    aer_res = polint/(vk*max(ustar,1.0e-4))

  end subroutine depvel_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module wrf_pmc_dep
