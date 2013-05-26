
!void FORT_INIT_PHI(Real* data, const int *low, const int *high,
!        const int *nghost, const Real *dx,
!        const Real *problem_low, const Real *problem_high);
subroutine init_phi(phi, low, high, nghost, dx, prob_low, prob_high)
    implicit none

    !   param
    integer :: low(2), high(2), nghost
    double precision :: phi(low(1)-nghost : high(1)+nghost, low(2)-nghost : high(2)+nghost)
    double precision :: dx(2)
    double precision :: prob_low(2)
    double precision :: prob_high(2)

    !   local
    integer :: i, j
    double precision :: x, y, r2

    do j = low(2), high(2)
        y = prob_low(2) + (dble(j)+0.5) * dx(2)
        do i = low(1), high(1)
            x = prob_low(1) + (dble(i)+0.5) * dx(1)
            r2 = (x-0.25)**2 + (y-0.25)**2
            phi(i,j) = 1.0 + exp(-r2 * 100.0)
        enddo
    enddo


end subroutine init_phi


!void FORT_COMPUTE_FLUX(Real *phi, const int *nghostPhi,
!        Real *fluxx, Real *fluxy, const int *nghostFlux,
!        const int *low, const int *high, const Real *dx);
subroutine compute_flux(phi, ngp, fluxx, fluxy, ngf, lo, hi, dx)
    implicit none

    !   param
    integer :: lo(2), hi(2), ngp, ngf
    double precision :: phi(lo(1)-ngp : hi(1)+ngp, lo(2)-ngp : hi(2)+ngp)
    ! flux is staggered on its direction, so +1
    double precision :: fluxx(lo(1)-ngf : hi(1)+ngf+1, lo(2)-ngf : hi(2)+ngf)
    double precision :: fluxy(lo(1)-ngf : hi(1)+ngf, lo(2)-ngf : hi(2)+ngf+1)
    double precision :: dx(2)

    !   local
    integer :: i, j

    ! x-flux
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
            fluxx(i,j) = (phi(i,j) - phi(i-1,j)) / dx(1)
        enddo
    enddo

    ! y-flux
    do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
            fluxy(i,j) = (phi(i,j) - phi(i,j-1)) / dx(2)
        enddo
    enddo

end subroutine compute_flux



!void FORT_UPDATE_PHI(Real *phiOld, Real *phinew, const int *nghostPhi,
!        Real *fluxx, Real *fluxy, const int *nghostFlux,
!        const int *low, const int *high, const Real *dx, const Real *dt);
subroutine update_phi(phiold, phinew, ngp, fluxx, fluxy, ngf, lo, hi, dx, dt)
    implicit none

    !   param
    integer :: lo(2), hi(2), ngp, ngf
    double precision :: phiold(lo(1)-ngp : hi(1)+ngp, lo(2)-ngp : hi(2)+ngp)
    double precision :: phinew(lo(1)-ngp : hi(1)+ngp, lo(2)-ngp : hi(2)+ngp)
    ! flux is staggered
    double precision :: fluxx(lo(1)-ngf : hi(1)+ngf+1, lo(2)-ngf : hi(2)+ngf)
    double precision :: fluxy(lo(1)-ngf : hi(1)+ngf, lo(2)-ngf : hi(2)+ngf+1)
    double precision :: dx(2), dt

    !   local
    integer :: i, j
    double precision :: dpx, dpy

    ! update with flux divergence
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            dpx = (fluxx(i+1,j) - fluxx(i,j)) / dx(1)
            dpy = (fluxy(i,j+1) - fluxy(i,j)) / dx(2)

            phinew(i,j) = phiold(i,j) + dt * (dpx + dpy)
        enddo
    enddo

end subroutine update_phi







!BL_FORT_PROC_DECL(PHI_FILL, phi_fill) (
!        BL_FORT_FAB_ARG(state),
!        const int dlo[], const int dhi[],
!        const Real dx[], const Real glo[],
!        const Real *time, const int bc[]
!);
subroutine phi_fill(phi, phi_l1, phi_l2, phi_h1, phi_h2, &
    domlo, domhi, dx, xlo, time, bc)

    implicit none

    !   param
    integer :: phi_l1, phi_l2, phi_h1, phi_h2
    integer :: bc(2,2,*)
    integer :: domlo(2), domhi(2)
    double precision :: dx(2), xlo(2), time
    double precision :: phi(phi_l1:phi_h1, phi_l2:phi_h2)

    call filcc(phi, phi_l1, phi_l2, phi_h1, phi_h2, domlo, domhi, dx, xlo, bc)

end subroutine phi_fill


!BL_FORT_PROC_DECL(INIT_DATA, init_data) (
!        const int &level, const Real &time,
!        const int lo[], const int hi[],
!        const int &numState,
!        BL_FORT_FAB_ARG(state),
!        const Real dx[],
!        const Real xlow[], const Real xhi[]
!);
subroutine init_data(level, time, lo, hi, nstate, &
    state, state_l1, state_l2, state_h1, state_h2, &
    dx, xlo, xhi)

    implicit none
!
    integer :: level, nstate
    integer :: lo(2), hi(2)
    integer :: state_l1, state_l2, state_h1, state_h2
    double precision :: xlo(2), xhi(2), time, dx(2)
    double precision :: state(state_l1:state_h1, state_l2:state_h2, *)
!
    integer :: i, j
    double precision :: x, y, r2, val

    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            x = xlo(1) + (dble(i)-lo(1)+0.5) * dx(1)
            y = xlo(2) + (dble(j)-lo(2)+0.5) * dx(2)
            r2 = (x-0.25)**2 + (y-0.25)**2
            val = 1.0 + exp(-r2*100.0)

            state(i,j,1) = val
        enddo
    enddo

end subroutine init_data

!BL_FORT_PROC_DECL(CALC_FLUX, calc_flux) (
!        BL_FORT_FAB_ARG(old),
!        BL_FORT_FAB_ARG(fluxx), BL_FORT_FAB_ARG(fluxy),
!        const int lo[], const int hi[],
!        const Real dx[]
!);
subroutine calc_flux(state, state_l1, state_l2, state_h1, state_h2, &
    fluxx, fluxx_l1, fluxx_l2, fluxx_h1, fluxx_h2, &
    fluxy, fluxy_l1, fluxy_l2, fluxy_h1, fluxy_h2, &
    lo, hi, dx)

    implicit none

    ! param
    integer :: state_l1, state_l2, state_h1, state_h2
    double precision :: state(state_l1:state_h1, state_l2:state_h2)
    integer :: fluxx_l1, fluxx_l2, fluxx_h1, fluxx_h2
    double precision :: fluxx(fluxx_l1:fluxx_h1, fluxx_l2:fluxx_h2)
    integer :: fluxy_l1, fluxy_l2, fluxy_h1, fluxy_h2
    double precision :: fluxy(fluxy_l1:fluxy_h1, fluxy_l2:fluxy_h2)
    integer :: lo(2), hi(2)
    double precision :: dx(2)

    !   local
    integer :: i, j

    ! x-flux
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)+1
            fluxx(i,j) = (state(i,j) - state(i-1,j)) / dx(1)
        enddo
    enddo

    ! y-flux
    do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
            fluxy(i,j) = (state(i,j) - state(i,j-1)) / dx(2)
        enddo
    enddo

end subroutine calc_flux


!BL_FORT_PROC_DECL(UPDATE_STATE, update_state) (
!        BL_FORT_FAB_ARG(oldData), BL_FORT_FAB_ARG(newData),
!        BL_FORT_FAB_ARG(fluxx), BL_FORT_FAB_ARG(fluxy),
!        const int lo[], const int hi[],
!        const Real dx[], const Real *dt
!);
subroutine update_state(old, old_l1, old_l2, old_h1, old_h2, &
    new, new_l1, new_l2, new_h1, new_h2, &
    fluxx, fluxx_l1, fluxx_l2, fluxx_h1, fluxx_h2, &
    fluxy, fluxy_l1, fluxy_l2, fluxy_h1, fluxy_h2, &
    lo, hi, dx, dt)

    implicit none

    ! param
    integer :: old_l1, old_l2, old_h1, old_h2
    double precision :: old(old_l1:old_h1, old_l2:old_h2)
    integer :: new_l1, new_l2, new_h1, new_h2
    double precision :: new(new_l1:new_h1, new_l2:new_h2)
    integer :: fluxx_l1, fluxx_l2, fluxx_h1, fluxx_h2
    double precision :: fluxx(fluxx_l1:fluxx_h1, fluxx_l2:fluxx_h2)
    integer :: fluxy_l1, fluxy_l2, fluxy_h1, fluxy_h2
    double precision :: fluxy(fluxy_l1:fluxy_h1, fluxy_l2:fluxy_h2)
    !
    integer :: lo(2), hi(2)
    double precision :: dx(2), dt

    !   local
    integer :: i, j
    double precision :: dpx, dpy

    ! update with flux divergence
    do j = lo(2), hi(2)
        do i = lo(1), hi(1)
            dpx = (fluxx(i+1,j) - fluxx(i,j)) / dx(1)
            dpy = (fluxy(i,j+1) - fluxy(i,j)) / dx(2)

            new(i,j) = old(i,j) + dt * (dpx + dpy)
        enddo
    enddo

end subroutine update_state





!            FORT_PROBINIT(&init,
!                          probin_file_name.dataPtr(),
!                          &probin_file_length,
!                          Geometry::ProbLo(),
!                          Geometry::ProbHi());
subroutine probinit(init, name, namelen, problo, probhi)
    implicit none
!
    integer :: init, namelen
    integer :: name(namelen)
    double precision :: problo(2), probhi(2)

!
!    stop "PROBINIT"

end subroutine probinit


