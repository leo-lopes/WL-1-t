module ising_wl_rot
    implicit none
    
    integer                 ::  L, L2, L2m1, energ, lambda, lam2, istate, jstate, dE, Emax
    integer, parameter      ::  mc_steps = 1000
    integer, allocatable    ::  k(:,:), E(:), H(:), inm(:), inp(:), indE(:)

    real(8)                 ::  t, avgE, avgE2, Cv, Z,  iL2, iL2m1, f, err
    real(8), parameter      ::  f0 = 1.d0, fstop = 10.d-7, hist_pr=0.8d0, hf=1.d0/2
    real(8), allocatable    ::  S(:), g(:)

contains

    subroutine alloc()
        allocate(k(L,L), inm(L), inp(L), E(L2m1), H(L2m1), g(L2m1), S(L2m1))
    end subroutine

    subroutine dealloc()
        deallocate(k, E, H, g, S)
    end subroutine

    subroutine initialize(ic)
        character(len=4), intent(in) :: ic
        integer    ::   i, j, clock
        real       ::   r
        call alloc()

        g=1.d0
        H=0
        S=0.d0
        k=1
        t=0.d0

        do i = 1, L
            inm(i) = i-1
            inp(i) = i+1
        enddo
        inm(1) = L
        inp(L) = 1

        if (ic=='rand') then
            do i = 1,L
                do j = 1,L
                    call system_clock(count=clock)
                    call random_seed(clock)
                    call random_number(r)
                    if (r<= 0.5) then
                        k(i,j)=-1
                    endif
                enddo
            enddo
        else if (ic=="allp") then
            k=1
        else if (ic=="allm") then
            k=-1
        endif

        E(1)=-2*L2
        E(L2m1)=-E(1)
        do i = 2,L2m1-1
            E(i) = 4*i - 2*L2
        enddo   

        call calculate_energy()

        Emax = E(L2m1)

        allocate(indE(2*Emax + 1))
        indE = -1

        do i = 1, L2m1
            indE(E(i)+Emax + 1)=i
        enddo

    end subroutine

    subroutine calculate_energy()
        integer     ::  i, j, ip, jp, im, jm, kij

        energ = 0 
        do i= 1,L
            do j=1,L
                kij = k(i,j)
                ip = inp(i)
                im = inm(i)
                jp = inp(j)
                jm = inm(j)

                energ = energ -kij*(k(ip,j)+k(i,jp)+k(im,j)+k(i,jm))
            enddo
        enddo
        energ = energ/2
    end subroutine

    subroutine calculate_spin_energy(i,j)
        integer, intent(in) :: i,j
        integer     ::  ip, jp, im, jm, kij
        dE = 0

        kij = k(i,j)
        ip = inp(i)
        im = inm(i)
        jp = inp(j)
        jm = inm(j)

        dE = 2*kij*(k(ip,j)+k(i,jp)+k(im,j)+k(i,jm))
    end subroutine

    subroutine test_flatness(result)
        logical, intent(out) :: result
        integer ::  i
        real(8) ::  avgH, inv_avgH, maxH, minH

        avgH=0
        do i = 1, L2m1
            avgH = avgH + H(i)
        enddo
        avgH=avgH/L2m1
        inv_avgH=1.d0/avgH

        maxH= maxval(H)
        minH= minval(H)

        ! print *, avgH*hist_pr, maxH, minH
        if (avgH*hist_pr<minH) then
            result = .true.
        else
            result = .false.
        endif 

    end subroutine 

    subroutine wl_mc_step()
        integer     ::     state, new_state, clock, new_E
        real(8)     ::     p, dS, a, r(2)

        ! select initial state
        call system_clock(count=clock)
        call random_seed(clock)
        call random_number(r)
        istate = floor(r(1)*L) + 1
        jstate = floor(r(2)*L) + 1
        call calculate_spin_energy(istate,jstate)

        ! print *, istate, jstate

        new_E = energ + dE ! new energy if you flip the spin

        state = indE(energ+Emax+1) ! finding the index of the energy
        new_state = indE(new_E+Emax+1) ! finding the index of the new energy

        dS = S(new_state) - S(state) ! calculate difference in entropy

        if (1.d0 < dexp(-dS)) then
            p = 1.d0
        else if (dexp(-dS)<=1.d0) then
            p = dexp(-dS)
        endif
        ! p= dexp(-dS)

        call system_clock(count=clock)
        call random_seed(clock)
        call random_number(a)
        
        if (a<p) then
            energ = new_E
            k(istate,jstate)=-k(istate,jstate)
            state=new_state
        endif

        S(state) = S(state) + f
        H(state) = H(state) + 1
        t = t + iL2
    end subroutine

    subroutine normalize()
        integer             :: i
        real(8)             :: S_temp

        if (S(L2m1)<S(1)) then
            S_temp = S(1) + dlog(1 + dexp(S(L2m1)-S(1))) - dlog(4.d0)
        else
            S_temp = S(L2m1) + dlog(1 + dexp(S(1)-S(L2m1))) - dlog(4.d0)
        endif
        S = S - S_temp
        do i = 1, L2m1
            if (S(i)<0) then 
            S(i)=0
            endif
        enddo
    end subroutine 

    subroutine error()
        integer     ::  i
        real(8)     ::  avg_S, avg_S2

        avg_S=0.d0
        avg_S2=0.d0
        err = 0.d0
        do i = 1, L2m1
            avg_S = avg_S + S(i)
            avg_S2 = avg_S2 + S(i)*S(i)
        enddo
        avg_S = avg_S*iL2m1
        avg_S2 = avg_S*iL2m1

        do i = 1,L2m1
            err = err + (S(i)-avg_S)**2
        enddo

        err = err*iL2m1
        err = dsqrt(err)
    end subroutine

    subroutine resampling(beta)
        integer              ::  i
        real(8)              ::  w
        real(8), intent(in)  :: beta

        avgE=0.d0
        avgE2=0.d0
        Cv=0.d0
        Z=0.d0
        w=0.d0
        do i=1,L2m1
            w = dexp(S(i)-S(1)-beta*(E(i)+Emax))
            avgE = avgE + w*E(i)
            avgE2 = avgE2 + w*E(i)*E(i)
            Z = Z + w
        enddo
        avgE = avgE/Z
        avgE2 = avgE2/Z
        Cv = (avgE2 - avgE**2)*beta**2
        
        avgE = avgE*iL2
        Cv   = Cv*iL2
    end subroutine

end module ising_wl_rot

program ising_wl
    use ising_wl_rot
    implicit none

    real(8), parameter    ::    dT = 1.d-4, Tf = 5.d0
    integer(8), parameter ::    Tsize = floor((Tf-0.5d0)/dT)

    character(len=5)      ::    ic
    character(len=2)      ::    method, Ls
    integer               ::    i, MC_total, u_step, iu_step
    real(8)               ::    Temp, beta, beta_ar(Tsize), invt, time
    logical               ::    is_flat

    L=10
    Temp = 0.5
    beta_ar(1) = 1.d0/Temp

    do i = 2, Tsize
        Temp = Temp + dT
        beta_ar(i)= 1.d0/Temp
    enddo

    L2=L*L
    L2m1 = L2 - 1
    iL2 = 1.d0/L2
    iL2m1 = 1.d0/L2m1

    MC_total = L2*mc_steps
    u_step = 100
    iu_step = u_step*L2
    ic="rand"
    method="wl" ! or it

    write (Ls , '(I0)') L

    call initialize(ic)
    
    f=f0
    print *, "initialization done"
    i=0
    time=0.d0
    
    open(unit=3,file='./data/L'//Ls//'/'//method//'_err_vs_t.dat',status='unknown',action='write')    
    if (method=="wl") then
        do while (f>fstop)
            do while (t<= mc_steps)
                call wl_mc_step()
                i=i+1
                time=time+t
                if (mod(i,iu_step) == 0) then
                    call error()
                    write(3,*) time, err, f
                endif
            enddo
            t=0
            call test_flatness(is_flat)
            if (is_flat) then
                H=0
                f = f*hf
                print *, "is flat, f= ", f
            endif
        enddo
    else if (method=="it") then 
1       call wl_mc_step()
        i=i+1
        if (mod(i,MC_total) == 0) then
            call test_flatness(is_flat)
            if (is_flat) then
                H=0
                f = f*hf
                print *, "is flat, f= ", f
                invt = 1.d0/t
            endif
        endif
        if (mod(i,iu_step) == 0) then
            call error()
            write(3,*) t, err, f
        endif
        if (f>invt) then 
            go to 1
        else
            print *, "f -> 1/t"
            i=0
            f = invt
            do while (f>fstop)
                call wl_mc_step()
                i=i+1
                f = 1.d0/t
                if (mod(i,1000000)==0) then
                    print *, "f is now", f
                endif
            enddo
            if (mod(i,iu_step) == 0) then
                call error()
                write(3,*) t, err, f
            endif
        endif
    endif

        open(unit=1,file='./data/L'//Ls//'/'//method//'_avgE_Cv.dat',status='unknown',action='write') ! save results
        open(unit=2,file='./data/L'//Ls//'/'//method//'_S.dat',status='unknown',action='write')
        Temp=0.d0 
        print *, "Saving"
        call normalize()
        g = dexp(S)

        do i = 1,L2m1
            write (2,*) S(i), g(i),E(i)
        enddo
        do i = 1, Tsize
            beta = beta_ar(i)
            Temp = 1.d0/beta
            call resampling(beta)
            write (1,*) avgE, Cv, Temp
        enddo


    close(1)
    close(2)
    close(3)

end program ising_wl
