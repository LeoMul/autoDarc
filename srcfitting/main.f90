program main
    use kind_param
    use factorials
    use trapz
    use data_for_funk
    use minimization
    implicit none  
    integer,parameter :: isize = 2
    real(rkind) :: a(isize,isize)
    real(rkind) :: b(isize)
    real(rkind) :: c(isize)
    real(rkind) :: alpha(isize)
    real(rkind) :: chisquared
    integer :: ipiv(isize)
    integer :: info 
    integer :: ii 
    integer :: l 

    real(rkind),allocatable :: rtest(:),ftest(:)
    integer :: npoints
    real*8 :: tt 
    open(1,file='inptest')
    read(1,*) npoints
    allocate(rtest(npoints),ftest(npoints))
    do  ii = 1,npoints
        read(1,*) rtest(ii), ftest(ii)
    end do
    close(1)
    mynpoints = npoints

    call init_factorials


    alpha(:) = 0.5d0
    alpha(2) = 5.0d0  
    l = 0 

    myl = l 
    myrad = rtest 
    myorb = ftest
    allocate(slaterorb(npoints))
    rlp1 = myrad ** (myl + 1)
    !call linlsq(alpha,c,isize,l,ftest,rtest,npoints)

    chisquared = chisq(alpha)

    call minimize(alpha,isize,chisq)

    print*,alpha

    contains 

    function chisq(alpha)
        use data_for_funk
        real*8 :: chisq 
        real*8 :: alpha(isize)
        real(rkind) :: c(isize)
        integer :: npoints , ii 
        real(rkind) :: dd 
        npoints = size(myrad)

        call linlsq(alpha,c,isize,l,myorb,myrad,npoints)
        !print*,c
        call construct_slater_fit(alpha,c,isize)
        chisq = 0.0 
        do ii = 1,npoints 
            dd = myorb(ii) - slaterorb(ii)
            chisq = chisq + dd*dd
        end do 


    end function

end program 