program main
    use kind_param
    use factorials
    use trapz
    use data_for_funk
    use minimization
    use asorbs
    implicit none  
    integer :: isize
    real(rkind),allocatable :: alpha(:)
    real(rkind) :: chisquared

    integer :: l ,index , ii

    real(rkind),allocatable :: rtest(:),ftest(:)

    integer :: npoints
    real*8 :: tt 

    call init_factorials


    !open(1,file='inptest')
    !read(1,*) npoints
    !allocate(rtest(npoints),ftest(npoints))
    !do  ii = 1,npoints
    !    read(1,*) rtest(ii), ftest(ii)
    !end do
    !close(1)
    call read_radout
    index = 2
    isize = princN(index) + 2
    print*,isize

    allocate(lv(isize))
    do ii = 1,isize 
        lv(ii) = myl + ii - 1
    end do
    lv(:) = myl
    allocate(alpha(isize))

    npoints = num_points
    myrad  = radial
    myorb  = orbitals(index,:)
    print*, myorb(1:4)

    alpha(1)  = nzed 
    !alpha(2:) = 5.0d0 

    call random_number(alpha)
    alpha = alpha * 10.0d0
    alpha(1)  = nzed 

    myl = orbL(index)
    allocate(slaterorb(npoints))
    rlp1 = myrad ** (myl + 1)
    chisquared = chisq(alpha)
    print*,myl

    call flush 
    print*,chisquared 

    call minimize(alpha,isize,chisq)
    
    print*,alpha
    do ii = 1,npoints 
        write(100,*) myrad(ii),myorb(ii) ,slaterorb(ii)
    end do 

    stop 

    mynpoints = npoints



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

    !call minimize(alpha,isize,chisq)

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

        call linlsq(alpha,c,isize,myl,myorb,myrad,npoints)
        !print*,c
        call construct_slater_fit(alpha,c,isize)
        chisq = 0.0 
        do ii = 1,npoints 
            dd = myorb(ii) - slaterorb(ii)
            chisq = chisq + dd*dd
        end do 


    end function

end program 