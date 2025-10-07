module data_for_funk 
    use kind_param
    real(rkind), allocatable :: myorb(:)
    real(rkind), allocatable :: myrad(:)
    real(rkind), allocatable :: slaterorb(:)
    real(rkind), allocatable :: rlp1(:)
    integer(ikind)           :: mynpoints
    integer(ikind)           :: myl 
    integer(ikind)  , allocatable         :: lv(:) !vector of l's in expansion
end module 

subroutine make_a_matrix(amatrix,alphaarray,isize,l)
    use kind_param
    use factorials
    use data_for_funk
    implicit none
    !inputs
    integer(ikind),intent(in   ) :: isize  
    integer(ikind),intent(in   ) :: l  

    real   (rkind),intent(inout) :: amatrix(isize,isize)
    real   (rkind),intent(inout) :: alphaarray(isize)

    !processing
    integer(ikind) :: ii,jj
    real(rkind) :: den 
    real(rkind) :: a1,a2,n1,n2,fac 
    !external
    real(rkind),external ::norm_slater
    !print*,'in make a matrix'
    fac = factorials_stored(2*l+2)

    !print*,alphaarray(1)
    !call flush
    do ii = 1,isize 
        a1 = alphaarray(ii)
        n1 = norm_slater(lv(ii),a1)
        !print*,n1
        do jj = 1,isize 
            a2 = alphaarray(jj)
            n2 = norm_slater(lv(jj),a2)
            den = a1 + a2 
            den = den** (lv(ii) + lv(jj) + 3)
            amatrix(ii,jj) = n1 * n2 * fac / den
        end do 
    end do 


end subroutine

subroutine make_b_vector(bvector,alphaarray,isize,l,orbital,radial,npoints)
    use kind_param
    use factorials
    use trapz
    use data_for_funk

    implicit none
    !inputs
    integer(ikind),intent(in   ) :: isize  ,npoints
    integer(ikind),intent(in   ) :: l  

    real   (rkind),intent(inout) :: bvector(isize)
    real   (rkind),intent(inout) :: alphaarray(isize)
    real   (rkind),intent(in ) :: radial (npoints)
    real   (rkind),intent(in ) :: orbital(npoints)

    !processing
    integer(ikind) :: ii,jj
    real(rkind) :: den 
    real(rkind) :: a1,a2,n1,n2,fac
    
    real(rkind),allocatable :: trialorb(:)
    real(rkind),allocatable :: rl(:)
    integer(ikind) :: sizearray
    real(rkind),external ::norm_slater

    sizearray = size(radial)
    allocate(trialorb(sizearray))


    do ii = 1,isize 
        rl = radial ** (lv(ii) + 1)
        a1 = alphaarray(ii)
        trialorb = rl * norm_slater(lv(ii),a1) * exp(-a1 * radial)

        bvector(ii) = trapz_int(radial,orbital*trialorb)

    end do 
    print*,'---------------------------------------------'
    print*, ' in b vector'
    print*,'a= ',alphaarray
    print*,'b= ',bvector
    print*,'---------------------------------------------'

    !do ii = 1,size(radial)
    !    write(100,*) radial(ii),trialorb(ii)
    !end do 



end subroutine

subroutine get_c_vector(cvector,amatrix,bvector,ndim)
    use kind_param
    use factorials
        use data_for_funk

    implicit none 
    integer(ikind) :: ndim
    real(rkind) , intent(inout) :: cvector(ndim) 
    real(rkind) :: bvector(ndim) 
    real(rkind) :: amatrix(ndim,ndim),a(ndim,ndim)
    integer(ikind) :: ipiv(ndim)
    integer :: info
    real(rkind) :: cnorm,c1,c2
    integer(ikind) :: ii ,jj 
    cvector = bvector
    a = amatrix
    print*,'---------------------------------------------'

    call DGETRF(ndim,ndim,A,ndim,IPIV,INFO)
    call DGETRS('N', ndim, 1, amatrix, ndim, ipiv, cvector, &
               ndim, info)
    print*,'ndim = ',ndim 
    print*,'info = ',info
    print*,'raw c = ',cvector
    cnorm = 0.0d0 
    do ii = 1,ndim 
        c1 = cvector(ii)
        do jj = 1,ndim 
            c2 = cvector(jj)
            cnorm = cnorm + c1 * c2 * A(ii,jj)
            print*,cnorm

        end do 
    end do 
    cnorm = sqrt(cnorm)
    cvector = cvector / cnorm
    print*,'In c vector'
    print*,'A=',amatrix
    print*,'b=',bvector
    print*,'corm=',cnorm
    print*,'c=',cvector
    print*,'---------------------------------------------'

end subroutine

subroutine linlsq(alpha,c,isize,l,orbital,radial,np)
    use kind_param
    use factorials
    implicit none 
    integer :: isize,l
    real(rkind),intent(inout) :: alpha(isize)
    real(rkind) :: a(isize,isize)
    real(rkind) :: b(isize)
    real(rkind),intent(inout) :: c(isize)
    real(rkind),intent(in) :: orbital(np),radial(np)
    integer(ikind) :: np 
    call make_b_vector(b,alpha,isize,l,orbital,radial,np)

    call make_a_matrix(a,alpha,isize,l)
    !print*,'going into c'
    call get_c_vector(c,a,b,isize)

end subroutine

subroutine construct_slater_fit(alpha,c,isize)
    use kind_param 
    use data_for_funk 
    implicit none 

    real*8 :: alpha(isize)
    real*8 :: c(isize)
    real*8 :: a 
    integer :: isize  
    integer :: ii 
    real(rkind),external :: norm_slater
    slaterorb = 0.0d0 
    
    do ii = 1, isize 
        a = alpha(ii)
        slaterorb =  slaterorb + myrad**(lv(ii)+1) * &
                     c(ii) * norm_slater(lv(ii),a) * exp( - a * myrad)
    end do 
    print*,'slater sum = ', sum(slaterorb)
    !slaterorb = slaterorb * rlp1

end subroutine
