module minimization 
    contains
    subroutine minimize(startingpoint,ndim,funk)
        !minimizes a function.
        integer :: ndim 
        real*8  :: startingpoint(ndim)
        real*8 :: length = 0.001d0 
        real*8 :: xx(ndim)
        real*8 :: x0(ndim)
        real*8 :: yy(ndim+1)
        real*8 :: p(ndim+1,ndim)
        real*8 :: ftol = 1e-7 
        real*8 :: funk
        integer :: iter = 0 
        integer :: np  
        integer :: mp 
        integer :: ii 
        real*8 :: pav(ndim)
        np = ndim 
        mp = ndim + 1
        xx = startingpoint

        ii=1
        x0 = xx 
        p(ii,:) = x0 
        !print*,p(ii,:)
        yy(ii) = funk(x0)
            
        do ii = 2,ndim+1 
            x0 = xx 
            x0(ii-1) = xx(ii-1) + length
            p(ii,:) = x0 
            !print*,p(ii,:)
            yy(ii) = funk(x0)
        end do 

        !print*,yy(:)
        call amoeba(p,yy,mp,np,ndim,ftol,funk,iter)
        pav = 0.0d0 
        do ii = 1,ndim+1 
            pav(:) = pav(:) + p(ii,:)
        end do 
        pav = pav / (ndim+1)
        startingpoint = pav 
    end subroutine

end module 