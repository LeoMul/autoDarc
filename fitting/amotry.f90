FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
    implicit none
    INTEGER ihi,mp,ndim,np,NMAX
    real*8 amotry,fac,p(mp,np),psum(np),y(mp),funk
    PARAMETER (NMAX=20)
    EXTERNAL funk
    ! USES funk
    !Extrapolates by a factor fac through the face of the simplex across from the high point,
    !tries it, and replaces the high point if the new point is better.
    INTEGER j
    real*8 fac1,fac2,ytry,ptry(NMAX)
    fac1=(1.-fac)/ndim
    fac2=fac1-fac
    do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
    11 enddo
    ytry=funk(ptry) 
    if (ytry.lt.y(ihi)) then 
        y(ihi)=ytry
        do 12 j=1,ndim
            psum(j)=psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j)=ptry(j)
        12 enddo 
    endif
    amotry=ytry
    return
END