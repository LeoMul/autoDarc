module trapz 
    contains
    function trapz_int(rgrid,func)

        use kind_param 
        implicit none 
        integer :: n
        real(rkind),dimension(:) :: rgrid , func 
        integer(ikind) :: ii
        real(rkind)    :: trapz_int 
        real(rkind)    :: h 
        real(rkind)    :: ff 

        integer(ikind) :: size1,size2,maxp
        size1 = size(rgrid)
        size2 = size(func)
        maxp = min(size1,size2)
        trapz_int = 0.0d0 
        !print*,maxp
        if (maxp .gt. 2) then 
            do ii = 2,maxp 
                h  = rgrid(ii) - rgrid(ii-1) 
                ff =  func(ii) +  func(ii-1)
                trapz_int = trapz_int + 0.5d0 * h * ff 
            end do 
        end if 

    end function trapz_int
end module trapz 