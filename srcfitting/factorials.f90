module factorials 
    use kind_param
    integer,parameter :: max_factorial = 100 

    real(rkind) :: factorials_stored(0:max_factorial)

    contains 

    subroutine init_factorials
        implicit none  
        integer :: ii 
        factorials_stored = 1.0  

        do ii = 1,max_factorial
            factorials_stored(ii) = factorials_stored(ii-1) * real(ii,rkind)
        end do 

    end subroutine

end module 