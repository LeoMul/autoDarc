function norm_slater(l,a)
    !calculates norm of a slater orbital 
    use factorials
    use kind_param
    real(rkind)    :: norm_slater  
    real(rkind)    :: a 
    integer(ikind) :: l 
    norm_slater = (2.0*a) ** (2*l+3)
    norm_slater = norm_slater / factorials_stored(2*l+2)
    norm_slater = sqrt(norm_slater)
    print*,'ns = ',norm_slater
end function