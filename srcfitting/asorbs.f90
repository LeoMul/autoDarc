module     asorbs 
    !module for AS stuff.
    character*2 :: angular_symbols_neg(5) = (/'s ','p ','d ','f ','g ' /)

    !Integers:
    integer             :: num_orbs
    integer             :: num_points
    integer,allocatable :: princN(:)
    integer,allocatable :: orbL(:)

    !Reals:
    real*8,allocatable  :: orbitals(:,:)
    real*8,allocatable  :: potential(:,:)
    real*8,allocatable  :: orbitalParams(:)
    real*8,allocatable  :: orbitalEnergy(:)
    real*8,allocatable  :: orbitalFloat(:)
    real*8,allocatable  :: radial(:),charge(:)

    !parameters:
    real*8              :: nzed 
    integer             :: nelec

    contains 

        subroutine read_radout 

        !Reads the 'radout' (formatted) output file from autostructure.
        !The (kappa-averaged? not always) orbitals are stored in 
        !the array orbitals.
        !which has dimension num_orbs X num_points.
        !other quantities are stored.
        !dummy quantities such as point number are forgotten about.
        IMPLICIT none
        character*6,parameter :: atom = '      '
        character*6,parameter :: term = '      '
        !filename.
        character*6,parameter :: file ='radout'
        !dummy variables
        real*8                :: DA(5)
        character*4           :: tlbl       
        integer               :: j 
        integer               :: i
        integer               :: dummy 
        integer               :: iostat
        integer               :: key,mb
        write(*,10)
        open(1,file=file,form='formatted')

        read(1,1,iostat=iostat) key,num_orbs,mb,num_points,nelec,nzed,(da(i),i=1,5),tlbl

        
        call asorbs_alloc


        !read points. 
        do i = 1,num_points,2 
            read(1,2) key,dummy,radial(i),charge(i),dummy,radial(i+1),charge(i+1)
        end do 


        do i = 1,num_orbs
            !read header first
            read(1,3) key,dummy,princN(i),orbL(i),orbitalFloat(i) &
                    ,orbitalEnergy(i),dummy,tlbl
            !then read all the corresponding points. 
            do j = 1,num_points,2
                read(1,2) key,dummy,orbitals(i,j),potential(i,j),& 
                          dummy,orbitals(i,j+1),potential(i,j+1)
            end do 
            
            !extract near nucleus expansion parameter, set first orbital entry to zero.
            orbitalParams(i) = orbitals(i,1) 
            orbitals(i,1)  = 0.0d0 
            write(*,9) princN(i),angular_symbols_neg(orbL(i)+1),orbitalEnergy(i)* (-0.5d0)
            !lpm 26.03.25 
            call estimateAZ(int(i,8),orbitalParams(i))
            !print*,orbitalParams(i)

        end do 
        close(1)
        
        !correct the units of orbial energy to Hartree atomic units as is required by most codes. 
        !also AS gives the orbital eigenvalue as negative
        orbitalEnergy = orbitalEnergy * (-0.5d0)

        write(*,11) num_orbs,num_points

        !write(*,5)
        !header (10030 in autostructure.)
        1 format(3I5,4X,I9,I4,F5.1,5F7.3,A4)
        !entry for radial + charge arrays (and orbital arrays)
        2 format(I5,2(I4,2(1PE14.7)))
        !header for each orbital.
        3 format(3I5,I3,F5.1,2X,F12.6,I6,29X,A4)
!        4 format('Write out output:'
     9    format(1X,I2,A1,1X,'E= ',ES13.6)
    10    format(' Reading radout from AUTOSTRUCTURE')
    11    format(' Found ',I4,' non-relativistic orbitals with ',I5,' radial points.')
    end subroutine

    subroutine asorbs_alloc
        !allocates the asorbs
        allocate(radial(num_points))
        allocate(charge(num_points))
        allocate(orbitals(num_orbs,num_points))
        allocate(potential(num_orbs,num_points))
        allocate(orbitalParams(num_orbs))
        allocate(princN(num_orbs))
        allocate(orbL(num_orbs))
        allocate(orbitalEnergy(num_orbs))
        allocate(orbitalFloat(num_orbs))
    end subroutine

    subroutine asorbs_dlloc
        !deallocates all of the asorbs
        deallocate(radial)
        deallocate(charge)
        deallocate(orbitals)
        deallocate(potential)
        deallocate(orbitalParams)
        deallocate(princN)
        deallocate(orbL)
        deallocate(orbitalEnergy)
        deallocate(orbitalFloat)
    end subroutine

    subroutine estimateAZ(orbital_index,az)

        !calculate the best az to use
        !i.e - the best value of the derivative of the wavefunction at zero.
        !only makes sense to use this for s orbitals, really 
        !as strictly (in a non-relativistic setting) the derivative is zero 
        !for l>= 1 at r=0. 

        real*8,intent(out) :: az
        real*8             :: int1,int2 
        real*8             :: rmax 
        real*8,allocatable ::table0(:), table1(:),table2(:)
        integer            :: ii,cutoff 
        integer*8          :: orbital_index
        integer            :: l 

        !intialize
        az = 0.0d0 
        rmax = 1.e-3 
        l = orbL(orbital_index)

        cutoff = minloc(abs(radial - rmax ),1)
        cutoff = 12
        allocate(table0(cutoff))
        allocate(table1(cutoff))
        allocate(table2(cutoff))

        int1 = 0.d0 
        int2 = 0.d0

        table0 = radial**(l+1) * (1.0d0 - nzed * radial)

        do ii = 1,cutoff
            table1(ii) = table0(ii) * orbitals(orbital_index,ii)
        end do
        table2 = table0 * table0

        do ii = 1,cutoff-1
            int1 = int1 + (table1(ii+1)+table1(ii)) * (radial(ii+1)-radial(ii))
            int2 = int2 + (table2(ii+1)+table2(ii)) * (radial(ii+1)-radial(ii))
        end do 

        az = int1/int2 
        !cleanup
        deallocate(table0)
        deallocate(table1)
        deallocate(table2)
    end subroutine
    
end module asorbs