!Program to parse autostructure data to the 
!darc R-matrix codes. 

!The code right now assumes an IC(R) input, producing a list of Dirac-orbitals.
!The large components are taken as the AS orbital: 
!   P_nk = P_nl .
!The small components are taken to be: 
!   Q_nk = 0.5 * \alpha * (d/dr + \kappa/r) * P_nl + O(\alpha^2)

!(Aside:
!The above expansion is is a truncation of the full: 

!   Q_nk = 0.5 * \alpha * (d/dr + \kappa/r) * P_nl / [1+ 0.25 alpha^2 * (e-V)]
!The correction only really matters close to the nucleus - where the potential
!is close to being a purely hydrogenic -Z/r where Z is the full nuclear charge.
!)

!which are then normalised according to:

!\int dr (PP* + QQ*) = 1.

!Notes:
!Using these non-relativistic orbitals from optimizing a 
!Breit-Pauli-Schrodinger Hamiltonian in a Dirac Hamiltonian
!is equivalent to the so-called RCI-P method, as breifly introduced in 
!Charlotte Froese Fischer et al 2016 J. Phys. B: At. Mol. Opt. Phys. 49 182004
!DOI 10.1088/0953-4075/49/18/182004

!This provides a way to use the much more efficient jj-coupled darc r-matrix codes.
!Instead of first doing an LS run and doing an expensive re-couple using stgjk. 

!To do: 
!Have each of the routines called optionally.
!Replace the optional input routine with command line options.
!I am of the opinion that utility codes such as these should be
!simply ran from the command line instead of worying about inputs 
!and optimizing runs.
!The above points go hand in hand.
!https://cyber.dabamos.de/programming/modernfortran/command-line-arguments.html


program autodarc 
    IMPLICIT none

    !fine structure constant
    real*8,parameter :: fineStruc = 0.0072973525643d0 

    !global Variables 
    real*8,allocatable :: orbitals(:,:)
    real*8,allocatable :: orbitalsLarge(:,:)
    real*8,allocatable :: orbitalsSmall(:,:)
    real*8,allocatable :: darcRadial(:)
    character*2,allocatable:: angular_string(:)

    real*8,allocatable :: potential(:,:)
    real*8,allocatable :: orbitalParams(:)
    integer,allocatable :: princN(:)
    integer,allocatable :: orbL(:)
    real*8,allocatable :: orbitalEnergy(:)
    real*8,allocatable :: orbitalFloat(:)
    real*8,allocatable :: radial(:),charge(:)
    logical :: inputExists
    logical :: callOrthog 
    
    real*8,allocatable :: genOcc(:)

    integer :: key,num_orbs,mb,nelec
    integer*4 :: num_points
    real*8 :: nzed 
    character*2 :: angular_symbols_neg(5) = (/'s ','p ','d ','f ','g ' /)
    character*2 :: angular_symbols_pos(4) =      (/'p-','d-','f-','g-' /)
    integer,allocatable :: kappaArray(:),princNREL(:),orbMap(:)
    integer :: numRelOrbs 

    real*8 :: rnt
    real*8 :: hff 

    call optionalInput

    !Call the program.
    write(*,1)
    call read_radout
    write(*,1)

    !I cant remember why i coded this with input variables. 
    call findROrbitals(orbL,princN,kappaArray,princNREL,numRelOrbs,orbMap)
    write(*,1)

    call readConfig
    write(*,1)

    call interpolateToDarc
    write(*,1)
    if(callOrthog) then
        call orthogonaliseGramSchmidt
    else 
        write(*,*) 'Orthogonalization turned off.'
    end if 
    write(*,1)
    call writeGrasp2018
    write(*,1)

    call writeTargetInp
    write(*,1)

    !call testingSmallR

    call writeORBINGrasp0
    call writeOutAutoOrbitals
    1 format('------------------------------------------------------------------------')

    contains

    subroutine optionalInput 
        namelist/autoDarc/ callOrthog

        callOrthog = .true.
        inquire(file='autoDarcInput',exist = inputExists)

        if(inputExists) then 
            open(1,file='autoDarcInput')
            read(1,autoDarc)
            close(1)
        end if 


    end subroutine

    subroutine writeOutAutoOrbitals 
        !writes out the Autostructure orbitals in a plott-able way.
        integer :: i ,j 

        radial(1) = 0.0d0 
        open(1,file='asorbitals.dat')
        write(1,3030) '  #Radial Grid',(i*2,princN(i),angular_symbols_neg(orbL(i)+1),'Orb  ',i*2+1,princN(i),angular_symbols_neg(orbL(i)+1),'Pot  ',i=1,num_orbs)
        do i = 1,num_points
            WRITE (1,3050) radial(i), (orbitals(j,i),potential(j,i),j=1,num_orbs)
        end do 

        close(1)
        3030 FORMAT(A16,1x,1000(I3,I3,A2,A5,3X))
        3050 FORMAT (1X,1P,1000E16.8)
    end subroutine


    subroutine orthogonaliseGramSchmidt
        !Implements the (modified) GramSchmidt procress as described on 
        !https://en.wikipedia.org/wiki/Gram–Schmidt_process

        !The code only orthogonalises radial orbitals of the same kappa.
        !This is because orbtials of different kappa are automatically 
        !orthogonal by the angular component. 
        !I briefly forgot this basic fact of Quantum Mechanics and lost some sleep.

        !I assume the first vector has already been normalised 
    


        integer :: ii,jj,numThisKappa,iii,jjj,kk
        real*8 :: overlap

        integer,allocatable :: orbMapKappa(:)
        integer :: targetKappa,numUniqueKappa
        integer,allocatable :: uniqueKappa(:)
        logical :: found 


        allocate(uniqueKappa(numRelOrbs))
        allocate(orbMapKappa(numRelOrbs))   
        
        write(*,1)

        orbMapKappa = 0 
        targetKappa = -1 
        numThisKappa=0
        
        numUniqueKappa = 1 
        uniqueKappa = kappaArray(1)
        !part of the code below could be moved in here priobably
        do ii = 2, numRelOrbs 
            targetKappa = kappaArray(ii)
            found = .false.
            do jj = 1,numUniqueKappa 
                if (uniqueKappa(jj).eq.targetKappa) then 
                    found = .true. 
                    exit 
                end if 
            end do 

            if (.not.found) then 
                numUniqueKappa = numUniqueKappa + 1
                uniqueKappa(numUniqueKappa) = targetKappa  
            end if 

        end do 
        !print*,'unique kappa',uniqueKappa(1:numUniqueKappa)
        !loop over all kappa groups and orthogonalise.
        write(*,2)
        do kk = 1, numUniqueKappa

            targetKappa = uniqueKappa(kk)
            numThisKappa = 0
            do ii = 1,numRelOrbs
                if(kappaArray(ii) .eq. targetKappa) then 
                    numThisKappa=numThisKappa+1
                    orbMapKappa(numThisKappa) = ii 
                end if 
            end do 

            do ii = 2,numThisKappa 

                iii = orbMapKappa(ii)
                !print*,princNREL(iii),kappaArray(iii)
                do jj=1,numThisKappa-1 

                    jjj = orbMapKappa(jj) 
                    overlap = overlapRadial(darcRadial,orbitalsSmall(iii,:),orbitalsLarge(iii,:),orbitalsSmall(jjj,:),orbitalsLarge(jjj,:))
                    orbitalsSmall(iii,:) = orbitalsSmall(iii,:) - overlap * orbitalsSmall(jjj,:)
                    orbitalsLarge(iii,:) = orbitalsLarge(iii,:) - overlap * orbitalsLarge(jjj,:)

                end do 
                overlap = normFactor(darcRadial,orbitalsSmall(iii,:),orbitalsLarge(iii,:))
                orbitalsSmall(iii,:) = orbitalsSmall(iii,:)/sqrt(overlap)
                orbitalsLarge(iii,:) = orbitalsLarge(iii,:)/sqrt(overlap)
            end do 
            write(*,4)
            do ii=1,numThisKappa 
                do jj=ii,numThisKappa 
                    iii = orbMapKappa(ii)
                    jjj = orbMapKappa(jj) 
                    overlap = overlapRadial(darcRadial,orbitalsSmall(iii,:),orbitalsLarge(iii,:),orbitalsSmall(jjj,:),orbitalsLarge(jjj,:))
                    write(*,3) princNREL(iii),angular_string(iii),princNREL(jjj),angular_string(jjj),overlap
                end do 
            end do
                



        end do 
        1 format (' Orthonormalizing same-kappa orbitals.')
        2 format (' Orthonormalization data:')
        3 format ('   ',I2,A2,1X,I2,A2,1X,ES10.3)
        4 format ('     Shells   Overlap')

    end subroutine

    subroutine readConfig 
        !Reads the 'CONFIG.DAT' output file of AUTOSTRUCTURE.
        !It needs this to get an approximation for generalised occupatuon numbers.

        !which is reuqired for DARC. 
        !It does this by simply taking an average of the non-relativisitic occupations,
        !and for l.ne.0 it splits the average into the l and l- orbitals weighted by the 
        !corresponding kappa. 


        integer :: numClosed, numOrbsTotal 
        integer :: numValence 
        integer :: numConfigs 
        integer :: ii 
        integer :: jj 

        integer :: upperKap
        integer :: lowerKap 
        integer :: tot 

        integer,allocatable :: princ(:)
        integer,allocatable :: orbitalL(:)
        
        integer*8,allocatable :: currentConfig(:)
        real*8,allocatable :: averageConfig(:)
        real*8,allocatable :: averageConfigRel(:)
        integer :: numRelLocal 
        open(1,file='CONFIG.DAT',FORM='FORMATTED')
        read(1,1) numClosed, numOrbsTotal
        !print*,numClosed,numOrbsTotal

        numValence = numOrbsTotal - numClosed 
        allocate(princ(numValence))
        allocate(orbitalL(numValence))
        
        read(1,2) (princ(ii),orbitalL(ii),ii=1,numValence)
        read(1,1) numConfigs

        !print*,numConfigs

        allocate(averageConfig(numValence))
        allocate(currentConfig(numValence))
        allocate(averageConfigRel(2*numValence)) !not an optimal dimension but its fine.
        !skip the min and max because i dont need it. 
        read(1,*)
        read(1,*)
        averageConfig = 0 
        do ii = 1,numConfigs
            read(1,3) (currentConfig(jj),jj=1,numValence)
            averageConfig = averageConfig + currentConfig
        end do 
        averageConfig = averageConfig / numConfigs 
        
        close(1)
        write(*,4)
        write(*,6)
        do ii = 1,numValence
            write(*,5) princ(ii),orbitalL(ii),averageConfig(ii)
        end do 
        write(*,7)

        numRelLocal = 0 
        do ii = 1,numValence
            if (orbitalL(ii).eq.0) then 
                write(*,5) princ(ii),-1,averageConfig(ii)
                numRelLocal = numRelLocal + 1
                averageConfigRel(numRelLocal) = averageConfig(ii)

            else
                upperKap =  orbitalL(ii)
                lowerKap = -orbitalL(ii)-1
                tot = (abs(lowerKap) + abs(upperKap)) !technically a factor of 2 but cancels a factor of 2 down here, small optimization. 
                numRelLocal = numRelLocal + 1
                averageConfigRel(numRelLocal) = averageConfig(ii)*abs(upperKap)/tot

                numRelLocal = numRelLocal + 1
                averageConfigRel(numRelLocal) = averageConfig(ii)*abs(lowerKap)/tot


                write(*,5) princ(ii),upperKap,averageConfigRel(numRelLocal-1)
                write(*,5) princ(ii),lowerKap,averageConfigRel(numRelLocal  )
                
            end if 
        end do 
        
        !set the rest of the generalised occupation numbers.
        do ii = 1,numRelLocal 
            genOcc(numRelOrbs+1-ii) = averageConfigRel(numRelLocal-ii+1)
        end do 


        1 format(2I5)
        2 format(10000(I3,I2))
        3 format(10000(I2))
        4 format(' Using average config as generalised occupation numbers.')
        6 format(' Average (non-relativistic) state: ')
        7 format(' Average (    relativistic) state: ')

        5 format(I3,1X,I2,ES10.3)
    end subroutine


    subroutine read_radout 

        !Reads the 'radout' (formatted) output file from autostructure.
        !The (kappa-averaged...) orbitals are stored in orbitals.
        !which has dimension num_orbs X num_points.
        !other quantities are stored.
        !dummy quantities such as point number are forgotten about.
        character*6,parameter :: atom = '      '
        character*6,parameter :: term = '      '
    
        character*6,parameter :: file ='radout'



        real*8 :: DA(5)
        character*4 :: tlbl       

        integer :: i,j,dummy 

        write(*,10)
        open(1,file=file,form='formatted')

        read(1,1) key,num_orbs,mb,num_points,nelec,nzed,(da(i),i=1,5),tlbl

        allocate(radial(num_points))
        allocate(charge(num_points))

        !read points. 
        do i = 1,num_points,2 
            read(1,2) key,dummy,radial(i),charge(i),dummy,radial(i+1),charge(i+1)
            !write(*,2) key,dummy,radial(i),charge(i),dummy,radial(i+1),charge(i+1)
        end do 

        allocate(orbitals(num_orbs,num_points))
        allocate(potential(num_orbs,num_points))

        allocate(orbitalParams(num_orbs))
        allocate(princN(num_orbs))
        allocate(orbL(num_orbs))
        allocate(orbitalEnergy(num_orbs))
        allocate(orbitalFloat(num_orbs))

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

            !lpm 26.03.25 
            call estimateAZ(int(i,8),orbitalParams(i))
            print*,orbitalParams(i)

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
    10    format(' Reading radout from AUTOSTRUCTURE')
    11    format(' Found ',I4,' non-relativistic orbitals with ',I5,' radial points.')
    end subroutine

    subroutine findROrbitals(orbL, princN, kappa,relprincN,numRelOrbs,orbMap)
        !calculates which relativistic orbitals are needed by the code. 

        integer :: orbL(:),princN(:)
        integer,allocatable :: kappa(:),relprincN(:),orbMap(:)
        integer :: numRelOrbs , ii 
        numRelOrbs = 0 

        !first scan to get optimal dimensions 
        do ii = 1,size(orbL)
            numRelOrbs = numRelOrbs+1 
            if (orbL(ii).ne.0) then 
                numRelOrbs = numRelOrbs + 1
            end if  
        end do 

        allocate(kappa(numRelOrbs))
        allocate(relprincN(numRelOrbs))
        allocate(orbMap(numRelOrbs))

        numRelOrbs = 0 
        do ii = 1,size(orbL)

            if (orbL(ii).ne.0) then 
                numRelOrbs = numRelOrbs + 1
                kappa(numRelOrbs) = orbL(ii) 
                relprincN(numRelOrbs) = princN(ii)
                orbMap(numRelOrbs) = ii 
            end if  

            numRelOrbs = numRelOrbs+1 
            kappa(numRelOrbs) = -orbL(ii) - 1 
            relprincN(numRelOrbs) = princN(ii)
            orbMap(numRelOrbs) = ii 
        end do 
        write(*,1) numRelOrbs
        !write(*,5)

        allocate(genOcc(numRelOrbs))
        !set all gen occ to max - and the readConfigs routine will take care of the rest.
        do ii = 1,numRelOrbs
            genOcc(ii) = 2.0 * abs(kappa(ii))
        end do 



    1   format(' There are ',I4, ' corresponding relativistic orbitals.')

    end subroutine

    subroutine estimateAZ(orbital_index,az)

        !calculate the best az to use
        !i.e - the best value of the derivative of the wavefunction at zero.
        !only makes sense to use this for s orbitals, really 
        !as strictly (in a non-relativistic setting) the derivative is zero 
        !for l>= 1 at r=0. 

        real*8,intent(out) :: az
        real*8 :: int1,int2 
        real*8 :: rmax 
        real*8,allocatable ::table0(:), table1(:),table2(:)
        integer :: ii,cutoff 
        integer*8 :: orbital_index

        integer :: l 
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
        !table1 = table0 * orbitals(orbital_index,1:cutoff)
        do ii = 1,cutoff
            table1(ii) = table0(ii) * orbitals(orbital_index,ii)
        end do
        table2 = table0 * table0


        do ii = 1,cutoff-1
            int1 = int1 + (table1(ii+1)+table1(ii)) * (radial(ii+1)-radial(ii))
            int2 = int2 + (table2(ii+1)+table2(ii)) * (radial(ii+1)-radial(ii))
        end do 
        print*,cutoff,radial(cutoff)
        az = int1/int2 

    end subroutine

    subroutine interpolateToDarc 
        !interpolates the as non-relativistic orbitals
        !onto a log grid, which is then used to calcualte the 
        !small component

        !numerical derivativees for Q are approximated by 
        !interpolating a small distance away from each 
        !required point and using a central difference approximation.

        !This is also what is done in COMPAS' rwfnmchfmcdf program.
        !See https://github.com/compas/grasp.

        real*8 :: h 

        integer, parameter :: darcDim = 200

        real*8 :: nf 

        integer :: ii,jj,i
        real*8,allocatable :: orbForDeriv(:,:),derivative
        real*8,allocatable :: rhoArray(:)
        REAL*8 ,allocatable :: potTest(:)
        REAL*8 ,allocatable :: den(:)

            !cff grid
        real*8 :: rho
        real*8 :: hff = 1.d0/16.d0
        integer,parameter :: n0 = 220
        real*8 :: hfgrid(n0+1)
        real*8 :: selectedOrb(n0+1)
        real*8 :: selectedgrid(n0+1)
        real*8 :: yy2(n0+1)
        integer :: closestInd,numHFpoints

        allocate(angular_string(numRelOrbs))
        

        !rnt = EXP (-65.0D00/16.0D00) / nzed
        !rnt = 1e-2
        !this is the spacing in logspace. set it so that 
        !the grid covers up to (around) the last point output by 
        !autostructure. no need to go  any further
        !i lost some of my sanity to a numerically broken down orbitsl
        !because of this.
        !let it be a lesson kids, check stuff. i guess.
        ! - lpm 22.02.25 

        !to interpoalte - need a small number of points
        !look at the CFF grid - performs the best
        rho =  -4.d0 
        do i = 1,n0 
            hfgrid(i+1) = exp(rho)/nzed
            rho = rho + hff
        end do 



        rnt=1e-3
        hff = log(radial(size(radial))/rnt +1.0d0) / darcDim

        allocate(darcRadial(darcDim))

        !allocate(yy2(num_points))

        allocate(rhoArray(size(darcRadial)))
        do ii = 1,size(darcRadial)
            rhoArray(ii) = (ii-1) * hff 
        end do 


        do ii = 1,darcDim
            darcRadial(ii) = rnt*(exp(hff*(ii-1))-1.0d0)
        end do

        h = darcRadial(2)/10000.0d0
        !allocate(asp(num_points))
        !allocate(bsp(num_points))
        !allocate(csp(num_points))
        allocate(orbitalsLarge(numRelOrbs,darcDim))
        allocate(orbitalsSmall(numRelOrbs,darcDim))

        allocate(orbForDeriv(darcDim,2))
        !Initialise
        orbitalsLarge = 0.0d0 
        orbitalsSmall = 0.0d0 
        write(*,1)
        allocate(potTest(darcDim))
        allocate(den(darcDim))
        potTest = -nzed / darcRadial
        
       
        do jj = 1,numRelOrbs

            !Select the radial points closeest to a more 
            !sparse grid. - allows for better numerical stability §
            !near the nucleus. 
            do ii = 1,n0 
                if (hfgrid(ii+1).lt.radial(size(radial))) then
                    closestInd = minloc(abs(radial - hfgrid(ii+1)),1)
                    selectedOrb(ii+1) = orbitals(orbMap(jj),closestInd)
                    selectedgrid(ii+1) = radial(closestInd)
                else 
                    numHFpoints = ii+1
                    exit
                end if 
            end do 
            do ii = numHFpoints+1,n0 
                selectedgrid(ii+1) = hfgrid(ii+1)
            end do



            orbForDeriv = 0.0d0 

            !orbitalEnergy In hartrees already 
            !needs to be minus though - bound orbital. 
            den = -1.0d0*orbitalEnergy(orbmap(jj)) - potTest

            den = 1.0 + 0.25d0*fineStruc*fineStruc*den
            
            !For s orbitals - use the derivative at zero to get a better 
            !interpolation.
            if (kappaArray(jj).eq.-1) then 
                !print*,'using more accuate',orbitalParams(orbmap(jj))
                call spline(selectedgrid,selectedorb,numHFpoints,orbitalParams(orbmap(jj)),0.0d0,yy2)
            else 
                call spline(selectedgrid,selectedorb,numHFpoints,0.0d0,0.0d0,yy2)
            end if 

            !Large component
            do ii = 2,darcDim
                call splint(selectedgrid,selectedorb,yy2,numHFpoints,darcRadial(ii),orbitalsLarge(jj,ii))
                call splint(selectedgrid,selectedorb,yy2,numHFpoints,darcRadial(ii)-h,orbForDeriv(ii,1))
                call splint(selectedgrid,selectedorb,yy2,numHFpoints,darcRadial(ii)+h,orbForDeriv(ii,2))
            end do 

            !Small component
            !SKIP POINT 1 = 0 
            do ii = 2,darcDim
                derivative = (orbForDeriv(ii,2) - orbForDeriv(ii,1))/(2.0d0*h)
                orbitalsSmall(jj,ii) = derivative + kappaArray(jj)*orbitalsLarge(jj,ii)/darcRadial(ii)
            end do 

            orbitalsSmall(jj,:) = orbitalsSmall(jj,:) * 0.5d0 *fineStruc

            !Normalise 
            NF = SQRT(normFactor(darcRadial,orbitalsSmall(jj,:),orbitalsLarge(jj,:)))
            orbitalsSmall(jj,:) = orbitalsSmall(jj,:)/NF
            orbitalsLarge(jj,:) = orbitalsLarge(jj,:)/NF


        end do 
        write(*,2) rnt,darcDim+1,hff
        write(*,1000)
        write(*,1020)

        !Write out info
        do ii = 1,numRelOrbs
            if     (kappaArray(ii) .ge. 0) then 
                angular_string(ii) = angular_symbols_pos(kappaArray(ii))
            elseif (kappaArray(ii) .le. 0) then
                angular_string(ii) = angular_symbols_neg(abs(kappaArray(ii)))
            else   
                angular_string(ii) = 'NU'
        end if 

            write(*,1010) princNREL(ii),angular_string(ii),kappaArray(ii),orbitalEnergy(orbMap(ii)) &
                ,darcDim+1,genOcc(ii)
        end do 

        !print*,size(orbitalsSmall)
        1 format(' Interpolating to a darc-style grid and calculating small components.')
        2 format(' The resulting log-grid has:',/'  RNT =',ES10.3,/,'  NP  =',i4,/,'  h   =',ES10.3)
        1000 FORMAT(' Final relativistic data: ')
        1020 FORMAT ('   Subshell',1X,'Kappa',5X,'Eigenv(Ha)',2X,'NPoints',4X,'GenOcc')
        1010 FORMAT (3X,3X,I3,A2,3X,I3,3X,ES12.3,4X,I5,F10.4)

    end subroutine 

   ! subroutine calcSmallCompoment()

    
    subroutine writeTargetInp
        !WRites target.inp for the darc codes.
        integer :: num_points,iter_grid,iter_orbs
        num_points = size(darcRadial)
        open(2,file='TARGET.INP')
        print*,'Writing to TARGET.INP'

        WRITE (2,3040) numRelOrbs,num_points
        WRITE (2,3050) (genOcc(iter_orbs),iter_orbs=1,numRelOrbs)
        WRITE (2,3050) (darcRadial(iter_grid),iter_grid=1,num_points)

        do iter_orbs = 1,numRelOrbs
            WRITE (2,3060) princNREL(iter_orbs),kappaArray(iter_orbs)
            WRITE (2,3050) (orbitalsLarge(iter_orbs,iter_grid),iter_grid=1,num_points)
            WRITE (2,3060) princNREL(iter_orbs),kappaArray(iter_orbs)
            WRITE (2,3050) (orbitalsSmall(iter_orbs,iter_grid),iter_grid=1,num_points)
        end do 

        close(2)

        3040 FORMAT (1X,2I7)
        3050 FORMAT (1X,1P,4E16.8)
        3060 FORMAT (1X,I4,2X,I4)

    end subroutine

    subroutine writeGrasp2018
        !writes a rwfn.out file in the format of grasp2018 - so one can 
        !read this into grasp2018 to either use as a starting point
        !or as their actual orbitals.
        !mostly for testing, but the unformatted file is relatively small

        !in the future i will implement an option for this to be switched off. 

        INTEGER :: ii,jj ,num_points
        num_points = size(darcRadial)
        print*,'Writing to rwfn.out'
        open (unit=9, file='rwfn.out',status='unknown',form='unformatted')
        write(9) 'G92RWF'
        do ii = 1,numRelOrbs 
            write(9) princNREL(ii),kappaArray(ii),orbitalEnergy(orbMap(ii)),num_points
            write(9) orbitalParams(orbMap(ii)),(orbitalsLarge(ii,jj),jj=1,num_points),(orbitalsSmall(ii,jj),jj=1,num_points)
            write(9) (darcRadial(jj),jj=1,num_points)
        end do
        close(9)

    end subroutine

    subroutine writeORBINGrasp0
        integer :: i ,num ,j
        real*8 :: pzi ,qzi 

        num = size(darcRadial)

        open(2,file='ORBIN.DAT',form='UNFORMATTED')
 
        do i = 1,numRelOrbs
            pzi = orbitalParams(orbmap(i))
            qzi = qzfunc(pzi,kappaArray(i))

            WRITE(2) angular_string(i),princNREL(i),kappaArray(i),num,nzed,hff,rnt, &
                     orbitalEnergy(orbmap(i)),pzi,qzi
            WRITE(2) (orbitalsLarge(I,j),j=1,num),(orbitalsSmall(I,j),j=1,num)
        end do 
        close(2)

    end subroutine


    function normFactor(radial,small,large)
        !norm factor of a P,Q relativistic orbital .
        real*8 :: radial(:),small(:),large(:)
        integer :: np, ii  
        real*8 :: normFactor,contribution 
        normFactor = 0.0d0 
        np = size(radial)
        do ii = 1, np-1
            contribution = large(ii)*large(ii) + small(ii)*small(ii) 
            contribution = contribution + large(ii+1) * large(ii+1) + small(ii+1) * small(ii+1)
            normFactor = normFactor + contribution * (radial(ii+1)-radial(ii))
        end do  
        normFactor = normFactor * 0.5d0 
        !rint*,normFactor,np
    end function

    function overlapRadial(radial,small1,large1,small2,large2)
        !overlap of a pair of dirac radial orbitals (P1,Q1) and (P2,Q2) .
        real*8 :: radial(:),small1(:),large1(:),small2(:),large2(:)
        integer :: np, ii  
        real*8 :: overlapRadial,contribution 
        overlapRadial = 0.0d0 
        np = size(radial)
        do ii = 1, np-1
            contribution = large1(ii)*large2(ii) + small1(ii)*small2(ii) 
            contribution = contribution + large1(ii+1) * large2(ii+1) + small1(ii+1) * small2(ii+1)
            overlapRadial = overlapRadial + contribution * (radial(ii+1)-radial(ii))
        end do  
        overlapRadial = overlapRadial * 0.5d0 
    end function

    function qzfunc(pz,kappa)
        real*8 :: pz ,qzfunc,factor
        integer :: kappa
        !real*8 :: Z 
        real*8,parameter :: c = 137.03599d0
        real *8 :: gama 

        gama = sqrt(kappa*kappa - (nzed/C)**2 )

        if (kappa .ge. 0.0d0) then 
            factor =  (C/nzed) * 1.0d0 / (kappa + gama)
        else  
            factor =  (nzed/C) * 1.0d0 / (kappa - gama)
        end if 
        !print*,kappa,gama,Z/C,factor
        qzfunc = pz * factor

    end function


        SUBROUTINE spline(x,y,n,yp1,ypn,y2)
            INTEGER n,NMAX
            DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
            PARAMETER (NMAX=500)
            INTEGER i,k
            DOUBLE PRECISION p,qn,sig,un,u(NMAX)
            !print*,'hello',n,yp1,ypn
            if (yp1.gt..99d30) then
              y2(1)=0.d0
              u(1)=0.d0
            else
              y2(1)=-0.5d0
              u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
            endif
            do 11 i=2,n-1
              sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
              p=sig*y2(i-1)+2.d0
              y2(i)=(sig-1.d0)/p
              u(i)=(6.d0*((y(i+1)-y(i))/(x(i+ 1)                        &
             -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*   &
              u(i-1))/p
        11    continue
            if (ypn.gt..99d30) then
              qn=0.d0
              un=0.d0
            else
              qn=0.5d0
              un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
            endif
            y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
            do 12 k=n-1,1,-1
              y2(k)=y2(k)*y2(k+1)+u(k)
        12    continue
            return
            END SUBROUTINE
        
            SUBROUTINE splint(xa,ya,y2a,n,x,y)
            INTEGER n
            DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
            INTEGER k,khi,klo
            DOUBLE PRECISION a,b,h
            klo=1
            khi=n
        1     if (khi-klo.gt.1) then
              k=(khi+klo)/2
              if(xa(k).gt.x)then
                khi=k
              else
                klo=k
              endif
            goto 1
            endif
            h=xa(khi)-xa(klo)
            if (h.eq.0.d0) print *, ' bad xa input in splint'
            a=(xa(khi)-x)/h
            b=(x-xa(klo))/h
            y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* &
              (h**2)/6.d0
            return
            END  SUBROUTINE
end

