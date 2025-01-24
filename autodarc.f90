!Program to parse autostructure data to the 
!darc R-matrix codes. 

!The code right now assumes an IC input, producing a list of Dirac-orbitals.
!The large components are taken as the AS orbital: 
!   P_nk = P_nl .
!The small components are taken to be: 
!   Q_nk = 0.5 * \alpha * (d/dr + \kappa/r) * P_nl 

!Using these non-relativistic orbitals from optimizing a 
!Breit-Pauli-Schrodinger Hamiltonian in a Dirac Hamiltonian
!is equivalent to the so-called RCI-P method, as breifly introduced in 
!Charlotte Froese Fischer et al 2016 J. Phys. B: At. Mol. Opt. Phys. 49 182004
!DOI 10.1088/0953-4075/49/18/182004

!This provides a way to use the much more efficient jj-coupled darc r-matrix codes.
!Instead of first doing an LS run and doing an expensive re-couple using stgjk. 


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

    real*8,allocatable :: genOcc(:)

    integer :: key,num_orbs,mb,nelec,jj,kk
    integer*4 :: num_points
    real*8 :: nzed ,a,lambda
    
    integer,allocatable :: kappaArray(:),princNREL(:),orbMap(:)
    integer :: numRelOrbs , ii 


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

    call orthogonaliseGramSchmidt
    write(*,1)
    call writeGrasp2018
    write(*,1)

    call writeTargetInp
    write(*,1)

    1 format('------------------------------------------------------------------------')

    contains

    subroutine orthogonaliseGramSchmidt
        !Implements the (modified) GramSchmidt procress as described on 
        !https://en.wikipedia.org/wiki/Gram–Schmidt_process

        !The code only orthogonalises radial orbitals of the same kappa.
        !This is because orbtials of different kappa are automatically 
        !orthogonal by the angular component. 
        !I briefly forgot this basic fact of Quantum Mechanics and lost some sleep.

        !I assume the first vector has already been normalised 
    


        integer :: ii,jj,numThisKappa,iii,jjj,kk
        real*8 :: norm ,overlap

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
        character*3 :: el
    
        character*6,parameter :: file ='radout'

        real*8,allocatable :: smallComp(:)


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

        4 format('Write out output:')


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

    subroutine interpolateToDarc 
        !interpolates the as non-relativistic orbitals
        !onto a log grid, which is then used to calcualte the 
        !small component

        !numerical derivativees for Q are approximated by 
        !interpolating a small distance away from each 
        !required point and using a central difference approximation.

        !This is also what is done in COMPAS' rwfnmchfmcdf program.
        !See https://github.com/compas/grasp.


        real*8,allocatable :: asp(:),bsp(:),csp(:)
        real*8 :: sev
        real*8 :: h 

        integer, parameter :: darcDim = 250
        real*8 :: rnt ,NF 
        real*8 :: hff 
        integer :: ii,jj,kk
        real*8,allocatable :: newOrb(:),smallComp(:)
        real*8,allocatable :: orbForDeriv(:,:),derivative


        character*2 :: angular_symbols_neg(5) = (/'s ','p ','d ','f ','g ' /)
        character*2 :: angular_symbols_pos(4) =      (/'p-','d-','f-','g-' /)
        allocate(angular_string(numRelOrbs))
        

        rnt = EXP (-65.0D00/16.0D00) / nzed
        
        !this is the spacing in logspace. set it so that 
        !the grid covers up to (around) the last point output by 
        !autostructure. no need to go  any further
        !i lost some of my sanity to a numerically broken down orbitsl
        !because of this.
        !let it be a lesson kids, check stuff. i guess.
        ! - lpm 22.04.25 
        hff = log(radial(size(radial))/rnt +1.0d0) / (darcDim)


        allocate(darcRadial(darcDim+1))

        do ii = 1,darcDim+1
            darcRadial(ii) = rnt*(exp(hff*(ii-1))-1.0d0)
        end do

        h = darcRadial(2)/1000.0d0
        allocate(asp(num_points))
        allocate(bsp(num_points))
        allocate(csp(num_points))
        allocate(orbitalsLarge(numRelOrbs,darcDim+1))
        allocate(orbitalsSmall(numRelOrbs,darcDim+1))

        allocate(orbForDeriv(darcDim+1,2))
        !Initialise
        orbitalsLarge = 0.0d0 
        orbitalsSmall = 0.0d0 
        write(*,1)
        do jj = 1,numRelOrbs
            orbForDeriv = 0.0d0 
            call SPLINE(num_points,radial,orbitals(orbMap(jj),:),asp,bsp,csp)
            !Large component
            do ii = 2,darcDim+1
                orbitalsLarge(jj,ii) = SEVAL(num_points,darcRadial(ii),radial,orbitals(orbMap(jj),:),asp,bsp,csp)
                !print*,orbMap(jj)

                orbForDeriv(ii,1) = SEVAL(num_points,darcRadial(ii)-h,radial,orbitals(orbMap(jj),:),asp,bsp,csp)
                orbForDeriv(ii,2) = SEVAL(num_points,darcRadial(ii)+h,radial,orbitals(orbMap(jj),:),asp,bsp,csp)
            end do 
            !Small component
            do ii = 2,darcDim+1
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
        !write(*,5)
        write(*,1000)
        write(*,1020)
        do ii = 1,numRelOrbs
            if     (kappaArray(ii) .ge. 0) then 
                angular_string(ii) = angular_symbols_pos(kappaArray(ii))
            elseif (kappaArray(ii) .le. 0) then
                angular_string(ii) = angular_symbols_neg(abs(kappaArray(ii)))
            else   
                !in principal this should never be reached. only if kappa==0
                angular_string(ii) = 'NULL'
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


    !the following two routines are borrowed from the below reference.

    !*                                                      *
    !* ---------------------------------------------------- *
    !* Ref.: From Numath Library By Tuan Dang Trong in      *
    !*       Fortran 77 [BIBLI 18].                         *
    !*                                                      *
    !*                   F90 Release By J-P Moreau, Paris.  *
    !*                           (www.jpmoreau.fr)          *
    !********************************************************
    FUNCTION SEVAL (N,U,X,Y,B,C,D)
        !------------------------------------------------------------------------
        !     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
        !     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
        !     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
        !     BY THE SPLINE SUBROUTINE.
        !
        !     INPUTS:
        !     N       NUMBER OF POINTS OF CURVE Y = F(X)
        !     U       ABSCISSA OF POINT TO BE INTERPOLATED
        !     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
        !             OF CURVE F(X)
        !     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
        !             CUBIC SPLINE
        !
        !     OUTPUTS:
        !     SEVAL   INTERPOLATED VALUE
        !             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
        !             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
        !
        !     REFERENCE :
        !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
        !     COMPUTATIONS. PRENTICE-HALL,INC.
        !------------------------------------------------------------------------
              REAL *8 B(N),C(N),D(N),X(N),Y(N),U,DX,SEVAL
              integer :: n,I,j,K
              DATA I/1/

        !     BINARY SEARCH

              IF (I.GE.N) I = 1
              IF (U.LT.X(I)) GO TO 10
              IF (U.LE.X(I+1)) GO TO 30
           10 I = 1
              J = N+1
           20 K = (I+J)/2
              IF (U.LT.X(K)) J = K
              IF (U.GE.X(K)) I = K
              IF (J.GT.I+1) GO TO 20

        !     SPLINE EVALUATION

           30 DX = U-X(I)
              SEVAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
              RETURN
        END

        SUBROUTINE SPLINE (N,X,Y,B,C,D)
            !---------------------------------------------------------------------
            !     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
            !     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
            !
            !     INPUTS:
            !     N       NUMBER OF GIVEN POINTS
            !     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
            !             OF FUNCTION F(X)
            !
            !     OUTPUTS:
            !     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
            !             OF THE CUBIC SPLINE
            !
            !     REFERENCE:
            !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
            !     COMPUTATIONS. PRENTICE-HALL,INC.
            !---------------------------------------------------------------------
                  IMPLICIT REAL *8 (A-H,O-Z)
                  integer :: n
                  real*8:: B(N),C(N),D(N),X(N),Y(N)
                  integer :: i,l,j,k,nm1
                  NM1 = N-1
                  IF (N.LT.2) RETURN
                  IF (N.LT.3) GO TO 50
            
            !     BUILD THE TRIDIAGONAL SYSTEM
            !     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)
            
                  D(1) = X(2)-X(1)
                  C(2) = (Y(2)-Y(1))/D(1)
                  DO 10 I = 2,NM1
                  D(I) = X(I+1)-X(I)
                  B(I) = 2.D0*(D(I-1)+D(I))
                  C(I+1) = (Y(I+1)-Y(I))/D(I)
                  C(I) = C(I+1)-C(I)
               10 CONTINUE
            
            !     CONDITIONS AT LIMITS
            !     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES
            
                  B(1) = -D(1)
                  B(N) = -D(N-1)
                  C(1) = 0.D0
                  C(N) = 0.D0
                  IF (N.EQ.3) GO TO 15
                  C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
                  C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
                  C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
                  C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
            
            !     FORWARD ELIMINATION
            
               15 DO 20 I = 2,N
                  T = D(I-1)/B(I-1)
                  B(I) = B(I)-T*D(I-1)
                  C(I) = C(I)-T*C(I-1)
               20 CONTINUE
            
            !     BACK SUBSTITUTION
            
                  C(N) = C(N)/B(N)
                  DO 30 L = 1,NM1
                  I = N-L
                  C(I) = (C(I)-D(I)*C(I+1))/B(I)
               30 CONTINUE
            
            !     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL
            
                  B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
                  DO 40 I = 1,NM1
                  B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
                  D(I) = (C(I+1)-C(I))/D(I)
                  C(I) = 3.D0*C(I)
               40 CONTINUE
                  C(N) = 3.D0*C(N)
                  D(N) = D(NM1)
                  RETURN
            
            !     CAS N = 2
            
               50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
                  C(1) = 0.D0
                  D(1) = 0.D0
                  B(2) = B(1)
                  C(2) = 0.D0
                  D(2) = 0.D0
                  RETURN
        END
            
            ! end of file tseval.f90

end program