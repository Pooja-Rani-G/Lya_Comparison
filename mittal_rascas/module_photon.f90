!!Shikhar: This version of module_photon accounts for the Hubble flow (provided by Thibault Garel) and calculates rate of scattering, P_alpha.
module module_photon
 
  use module_gas_composition
  use module_mesh
  use module_constants
  use module_random
  use module_domain
  use module_utils, only: path, get_voigt_over_dvoigt

  implicit none

  ! 2 types for photons, one for the initial properties is called photon_init
  ! and one for the properties that evolve during the RT called photon_current

  type photon_current
     integer(kind=4)           :: ID
     integer(kind=4)           :: status       ! =0 if flying, =1 if escape, =2 if absorption (by dust)
     real(kind=8),dimension(3) :: xlast        ! coordinates of last interaction in box units
     real(kind=8),dimension(3) :: xcurr        ! current position of the photon in box units
     real(kind=8)              :: nu_ext       ! external frame frequency (Hz)
     real(kind=8),dimension(3) :: k            ! normalised propagation vector 
     integer(kind=4)           :: nb_abs       ! number of interactions before escape
     real(kind=8)              :: time         ! time in [s] from emission to escape/absorption        
     real(kind=8)              :: tau_abs_curr ! current optical depth (useful when photon change mesh domain)
     integer(kind=4)           :: iran         ! state of the random generator
!===========================================================================================
     real(kind=8)              :: total_OD     ! total optical depth (emission to escape)
!===========================================================================================
  end type photon_current
  ! Note: if you change something in photon_current, don't forget to update the mpi_photon_type in module_parallel_mpi.f90


  type photon_init
     integer(kind=4)           :: ID
     real(kind=8)              :: nu_em    ! emitted frequency in external frame
     real(kind=8),dimension(3) :: x_em
     real(kind=8),dimension(3) :: k_em
     integer(kind=4)           :: iran     ! state of the random generator
  end type photon_init

  
  public  :: MCRT, propagate, init_photons_from_file, dump_photons
  
contains

!!Shikhar: added 'lev', 'palpha_local' and 'one_plus_z' here. 'palpha_local' will be called back to module_worker.f90 in case of parallel
!! version and rascas-serial.f90 in case of serial version.
  subroutine MCRT(npp,photpacket,mesh_dom,compute_dom,palpha_local,lev,one_plus_z)

    ! this is the Monte Carlo Radiative transfer routine... but also a single loop over photons...

    integer(kind=4),intent(in)                        :: npp
    type(photon_current),dimension(npp),intent(inout) :: photpacket
    type(mesh),intent(in)                             :: mesh_dom
    type(domain),intent(in)                           :: compute_dom
    integer(kind=4)                                   :: i
    !=============================================================================================
    integer(kind=4),intent(in)    :: lev
    real(kind=8),intent(in)       :: one_plus_z
    real(kind=8),intent(inout)    :: palpha_local(:,:,:)
    !=============================================================================================
    
    do i=1,npp
       ! case of photpacket not fully filled...
       if (photpacket(i)%ID>0) then
          call propagate(photpacket(i),mesh_dom,compute_dom,palpha_local,lev,one_plus_z)
       endif
    enddo
    
  end subroutine MCRT

!!Shikhar: will store palpha for an individual photon in 'pa_single_photon'
  subroutine propagate(p,domesh,domaine_calcul,pa_single_photon,lev,one_plus_z)    
    type(photon_current),intent(inout)   :: p              ! photon 
    type(mesh),intent(in)                :: domesh         ! mesh
    type(domain),intent(in)              :: domaine_calcul ! computational domain in which photons are propagating
    
    type(gas)                            :: cell_gas                   ! gas in the current cell 
    integer(kind=4)                      :: icell, ioct, ind, ileaf, cell_level  ! current cell indices and level
    real(kind=8)                         :: cell_size, cell_size_cm, tau_abs, scalar, nu_cell, nu_ext, rtau
    real(kind=8),dimension(3)            :: ppos,ppos_cell             ! working coordinates of photon (in box and in cell units)
    real(kind=8)                         :: distance_to_border,distance_to_border_cm, d,distance_to_border_box_units
    real(kind=8)                         :: time
    integer(kind=4)                      :: scatter_flag, i, icellnew, iran, npush
    real(kind=8),dimension(3)            :: vgas, k, cell_corner, posoct, pcell
    logical                              :: cell_fully_in_domain, flagoutvol, in_domain, OutOfDomainBeforeCell
    real(kind=8)                         :: dborder, dborder_cm, error

    !===========================================================================================
    real(kind=8),intent(inout)  :: pa_single_photon(:,:,:)
    real(kind=8),intent(in)     :: one_plus_z
    integer(kind=4),intent(in)  :: lev
    integer(kind=4)             :: ix,jy,kz
    real(kind=8)                :: Vol, num_den, TOD, tau_substep
    !=========================================================================================== 

    !--CORESKIP--
    real(kind=8) :: delta_nu_doppler,a,xcw,nu_0
    !--PIKSEROC--
    
    ! HUBBLE-FLOW                                                                                                      
    real(kind=8)                         :: Hub_cms_cm, dmax_hub1, dmax_hub2, voigt_over_dvoigt
    real(kind=8)                         :: x_cell, true_distance_to_border_cm, sum_substeps_distance_cm, substep_distance_cm,remaining_distance_to_border_cm,remaining_distance_to_border,should_be_zero
    integer(kind=4)                      :: nsubsteps
    real(kind=8),parameter               :: ntmax = 10.0d0
    real(kind=8)                         :: Hub_kms_Mpc ! [km/s/Mpc]
    
    
    Hub_kms_Mpc = H_0 * sqrt(Omega_M*one_plus_z**3.0 + Omega_L)
    Hub_cms_cm  = Hub_kms_Mpc * 1.0d5 / mpc                             ! [cm/s/cm]       
    ! WOLF-HUBBLE

    ! initialise working props of photon
    ppos    = p%xcurr        ! position within full simulation box, in box units.
    time    = p%time
    tau_abs = p%tau_abs_curr
    iran    = p%iran
!=============================================================================================    
    TOD     = p%total_OD
!=============================================================================================

    ! find cell in which the photon is, and define all its indices
    icell = in_cell_finder(domesh,ppos)
    if(domesh%son(icell)>=0)then
       print*,'ERROR: not a leaf cell'
       stop
    endif
    ileaf = - domesh%son(icell)
    ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
    ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct
    flagoutvol = .false.
   
    ! propagate photon until escape or death ... 
    photon_propagation : do 
       
       ! gather properties properties of current cell
       cell_level   = domesh%octlevel(ioct)      ! level of current cell
       cell_size    = 0.5d0**cell_level          ! size of current cell in box units
       cell_size_cm = cell_size * box_size_cm    ! size of the current cell in cm
       cell_gas     = domesh%gas(ileaf)
       ! compute position of photon in current-cell units
       posoct(:)    = domesh%xoct(ioct,:)
       cell_corner  = get_cell_corner(posoct,ind,cell_level)   ! position of cell corner, in box units.
       ppos_cell    = (ppos - cell_corner) / cell_size         ! position of photon in cell units (x,y,z in [0,1] within cell)

       if((ppos_cell(1)>1.0d0).or.(ppos_cell(2)>1.0d0).or.(ppos_cell(3)>1.0d0).or. &
            (ppos_cell(1)<0.0d0).or.(ppos_cell(2)<0.0d0).or.(ppos_cell(3)<0.0d0))then
          print*,"ERROR: problem in computing ppos_cell"
          stop
       endif

!==========================================================================       
       ix=ppos(1)/0.5**lev+1
       jy=ppos(2)/0.5**lev+1
       kz=ppos(3)/0.5**lev+1
       Vol = cell_size_cm**3
       num_den = cell_gas%nHI
!==========================================================================
       vgas                = get_gas_velocity(cell_gas)
       scalar              = p%k(1) * vgas(1) + p%k(2) * vgas(2) + p%k(3) * vgas(3)
       nu_cell             = (1.0d0 - scalar/clight) * p%nu_ext

       ! define/update flag_cell_fully_in_comp_dom to avoid various tests in the following
       pcell = cell_corner + 0.5d0*cell_size
       cell_fully_in_domain = domain_fully_contains_cell(pcell,cell_size,domaine_calcul)
     
       !--CORESKIP--
       !if (HI_core_skip) then
       ! Needed for HUBBLE-FLOW  
       delta_nu_doppler = cell_gas%dopwidth/(1215.67d0/cmtoA)
       a    = 6.265d8/fourpi/delta_nu_doppler
       xcw  = 6.9184721d0 + 81.766279d0 / (log10(a)-14.651253d0)  ! Core/wing limit in x units ! Smith+15, Eq. 21
       nu_0 = clight /(1215.67d0/cmtoA)
       !end if
       !--PIKSEROC--
       
       propag_in_cell : do

          ! generate the opt depth where the photon is scattered/absorbed
          if (tau_abs <= 0.0d0) then
             rtau    = ran3(iran)
             tau_abs = -log(1.0d0-rtau+1.0d-30)
             !=====================================
             TOD = TOD + tau_abs
             !=====================================
          end if

          ! compute distance of photon to border of cell along propagation direction
          distance_to_border           = path(ppos_cell,p%k)                   ! in cell units
          distance_to_border_cm        = distance_to_border * cell_size_cm     ! cm
          distance_to_border_box_units = distance_to_border * cell_size        ! in box units
          ! if cell not fully in domain, modify distance_to_border to "distance_to_domain_border" if relevant
!          OutOfDomainBeforeCell = .False.
          if(.not.(cell_fully_in_domain))then
             dborder    = domain_distance_to_border_along_k(ppos,p%k,domaine_calcul)  ! in box units
             dborder_cm = dborder * box_size_cm                                       ! from box units to cm
             ! compare distance to cell border and distance to domain border and take the min
             if (dborder_cm < distance_to_border_cm) then
                OutOfDomainBeforeCell        = .True.
                distance_to_border_cm        = dborder_cm
                distance_to_border_box_units = dborder
                distance_to_border           = distance_to_border_cm / cell_size_cm
             end if
          endif


          if ((distance_to_border < 0.0d0) .or. (distance_to_border > sqrt(3.001d0))) then
             print*,'Dborder... ',distance_to_border
             stop
          end if

          ! HUBBLE-FLOW
          ! Propagation in Hubble flow is made with an adaptive scheme by advancing photons over substeps (dr = dv/H(z)) within cells.
          ! Probability of interaction with HI is a sharp function in x (i.e. Voigt profile, V(a,x), centered at x=0). We thus need that : 1/ dv to be small enough wrt to b (otherwise interaction in the core can be missed), and 2/ V(a,x) has a steep form, esp. in the core, so dv must be small enough to evalutate V(a,x) at correct x.
          ! 1st condition : dv < b / 5, i.e. dr < b / (5H(z))
          ! b is called dopwidth here, and is in cm/s 
          dmax_hub1 = cell_gas%dopwidth / (5.0 * Hub_cms_cm)  ! cm
          ! 2nd condition: dr < 0.01 b/H(z) * V(a,x)/V’(a,x), where V(a,x) is the Voigt function
          ! See line 306 below
          
!!$          delta_nu_doppler = cell_gas%dopwidth/(1215.67d0/cmtoA)
!!$          a    = 6.265d8/fourpi/delta_nu_doppler
!!$          ! Compute the wing/core transition xcw for Lya
!!$          xcw  = 6.9184721d0 + 81.766279d0 / (log10(a)-14.651253d0)  ! Smith+15, Eq. 21
!!$          nu_0 = clight /(1215.67d0/cmtoA)
          
          ! Initialize distances
          true_distance_to_border_cm      = distance_to_border_cm
          remaining_distance_to_border_cm = distance_to_border_cm
          
          sum_substeps_distance_cm = 0.0
          scatter_flag = 0
          nsubsteps = 0

          
          do while (scatter_flag .eq. 0 .and. sum_substeps_distance_cm < true_distance_to_border_cm)
             x_cell              = (nu_cell - nu_0) / delta_nu_doppler
             ! Evaluate V/V'
             voigt_over_dvoigt   = get_voigt_over_dvoigt(x_cell,xcw)         
             dmax_hub2           = 0.01 * cell_gas%dopwidth / Hub_cms_cm * voigt_over_dvoigt ! cm
             substep_distance_cm = min(remaining_distance_to_border_cm,dmax_hub1,dmax_hub2)

!=====================================================================================================================
             scalar  = Hub_cms_cm * substep_distance_cm
             p%nu_ext            = (1.0d0 - scalar/clight) * p%nu_ext
             scalar              = p%k(1) * vgas(1) + p%k(2) * vgas(2) + p%k(3) * vgas(3)
             !! compute photon's frequency in cell's moving frame
             nu_cell             = (1.0d0 - scalar/clight) * p%nu_ext
             
             ! check whether scattering occurs within cell or domain (scatter_flag > 0) or not (scatter_flag==0)
             ! The statement above will no longer be true with substepping
             !scatter_flag = gas_get_scatter_flag(cell_gas, substep_distance_cm, nu_cell_temp, tau_abs, iran)
             scatter_flag = gas_get_scatter_flag(cell_gas, substep_distance_cm, nu_cell, tau_abs, iran,tau_substep)
             pa_single_photon(ix,jy,kz) = pa_single_photon(ix,jy,kz) + tau_substep/(Vol*num_den)
!=====================================================================================================================
             
             ! Update total "travelled" distance" and "remaining" distance
             sum_substeps_distance_cm = sum_substeps_distance_cm + substep_distance_cm
             remaining_distance_to_border_cm    = true_distance_to_border_cm-sum_substeps_distance_cm             
             distance_to_border_cm    = sum_substeps_distance_cm ! -sum_substeps_distance_cm             

             !! remaining_distance_to_border_cm can be < 0 due to numerical precision
             !! If so set it to 0, set distance_to_border_cm to true_distance_to_border_cm and exit while loop
             if (remaining_distance_to_border_cm < 0.0d0) then
                ! Check that it is indeed due to num. prec. 
                remaining_distance_to_border = remaining_distance_to_border_cm / cell_size_cm ! cell units
                if (abs(remaining_distance_to_border) > 1.d-14) then
                   print*,'remaining_distance_to_border < 0 but nut due to num. prec. ... : ',remaining_distance_to_border,remaining_distance_to_border_cm
                   stop
                else
                   remaining_distance_to_border_cm = 0.0d0
                   distance_to_border_cm = true_distance_to_border_cm
                   exit
                end if
             endif
                          
             nsubsteps = nsubsteps + 1
             
          end do ! end substepping
          
          ! Recompute distance in other units...
          distance_to_border = distance_to_border_cm / cell_size_cm     ! cell units
          distance_to_border_box_units = distance_to_border * cell_size ! box units

          ! WOLF-HUBBLE
          
          if (scatter_flag == 0) then   ! next scattering event will not occur in the cell or in the domain

             ! Check that true_distance_to_border_cm=distance_to_border_cm if no scattering in this cell
             !if (distance_to_border_cm .ne. true_distance_to_border_cm) then
             should_be_zero = abs(true_distance_to_border_cm/distance_to_border_cm) - 1.0d0
             if (should_be_zero > 1.0d-14) then
                print*,'Problem with distance_to_border_cm due to substeps...: ',true_distance_to_border_cm,distance_to_border_cm,should_be_zero
                stop
             end if
             
             ! move photon out of cell or domain
             ppos = ppos + p%k * distance_to_border_box_units *(1.0d0 + epsilon(1.0d0))
             
             ! correct for periodicity
             do i=1,3
                if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
                if (ppos(i) .ge. 1.0d0) ppos(i)=ppos(i)-1.0d0
             enddo
             ! update travel time
             time = time + distance_to_border_cm/clight

!             if (OutOfDomainBeforeCell) then ! photon exits computational domain and is done 
                ! it may happen due to numerical precision that the photon is still in the domain despite epsilon above.
                ! -> check and issue warning if it is the case. The error should not be larger than a few times epsilon. 
!                in_domain = domain_contains_point(ppos,domaine_calcul)
!                if (in_domain) then
!                   if (domain_distance_to_border_along_k(ppos,p%k,domaine_calcul)>3.d0*epsilon(distance_to_border)) then
!Shikhar: I suppressed the following warning.
!                      print*,'WARNING : photon still in domain when it should not ... '
!                      error = nint(domain_distance_to_border_along_k(ppos,p%k,domaine_calcul)/epsilon(distance_to_border))
!            !          print*,'          (error ~ ',error,' times num. prec.) '
!                   end if
!                end if
!                p%status       = 1
!                p%xcurr        = ppos
!                p%time         = time
!                p%tau_abs_curr = tau_abs
!                p%iran         = iran
                !==================================
!                p%total_OD     = TOD - tau_abs
                !==================================
!                exit photon_propagation

!             else
                ! photon exits current cell -> find into which new cell it goes
                call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
                ! It may happen due to numerical precision that the photon is still in the current cell (i.e. icell == icellnew).
                ! -> give it an extra push untill it is out. 
                npush = 0
                do while (icell==icellnew)
                   npush = npush + 1
                   ppos(1) = ppos(1) + merge(-1.0d0,1.0d0,p%k(1)<0.0d0) * epsilon(ppos(1))
                   ppos(2) = ppos(2) + merge(-1.0d0,1.0d0,p%k(2)<0.0d0) * epsilon(ppos(2))
                   ppos(3) = ppos(3) + merge(-1.0d0,1.0d0,p%k(3)<0.0d0) * epsilon(ppos(3))
                   ! correct for periodicity
                   do i=1,3
                      if (ppos(i) < 0.0d0) ppos(i)=ppos(i)+1.0d0
                      if (ppos(i) .ge. 1.0d0) ppos(i)=ppos(i)-1.0d0
                   enddo
                   
                   call whereIsPhotonGoing(domesh,icell,ppos,icellnew,flagoutvol)
                end do
!=============================================================================
                !if (npush > 1) print*,'WARNING : npush > 1 needed in module_photon:propagate.'
                ! test whether photon was pushed out of domain with the extra pushes
                ! (and in that case, call it done). 
                !if (npush > 0) then 
                !   in_domain = domain_contains_point(ppos,domaine_calcul)
                !   if (.not. in_domain) then
                !      print*,'WARNING: pushed photon outside domain ... '
                !      p%status       = 1
                !      p%xcurr        = ppos
                !      p%time         = time
                !      p%tau_abs_curr = tau_abs
                !      p%iran         = iran
                !      p%total_OD     = TOD - tau_abs
                !      exit photon_propagation
                !   end if
                !end if
!=============================================================================
                ! check if the new cell is outside of the current cpu domain (flagoutvol)
                ! And if so send it back to master
                if(flagoutvol)then
                   p%xcurr        = ppos
                   p%time         = time
                   p%tau_abs_curr = tau_abs
                   p%iran         = iran
                   !==================================
                   p%total_OD     = TOD - tau_abs
                   !==================================
                   exit photon_propagation
                endif
                ! Finally, if we're here, the photon entered a cell within the current cpu domain so we go on. 
                icell = icellnew
                ileaf = - domesh%son(icell)
                ind   = (icell - domesh%nCoarse - 1) / domesh%nOct + 1   ! JB: should we make a few simple functions to do all this ? 
                ioct  = icell - domesh%nCoarse - (ind - 1) * domesh%nOct

                ! there has been no interaction in the cell, tau_abs has been updated in gas_get_scatter_flag
                ! -> move to next cell
                exit propag_in_cell

!             end if


          else
             ! Next event happens inside this cell and in the domain.

             ! length and time travelled by the photon before event
             d    = distance_to_border_cm   ! NB: at this point, distance_to_border became "distance_to_interaction" in gas_get_scatter_flag
             time = time + d/clight
             d    = d / cell_size_cm        ! in cell units

             ! update ppos_cell
             do i=1,3
                ppos_cell(i) = ppos_cell(i) + p%k(i) * d
             enddo
             ! update ppos according to ppos_cell
             ppos = ppos_cell * cell_size + cell_corner

             !------------
             ! scattering
             !------------

             p%nb_abs = p%nb_abs + 1     ! increment nb of scatterings
             p%xlast = ppos              ! memorize the location of "potential" last interaction
             ! a scattering event modifies nu_cell, k, and nu_ext
             ! so it needs to transport nu_cell, k, nu_ext, but not type p (because type p not known in gas)
             nu_ext = p%nu_ext
             k = p%k
             call gas_scatter(scatter_flag, cell_gas, nu_cell, k, nu_ext, iran)
             p%nu_ext = nu_ext             
             ! NB: for TEST case, to have photons propagating straight on, comment the following line
             p%k = k
             ! there has been an interaction -> reset tau_abs
             tau_abs = -1.0d0

!=====================================================================
!Shikhar: set here the condition for the medium to become transparent.
             x_cell = (nu_cell - nu_0) / delta_nu_doppler
             if (x_cell .le. -8.0d0) then
                p%status       = 1
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
                p%total_OD     = TOD
                exit photon_propagation
             end if
!=====================================================================
             ! scatter_flag allows to know the status (aborbed or not) of the photon in case of dust
             ! new convention: negative if absorbed
             if(scatter_flag<0)then
                ! photon has been absorbed by dust, photon done, nothing else to do
                p%status       = 2
                p%xcurr        = ppos
                p%time         = time
                p%tau_abs_curr = tau_abs
                p%iran         = iran
                !==================================
                p%total_OD     = TOD - tau_abs
                !==================================
                exit photon_propagation
             endif

          end if

       end do propag_in_cell

    end do photon_propagation   

    ! End of the photon propagation. There are 3 possible cases:
    !   1/ photon is out of the computational domain == escaped           -> in_domain=.false. && p%status=1
    !   2/ photon is out of mesh-cpu domain -> sent back to master, etc.  -> flagoutvol==.true.
    !   3/ photon is dead                                                 -> p%status=2
    
    ! some simple sanity checks
    if(.not.(flagoutvol).and.(p%status==0))then
       print*,'ERROR: problem 1 with photon propagation in module_photon.f90!',flagoutvol,p%status
       stop
    endif

    if(.not.(flagoutvol.or.(p%status==1).or.(p%status==2)))then
       print *,'ERROR: problem 2 with photon propagation in module_photon.f90!',flagoutvol,p%status
       stop
    endif

  end subroutine propagate


!========================================================================
!Shikhar: added the extra parameter total_flux in the following function
!========================================================================
  subroutine init_photons_from_file(file,pgrid,total_flux)

    character(2000),intent(in)                                 :: file
    type(photon_current),dimension(:),allocatable, intent(out) :: pgrid
    type(photon_init),dimension(:),allocatable                 :: pgridinit
    integer(kind=4)                                            :: i, n_photon, iseed
    real(kind=8)                                               :: knorm
!=====================================================================================    
    real(kind=8), intent(out)                                  :: total_flux
!=====================================================================================
    ! read ICs
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='read')
    read(14) n_photon
    read(14) total_flux    ! nb of real photons [# / s]
    allocate(pgridinit(n_photon))
    read(14) iseed
    read(14) (pgridinit(i)%ID,i=1,n_photon)
    read(14) (pgridinit(i)%nu_em,i=1,n_photon)
    read(14) (pgridinit(i)%x_em(:),i=1,n_photon)
    read(14) (pgridinit(i)%k_em(:),i=1,n_photon)
!=================================================
!    read(14) (pgridinit(i)%iran,i=1,n_photon)
!=================================================
    close(14)
!=================================================
    pgridinit%iran=-pgridinit%ID
!=================================================

    ! build photgrid current
    allocate(pgrid(n_photon))
    do i=1,n_photon
       pgrid(i)%ID           = pgridinit(i)%ID
       pgrid(i)%status       = 0
       pgrid(i)%xlast        = pgridinit(i)%x_em
       pgrid(i)%xcurr        = pgridinit(i)%x_em
       pgrid(i)%nu_ext       = pgridinit(i)%nu_em
       ! make sure k is normalised
       knorm                 = sqrt(pgridinit(i)%k_em(1)**2 + pgridinit(i)%k_em(2)**2 + pgridinit(i)%k_em(3)**2)
       pgrid(i)%k            = pgridinit(i)%k_em / knorm
       pgrid(i)%nb_abs       = 0
       pgrid(i)%time         = 0.0d0
       pgrid(i)%tau_abs_curr = -1.0d0
       pgrid(i)%iran         = pgridinit(i)%iran
       !=========================================================
       pgrid(i)%total_OD     = 0.0d0
       !=========================================================
    enddo
    deallocate(pgridinit)

  end subroutine init_photons_from_file



  subroutine restore_photons(file,pgrid)

    character(2000),intent(in)                                 :: file
    type(photon_current),dimension(:),allocatable, intent(out) :: pgrid
    integer(kind=4)                                            :: i, np

    ! restore photons from saved file
    open(unit=16, file=trim(file), status='old', form='unformatted', action='read')
    read(16) np

    allocate(pgrid(np))

    read(16) (pgrid(i)%ID,           i=1,np)
    read(16) (pgrid(i)%status,       i=1,np)
    read(16) (pgrid(i)%xlast(:),     i=1,np)
    read(16) (pgrid(i)%xcurr(:),     i=1,np)
    read(16) (pgrid(i)%nu_ext,       i=1,np)
    read(16) (pgrid(i)%k(:),         i=1,np)
    read(16) (pgrid(i)%nb_abs,       i=1,np)
    read(16) (pgrid(i)%time,         i=1,np)
    read(16) (pgrid(i)%tau_abs_curr, i=1,np)
    read(16) (pgrid(i)%iran,         i=1,np)

    close(16)

  end subroutine restore_photons



  subroutine save_photons(file,pgrid)

    character(2000),intent(in)                   :: file
    type(photon_current),dimension(:),intent(in) :: pgrid
    integer(kind=4)                              :: i,np

    np = size(pgrid)
    open(unit=16, file=trim(file), status='unknown', form='unformatted', action='write')
    write(16) np
    write(16) (pgrid(i)%ID,           i=1,np)
    write(16) (pgrid(i)%status,       i=1,np)
    write(16) (pgrid(i)%xlast(:),     i=1,np)
    write(16) (pgrid(i)%xcurr(:),     i=1,np)
    write(16) (pgrid(i)%nu_ext,       i=1,np)
    write(16) (pgrid(i)%k(:),         i=1,np)
    write(16) (pgrid(i)%nb_abs,       i=1,np)
    write(16) (pgrid(i)%time,         i=1,np)
    write(16) (pgrid(i)%tau_abs_curr, i=1,np)
    write(16) (pgrid(i)%iran,         i=1,np)
    close(16)

  end subroutine save_photons



  subroutine dump_photons(file,pgrid)

    character(2000),intent(in)                   :: file
    type(photon_current),dimension(:),intent(in) :: pgrid
    integer(kind=4)                              :: i,np

    np = size(pgrid)
    open(unit=14, file=trim(file), status='unknown', form='unformatted', action='write')
    write(14) np
    write(14) (pgrid(i)%ID,      i=1,np)
    write(14) (pgrid(i)%status,  i=1,np)
    write(14) (pgrid(i)%xlast(:),i=1,np)
    write(14) (pgrid(i)%nu_ext,  i=1,np)
    write(14) (pgrid(i)%k(:),    i=1,np)
    write(14) (pgrid(i)%nb_abs,  i=1,np)
    write(14) (pgrid(i)%time,    i=1,np)
!==================================================
    write(14) (pgrid(i)%total_OD,i=1,np)
!==================================================
    close(14)

  end subroutine dump_photons

end module module_photon
