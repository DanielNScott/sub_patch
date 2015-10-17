module ed_state_vars
implicit none
!=============================================================================================!
! Types
!=============================================================================================!
!---------------------------------------------------------------------------------------------!
type varpropstype
    !-----------------------------------------------------------------------------------------!
    ! This type exists to centralize data controlling how subroutines interact with different !
    ! variables based on their types and properties. Could this be done using polymorphism?   !
    ! Maybe not; Since fortran doesn't support multiple inheritance, the number of types would!
    ! have to grow combinatorially with the number of properties. Or maybe it could be, if a  !
    ! tree structure was used for typing?                                                     !
    !                                                                                         !
    ! In any case, the point is to have a lookup table specifying how variables should be     !
    ! manipulated, which is then interpreted across the model. Note that this already         !
    ! essentially exists in the form of type var_table in module ed_var_tables! As such, my   !
    ! proposal here is really to rework/add to it, and to start using it outside i/o.         !
    !-----------------------------------------------------------------------------------------!
    character(len=16)                   :: name
    logical                             :: is_avg
    logical                             :: is_cbvar
    logical                             :: is_dmean
    logical                             :: is_mmean
    logical                             :: is_C13
    logical                             :: scale_by_lai
    logical                             :: scale_by_nplant

    real                                :: init_val
    integer                             :: val_dim
    integer                             :: val_dim_len1
    integer                             :: val_dim_len2

end type varpropstype
!---------------------------------------------------------------------------------------------!



!---------------------------------------------------------------------------------------------!
type subpatype
    !-----------------------------------------------------------------------------------------!
    ! The sub-patch type is my solution to making the patch type into a traversable data      !
    ! structure. This has the potential to make big hunks of code MUCH smaller, and to make   !
    ! aspects of the model more conceptually clear, as described elsewhere. It contains the   !
    ! actual values of variables in the patch, as well as info about what to do with them,    !
    ! which should ensure their proper use/propagation.                                       !
    !                                                                                         !
    ! Note that the subpa is not what will be passed around to most routines; The patch type  !
    ! will still interface with most of ED, and as a result there shouldn't be a performance  !
    ! penalty incurred by adopting it.                                                        !
    !-----------------------------------------------------------------------------------------!
    real, pointer, dimension(:)     :: val_1d           ! dimension(  ncohorts), e.g. bleaf
    real, pointer, dimension(:,:)   :: val_2d           ! dimension(:,ncohorts), e.g. cb
    logical                         :: is_avg           ! control variable
    logical                         :: is_cbvar         ! control variable
    logical                         :: is_dmean         ! control variable
    logical                         :: is_C13           ! control variable
    logical                         :: scale_by_lai     ! control variable
    logical                         :: scale_by_nplant  ! control variable
    character(len=16)               :: name             ! Name in patch pointing to 'val'

    contains
        !-------------------------------------------------------------------------------------!
        ! These proceedures implement the cores of various functions of fuse_fiss_utils.      !
        ! 'Aliases' below refer to the subroutine names which these are used in.              !
        !                                                                                     !
        ! Why implement the core of these subroutines as type-bound procedures?               !
        !   - They can't be screwed up when you create a new variable                         !
        !   - It makes iterating over the sub-patch to apply them more intuitive              !
        !   - It makes it conceptually clear what is happening - patch variables need THIS    !
        !     (below) list of things to happen to them. From a design perspective, this is    !
        !     great, because I know exactly what's going on with this data.
        !-------------------------------------------------------------------------------------!
        procedure   :: create_clone ! Alias: clone_cohort
        procedure   :: fuse         ! Alias: fuse_2_cohorts
!       procedure   :: rescale      ! Alias: split_cohorts, rescale_patches
        !-------------------------------------------------------------------------------------!
        ! NOTE: All of these proceedures could be wrappers for a single proceedure, call it   !
        ! recombine(co1,co2,sc1,sc2,sc3) which scales cohort 1 by sc1, cohort2 by sc2, and    !
        ! their sum by sc3, and assigns this to cohort 1.                                     !
        !-------------------------------------------------------------------------------------!

        !-------------------------------------------------------------------------------------!
        ! Some implementations for average_utils                                              !
        !-------------------------------------------------------------------------------------!
        !procedure   :: zero_mo_vars ! Alias: zero_ed_monthly_output_vars


end type subpatype
!---------------------------------------------------------------------------------------------!



!---------------------------------------------------------------------------------------------!
type patchtype
    integer                       :: nvars = 3     ! number of patch vars other than this one.
    real, pointer, dimension(:)   :: bleaf
    real, pointer, dimension(:)   :: broot
    real, pointer, dimension(:,:) :: cb
end type patchtype
!---------------------------------------------------------------------------------------------!



!---------------------------------------------------------------------------------------------!
type sitetype

    type(patchtype), pointer    :: patch

end type sitetype
!---------------------------------------------------------------------------------------------!

!=============================================================================================!
! Module Namespace
!=============================================================================================!
type(varpropstype), allocatable, dimension(:)    :: varprops

contains
!=============================================================================================!
! Type Bound Proceedures
!=============================================================================================!


!---------------------------------------------------------------------------------------------!
    subroutine create_clone(subpa,isc,idt)
        implicit none
        !--- Arguments -----------------------------------------------------------------------!
        class(subpatype)             :: subpa   ! Sub-Patch
        integer                      :: isc     ! Index of "Source" cohort
        integer                      :: idt     ! Index of "Destination" cohort"
        !--- Local Vars ----------------------------------------------------------------------!
        integer                      :: imonth
        !-------------------------------------------------------------------------------------!

        if (subpa%is_cbvar) then
            do imonth = 1,3
                subpa%val_2d(imonth,idt) = subpa%val_2d(imonth,isc)
            end do
        else
            subpa%val_1d(idt) = subpa%val_1d(isc)
        end if

    end subroutine create_clone
!---------------------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------------------!
    subroutine fuse(subpa,donc,recc,nplant,lai,newn)
        implicit none
        !--- Arguments -----------------------------------------------------------------------!
        ! Note: Some of these are real in fuse_fiss, but seem like they should be ints...
        !-------------------------------------------------------------------------------------!
        class(subpatype)             :: subpa   ! Sub-Patch
        integer                      :: donc    ! Donating cohort.
        integer                      :: recc    ! Receptor cohort.
        real, dimension(:)           :: nplant  ! 'Current' nplant
        real, dimension(:)           :: lai     ! 'Current' nplant
        real                         :: newn    ! New nplant
        !--- Local Vars ----------------------------------------------------------------------!
        integer                      :: imonth
        real                         :: newni
        real                         :: newlaii
        !-------------------------------------------------------------------------------------!
        !------------------------------------------------------------------------------------!
        !    Find the scaling factor for variables that are not "extensive".                 !
        !  - If the unit is X/plant, then we scale by nplant.                                !
        !  - If the unit is X/m2_leaf, then we scale by LAI.                                 !
        !  - If the unit is X/m2_gnd, then we add, since they are "extensive".               !
        !------------------------------------------------------------------------------------!
        newni   = 1.0 / newn
        if (lai(recc) + lai(donc) > 0.0) then
         newlaii = 1.0 / (lai(recc) + lai(donc))
        else
         newlaii = 0.0
        end if
        !------------------------------------------------------------------------------------!

        write(*,*) subpa%name, subpa%is_cbvar

        if (subpa%is_cbvar) then
            do imonth = 1,3
                subpa%val_2d(imonth,recc) = (nplant(recc) *subpa%val_2d(imonth,recc)          &
                                           + nplant(donc) *subpa%val_2d(imonth,donc) ) * newni
            end do
        else
            subpa%val_1d(recc) = (nplant(recc) *subpa%val_1d(recc)                            &
                                + nplant(donc) *subpa%val_1d(donc) ) * newni
        end if

    end subroutine fuse
!---------------------------------------------------------------------------------------------!




!=============================================================================================!
! Module Subroutines
!=============================================================================================!
!---------------------------------------------------------------------------------------------!
! THIS SUBROUTINE, assoc_subpa should be one of only three places variables need to be        !
! manually added to the model (to achieve basic functionality), along with addition to the    !
! type construct, and the varprop assigment. It is really just a casing wrapper for the SR    !
! assoc_null() found below.                                                                   !
!---------------------------------------------------------------------------------------------!
    subroutine assoc_subpa(subpa,patch,assoc)
        implicit none
        type(subpatype), pointer, dimension(:)  :: subpa
        type(patchtype)                         :: patch
        logical                                 :: assoc
        integer                                 :: i

        do i = 1,size(subpa)
            select case(subpa(i)%name)
            case('bleaf'); call assoc_null(subpa(i),assoc, d1var = patch%bleaf          )
            case('broot'); call assoc_null(subpa(i),assoc, d1var = patch%broot          )
            case('cb'   ); call assoc_null(subpa(i),assoc, d2var = patch%cb             )
            case default
                write (*,*) 'Error: Nothing associated in subpatype...'
            end select
        end do

        !-------------------------------------------------------------------------------------!
        ! Deallocating subpatch can be done here if there is no instance in which we want to  !
        ! keep it around while dropping patchtype pointers, otherwise can be done elsewhere   !
        !-------------------------------------------------------------------------------------!
         if ( .not. assoc) then
             deallocate(subpa)
         end if

    end subroutine assoc_subpa
!---------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------------!
! Associates individual patch variables with sub-patch values.
!---------------------------------------------------------------------------------------------!
    subroutine assoc_null(subpa,assoc,d1var,d2var)
        implicit none
        type(subpatype), target                 :: subpa     ! subpatch variable
        logical                                 :: assoc     ! Associating, or nullifying?
        real, pointer, optional, dimension(:)   :: d1var     ! Dimension 1 patch variable
        real, pointer, optional, dimension(:,:) :: d2var     ! Dimension 2 patch variable


        if (assoc) then
            if (present(d1var)) then
                d1var => subpa%val_1d
            else if (present(d2var)) then
                d2var => subpa%val_2d
            end if
        else
            if (present(d1var)) then
                nullify(d1var)
            else if (present(d2var)) then
                nullify(d2var)
            end if
        end if

    end subroutine assoc_null
!---------------------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------------------!
    subroutine init_subpa(subpa,subpa_dim)
        implicit none
        type(subpatype),pointer, dimension(:)   :: subpa
        integer                                 :: i
        integer                                 :: j
        integer                                 :: k
        integer                                 :: subpa_dim

        allocate(subpa(subpa_dim))
        do i = 1,subpa_dim
            !--------------------------------------------------------------------------------!
            ! Allocate dimensions of sub patch type depending on the variable's spec. in     !
            ! the 'varprop' variable.                                                        !
            !--------------------------------------------------------------------------------!
            subpa(i)%name     = trim(varprops(i)%name)
            subpa(i)%is_cbvar = varprops(i)%is_cbvar

            select case(varprops(i)%val_dim)
            case(1)
                allocate(subpa(i)%val_1d(varprops(i)%val_dim_len1))

                do j = 1,varprops(i)%val_dim_len1
                    subpa(i)%val_1d(j) = varprops(i)%init_val
                end do
            case(2)
                allocate(subpa(i)%val_2d(varprops(i)%val_dim_len1,varprops(i)%val_dim_len2))

                do k = 1,varprops(i)%val_dim_len1
                    do j = 1,varprops(i)%val_dim_len2
                        subpa(i)%val_2d(k,j) = varprops(i)%init_val
                    end do
                end do
            end select
        end do

    end subroutine init_subpa
!---------------------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------------------!
    subroutine print_subpa(subpa)
        implicit none
        type(subpatype), dimension(:)       :: subpa
        integer                             :: i
        integer                             :: j

        do i = 1,size(subpa)
            if (associated(subpa(i)%val_1d)) then
                write(*,*) 'subpa(i)%name, subpa(i)%val_1d(:)    : ', subpa(i)%name, subpa(i)%val_1d(:)
            else if (associated(subpa(i)%val_2d)) then
                do j = 1,size(subpa(i)%val_2d(:,1))
                    if (j == 1) then
                        write(*,*) 'subpa(i)%name, subpa(i)%val_2d(:,:)  : ', subpa(i)%name,  subpa(i)%val_2d(j,:)
                    else
                        write(*,*) '                                     : ', subpa(i)%name,  subpa(i)%val_2d(j,:)
                    end if
                end do
            end if
        end do

    end subroutine print_subpa
!---------------------------------------------------------------------------------------------!



end module ed_state_vars
