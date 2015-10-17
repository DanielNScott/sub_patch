!---------------------------------------------------------------------------------------------!
! This program is designed to provide a proof of concept for the idea that ed                 !
! variables can be reworked so that all of the existing code functions properly AND meets the !
! goals that...                                                                               !
!                                                                                             !
! - adding variables should be conceptually straight-forward and require a minimum of coding  !
!   to ensure proper functionality in all book-keeping respects. (e.g. if a var is a patch    !
!   variable it should automatically be manipulated in fussion/fission)                       !
!                                                                                             !
! - we should be capable of iterating over ed structures, which would greatly reduce code     !
!   size, while improving interpretability, maintainability, and readability                  !
!                                                                                             !
! - a variable should contain information about it's use and manipulation, obviating the need !
!   for a developer to (correctly) make those choices throughout aspects of the model which   !
!   are not actually of scientific or project-related interest.                               !
!---------------------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------------------!
! The program 'main' is a series of tests/implementations of the toy version of certain ed    !
! subcomponents built in ed_state_vars.
!---------------------------------------------------------------------------------------------!
program main

    use ed_state_vars
    implicit none

    !--- General Local Vars ------------------------------------------------------------------!
    type(subpatype), pointer, dimension(:)  :: subpa            ! Sub-patch
    type(patchtype)                         :: patch            ! Local toy patch structure
    logical                                 :: assoc = .true.   ! Flag for init_subpa
    integer                                 :: i                ! Loop index
    integer                                 :: j                ! Loop index
    integer                                 :: ncohorts = 3     ! Number of cohorts
    !--- Local Fusion/Fission Vars -----------------------------------------------------------!
    integer                         :: donc     ! Donating cohort. (std ED name)
    integer                         :: recc     ! Receptor cohort. (std ED name)
    real, allocatable, dimension(:) :: nplant   ! 'Current' nplant (std ED name)
    real, allocatable, dimension(:) :: lai      ! 'Current' nplant (std ED name)
    real                            :: newn     ! New nplant       (std ED name)
    !-----------------------------------------------------------------------------------------!
    ! NOTES:                                                                                  !
    !   - SUBPATYPE is defined in ed_state_vars. It is my solution to the problem of making   !
    !     the 'patchtype' data structure that currently lives in ED amenable to being         !
    !     iterated over.                                                                      !
    !                                                                                         !
    !   - PATCHTYPE is a mock-up of the standard patchtype in ED                              !
    !                                                                                         !
    !   - ASSOC is a flag to tell init_subpa if it should associate or nullify patch pointers !
    !     to the sub-patch data structure.                                                    !
    !-----------------------------------------------------------------------------------------!


    !-----------------------------------------------------------------------------------------!
    ! Allocate and set varprops variable which localizes information about how a variable     !
    ! gets manipulated. It can be deallocated when not required. This particular              !
    ! implenentation is clumsy but illustrates the point. A character matrix that gets        !
    ! tokenized would probably be the most intuitive, easy to manipulate approach.            !
    !-----------------------------------------------------------------------------------------!
    allocate(varprops(patch%nvars))

    varprops(:)%name         = ['bleaf','broot','cb   ']
    varprops(:)%init_val     = 0
    varprops(:)%is_cbvar     = [.false., .false., .true.]
    varprops(:)%val_dim      = [1,1,2]
    varprops(:)%val_dim_len1 = ncohorts
    varprops(3)%val_dim_len1 = ncohorts
    varprops(3)%val_dim_len2 = 4

    !-----------------------------------------------------------------------------------------!
    ! Initialize subpa and assign values to it, then print it.                              !
    !-----------------------------------------------------------------------------------------!
    write(*,*) ''
    write(*,*) '-------------------------------------------------------'
    write(*,*) '        Testing init_subpa and print_subpa.'
    write(*,*) 'Iterating through subpa, allocating, valuing, printing'
    write(*,*) '-------------------------------------------------------'

    call init_subpa(subpa, patch%nvars)

    do i = 1,patch%nvars
        if ( varprops(i)%val_dim == 1 ) then
            subpa(i)%val_1d(1) = 15.0 - i
        else
            subpa(i)%val_2d(:,1) = i
            subpa(i)%val_2d(:,2) = [1,2,3]
        end if
    end do

    call print_subpa(subpa)

    !-----------------------------------------------------------------------------------------!
    ! Test assoc_subpa fuction which associates a patchtype variable with a sub-patch. This   !
    ! way of associating a patch compresses all the allocate, deallocate, nullify code at the !
    ! patch level by a factor of 3 and, with the aid of varprops for centrally storing        !
    ! information about variables makes organizing allocations by conditionals unneccessary.  !
    !-----------------------------------------------------------------------------------------!
    write(*,*) ''
    write(*,*) '----------------------------------------------------'
    write(*,*) '    Testing assoc_subpa association function'
    write(*,*) 'Assigning patch%bleaf => subpa(-)%val variables'
    write(*,*) '----------------------------------------------------'
    call assoc_subpa(subpa,patch,assoc)

    write(*,*) 'patch%bleaf       : ', patch%bleaf
    write(*,*) 'patch%broot       : ', patch%broot
    write(*,*) 'patch%cb          : ', patch%cb(1,:)
    write(*,*) '                  : ', patch%cb(2,:)
    write(*,*) '                  : ', patch%cb(3,:)
    write(*,*) ''

    do i=1,patch%nvars
        if (varprops(i)%val_dim == 1) then
            write(*,*) 'subpa(i)%name, subpa(i)%val_1d     : ', subpa(i)%name, subpa(i)%val_1d
        
        else if (varprops(i)%val_dim == 2) then
             write(*,*) 'subpa(i)%name                      : ', subpa(i)%name

            do j=1,3
                 write(*,*) 'subpa(i)%val_2d                   : ', subpa(i)%val_2d(j,:)
            end do
        end if
    end do

    !-----------------------------------------------------------------------------------------!
    ! Test cohort assignments with patch vars are equiv to assigments at sub-patch, just      !
    !  to make sure everything's work properly...                                             !
    !-----------------------------------------------------------------------------------------!
    write(*,*) ''
    write(*,*) '----------------------------------------------------'
    write(*,*) '    Testing cohort assigment via patch vars'
    write(*,*) '----------------------------------------------------'

    patch%bleaf    (2:3)     = [1, 2]
    patch%broot    (2:3)     = [3, 4]
    call print_subpa(subpa)


    !-----------------------------------------------------------------------------------------!
    ! Test sub-patch level implementation of cohort fusion.
    !-----------------------------------------------------------------------------------------!
    write(*,*) ''
    write(*,*) '-------------------------------------------------------'
    write(*,*) '    Testing cohort fusion loop'
    write(*,*) 'This performs (nplant1 *val1 + nplant2 *val2) *newni'
    write(*,*) '-------------------------------------------------------'

    allocate(lai(ncohorts))
    allocate(nplant(ncohorts))

    donc    = 2
    recc    = 3
    nplant  = [2, 1, 1]
    lai     = [0.387, 0.563, 0.01]
    newn    = 2

    Write(*,*) 'donc    :', donc
    Write(*,*) 'recc    :', recc
    Write(*,*) 'nplant  :', nplant(:)
    Write(*,*) 'lai     :', lai(:)
    Write(*,*) 'newn    :', newn

    do i=1,patch%nvars
        call subpa(i)%fuse(donc,recc,nplant,lai,newn)
    end do
    call print_subpa(subpa)


    !-----------------------------------------------------------------------------------------!
    ! Test assoc_subpa nullify function
    !-----------------------------------------------------------------------------------------!
    write(*,*) ''
    write(*,*) '----------------------------------------------------'
    write(*,*) '    Testing assoc_subpa nullify function'
    write(*,*) ' For nullifying the patch var names, which should be'
    write(*,*) ' nullified before nullifying (via looping) the subpa'
    write(*,*) '----------------------------------------------------'
    assoc = .false.
    call assoc_subpa(subpa,patch,assoc)

    if (associated(patch%bleaf    )) write(*,*) 'patch%bleaf       : ', patch%bleaf
    if (associated(patch%broot    )) write(*,*) 'patch%broot       : ', patch%broot

    if (.not. associated(patch%bleaf    )) write(*,*) 'patch%bleaf       : null'
    if (.not. associated(patch%broot    )) write(*,*) 'patch%broot       : null'

    if (associated(subpa)) then
        do i=1,patch%nvars
            write(*,*) 'i, subpa(i)%name, subpa(i)%val_1d     : ', i, subpa(i)%name, subpa(i)%val_1d
        end do
    else
        write(*,*) 'subpa             : null'
    end if


end program main
