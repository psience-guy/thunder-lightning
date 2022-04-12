! copyright info:
!
!                             @Copyright 2013
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! Caltech - Brandon Keith
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Program Description
! ==========================================================================
!> This is the main driver for the LIGHTNING version of FIREBALL.
! ==========================================================================
! Code written by:
!> @author Ning Ma
!! @author James P. Lewis
! Box 6315, 209 Hodges Hall
! Department of Physics
! West Virginia University
! Morgantown, WV 26506-6315
!
! (304) 293-3422 x1409 (office)
! (304) 293-5732 (FAX)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        program lightning

! /SYSTEM
        use M_species
        use M_configuraciones
        use M_neighbors
        use M_neighbors_PP
        use M_atom_functions

! /FDATA
        use M_Fdata_1c
        use M_Fdata_2c

! /GRID
        use M_grid
        use M_isosurfaces

! /ASSEMBLERS
        use M_assemble_2c
        use M_assemble_3c
        use M_assemble_ewald
        use M_assemble_usr
        use M_assemble_vxc
        use M_assemble_PP_2c
        use M_assemble_PP_3c

! /SOLVESH
        use M_kspace
        use M_density_matrix

! /SCF
        use M_charges

        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Parameters and Data Declaration
! ===========================================================================
! None

! Variable Declaration and Description
! ===========================================================================
        integer iscf_iteration
        integer istructure, iseparate

        real sigma

        character (len = 25) :: slogfile

! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        real timei, timef

! Energies
        real ebs                                     ! band-structure energy
        real uii_uee, uxcdcc                         ! short-range energies

! Interfaces
        interface
          subroutine Qmixer (t, iscf_iteration, sigma)
          use M_configuraciones
            type(T_structure), target :: t
            integer, intent (in) :: iscf_iteration
            real, intent (inout) :: sigma
          end subroutine Qmixer
        end interface

        interface
          subroutine writeout_energies (t, ebs, uii_uee, uxcdcc)
          use M_configuraciones
            type(T_structure), target :: t
            real, intent (in) :: ebs
            real, intent (in) :: uii_uee
            real, intent (in) :: uxcdcc
          end subroutine writeout_energies
        end interface

        interface
          subroutine writeout_xyz (t, ebs, uii_uee, uxcdcc)
          use M_configuraciones
            type(T_structure), target :: t
            real, intent (in) :: ebs
            real, intent (in) :: uii_uee
            real, intent (in) :: uxcdcc
          end subroutine writeout_xyz
        end interface

        interface
          subroutine absorption (t)
          use M_configuraciones
            type(T_structure), target :: t
          end subroutine absorption
        end interface

        interface
          subroutine dos (t)
          use M_configuraciones
            type(T_structure), target :: t
          end subroutine dos
        end interface

! Allocate Arrays
! ===========================================================================
! None

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call cpu_time (time_begin)
        iseparate = 1
        open (unit = ilogfile, file = 'output.log', status = 'replace')
        call welcome_lightning

! ===========================================================================
! ---------------------------------------------------------------------------
!             R E A D   I N   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
        call read_Fdata_location
        call read_info
        call read_Fdata_1c
        call read_Fdata_2c
        call read_Fdata_3c

! Read in the wavefunctions
        call read_wavefunctions

! Read parameters from structures.inp file
        write (ilogfile,*)
        write (ilogfile,*) ' Reading parameters from the structure input. '
        call read_parameters

        write (ilogfile, *)
        open (unit = 1, file = 'structures.inp', status = 'old')
        read (1, *) nstructures
        if (nstructures .gt. 999) then
          stop ' Cannot calculate more than 999 structures! '
        end if
        write (ilogfile, *) ' Number of structure calculating = ', nstructures
        allocate (structures (nstructures))

! Loop over all structures
! This loop can be made parallel if each subroutine in lightning
! is modified to take s as a parameter, rather than reference s directly
!!$omp parallel do private (istructure, s, slogfile, sigma, iscf_iteration)   &
!!$omp             private (timei, timef, ebs, uii_uee, uxcdcc, etot)
        do istructure = 1, nstructures
          s => structures(istructure)
          read (1, *) s%basisfile
          if (iseparate .eq. 1) then
            s%logfile = istructure + 1000
            s%inpfile = istructure + 2000
            slogfile = s%basisfile(:len(trim(s%basisfile)) - 4)
            slogfile = trim(slogfile)//'.log'
            open (unit = s%logfile, file = slogfile, status = 'replace')
          end if
          write (ilogfile,*)
          write (ilogfile, 100) slogfile
          write (s%logfile, *) ' Structure = ', istructure

          ! Read in the coordinates and parameters
          call read_positions (s)

          ! Set the charges
          call read_charges (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!           N E I G H B O R S   S Y S T E M   I N F O R M A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
          call driver_neighbors (s)
          call driver_neighbors_PP (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!              A S S E M B L E   T H E   H A M I L T O N I A N
! ---------------------------------------------------------------------------
! ===========================================================================
! Assemble the Hamiltonian matrix:
          write (s%logfile, *)
          write (s%logfile, *) ' Calling two-center non-charge dependent assemblers. '
          call assemble_S (s)
          call assemble_T (s)
          call assemble_dipole_z (s)
          call assemble_svnl (s)
          call assemble_vnl_2c (s)

          write (s%logfile,*) ' Calling three-center non-charge dependent assemblers. '
          call assemble_vnl_3c (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                  G R I D    I N I T I A L I Z A T I O N
! ---------------------------------------------------------------------------
! ===========================================================================
          if (iwriteout_ewf .ne. 0) then
            write (s%logfile, *)
            write (s%logfile, *) ' Initialize the grid. '
            call initialize_grid (s)
          end if

! ===========================================================================
! ---------------------------------------------------------------------------
!              S C F   L O O P
! ---------------------------------------------------------------------------
! Note that self-consistency is preformed regardless of method used.
! But, in Harris, we just do a single scf cycle.
! ===========================================================================
          sigma = 999.0d0
          iscf_iteration = 1
          do while (sigma .gt. scf_tolerance_set .and.                       &
      &             iscf_iteration .le. max_scf_iterations_set)
            write (s%logfile, *)
            write (s%logfile, *) ' Begin scf step = ', iscf_iteration
            write (s%logfile, *) ' Calling two-center charge dependent assemblers. '
            call assemble_vna_2c (s)
            call assemble_ewaldsr (s)
            call assemble_ewaldlr (s)

            write (s%logfile, *) ' Calling three-center charge dependent assemblers. '
            call assemble_vna_3c (s)
            call assemble_vxc (s)

! ===========================================================================
! ---------------------------------------------------------------------------
!                         D I A G O N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Calculating the overlap matrix in K space
            write (s%logfile, *) ' Calling kspace: '
            call cpu_time (timei)
            call driver_kspace (s, iscf_iteration)
            call cpu_time (timef)
            write (s%logfile, *)
            write (s%logfile, *) ' kspace time: ', timef - timei
            call density_matrix (s)
            if (iwriteout_density .eq. 1) call writeout_density (s)

            call calculate_charges (s)
            if (iwriteout_charges .eq. 1) call writeout_charges (s)
            call Qmixer (s, iscf_iteration, sigma)

! ===========================================================================
! ---------------------------------------------------------------------------
!                       T O T A L   E N E R G I E S
! ---------------------------------------------------------------------------
! ===========================================================================
            call calculate_ebs (s, ebs)
            call assemble_uee (s, uii_uee)
            call assemble_uxc (s, uxcdcc)
            call writeout_energies (s, ebs, uii_uee, uxcdcc)

! End scf loop
            if (sigma .gt. 0.0d0) then
              iscf_iteration = iscf_iteration + 1
            else
              exit
            end if
          end do

          if (iwriteout_populations .eq. 1) call calculate_populations (s)
          if (iwriteout_xyz .eq. 1) call writeout_xyz (s, ebs, uii_uee, uxcdcc)
          if (iwriteout_ewf .ne. 0) call project_orbitals_grid (s, 1)

! Destroy Hamiltonian matrix elements storage
          call destroy_assemble_2c (s)
          call destroy_assemble_PP_2c (s)
          call destroy_assemble_vxc_McWEDA (s)

          if (iwriteout_ewf .ne. 0) call destroy_grid

! Calculate the absorption spectrum.
          if (iwriteout_abs .eq. 1) call absorption (s)

! Destroy final arrays
          call destroy_kspace (s)
          call destroy_denmat (s)
          call destroy_charges (s)

          call destroy_neighbors (s)
          call destroy_neighbors_PP (s)

          ! destroy neighbors last
          deallocate (s%xl) ! where to put this?

! Calculate the electronic density of states.
! We do this after destorying some arrays so that we can optimize the
! memory usage.
          if (iwriteout_dos .eq. 1) call dos (s)

          if (iseparate .eq. 1) close (s%logfile)
        end do ! end loop over all structures
!!$omp end parallel do
        close (1)

! ===========================================================================
! ---------------------------------------------------------------------------
!               D E S T R O Y   A R R A Y S - F I N A L I Z E
! ---------------------------------------------------------------------------
! ===========================================================================
! Destroy datafile storage
        call destroy_Fdata_1C
        call destroy_Fdata_2C
        call destroy_Fdata_3c

! Destroy SYSTEM information.
        call destroy_positions
        call destroy_species

        call cpu_time (time_end)

        write (ilogfile,*) ' LIGHTNING RUNTIME : ', time_end-time_begin, '[sec] '
        write (*,*) ' LIGHTNING RUNTIME : ', time_end-time_begin, '[sec] '
        close (ilogfile)

! Deallocate Arrays
! ===========================================================================
! None

! Format Statements
! ===========================================================================
100     format (2x, ' Working on structure - ', a25)

! End Program
! ===========================================================================
        stop
        end program lightning
