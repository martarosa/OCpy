module interface_ocpy

    use constants
    use readio_medium
    use td_contmed

    implicit none


      real(dbl) 		        :: dt			        	! time step
      character(flg) 	        :: Fmdm				        ! kind of medium
      integer(i4b)       	    :: n_ci,n_ci_read	     	! number of CIS states
      real(dbl), allocatable    :: c_i(:),e_ci(:)	    	! CIS coefficients and energies
      real(dbl), allocatable    :: mut(:,:,:)	        	! CIS transition dipoles
      real(dbl) 		        :: mol_cc(3)			    ! molecule center
      real(dbl) 		        :: fmax(3,10),omega(10)		! field amplitude and frequency
      real(dbl)                 :: tdelay(10),pshift(10)    ! time delay and phase shift
      character(flg) 	        :: Ffld				        ! shape of impulse
      character(flg)            :: Fbin                     ! binary output
      character(flg)            :: Fopt                     ! OMP-optimized matrix/vector multiplication
      integer(i4b) 	        :: n_out,n_f		     	! auxiliaries for output
      integer(i4b)              :: nthr                     ! number of threads
      character(1)              :: medium_res              !restart for medium
      integer(i4b) :: n_res ! frequency for restart


    contains


    subroutine set_global_tdplas(oc_dt, oc_env, oc_ci, oc_en_ci, oc_mut, oc_n_ci)


        implicit none


        real(dbl)     , intent(in)           :: oc_dt				                    ! time step
        character(3)  , intent(in)           :: oc_env        	                        ! kind of medium

        integer(i4b)  , optional, intent(in) :: oc_n_ci

        complex(cmp)  , intent(in) :: oc_ci(:)                    ! CIS coefficients
        real(dbl)     , dimension(oc_n_ci), intent(in)             :: oc_en_ci 	        ! CIS energies
        real(dbl)     , dimension(3, oc_n_ci, oc_n_ci), intent(in) :: oc_mut(:,:,:)	    ! CIS transition dipoles



        dt = oc_dt
        Fmdm = oc_env
        n_ci = oc_n_ci
        c_i = oc_ci
        en_ci = oc_en_ci
        mut = oc_mut

        Fbin = 'non'

        n_ci_read = n_ci
        mol_cc = 0
        fmax = 0
        omega = 0
        tdelay = 0
        pshift = 0
        Ffld = 'non'
        Fopt = 'non'
        n_out = 1
        n_f = 1
        nthr = 1
        medium_res = 'n'
        n_res = 0

        return

    end subroutine set_global_tdplas



    subroutine read_tdplas_namelist
        call read_medium

        !MR there should be a check on the TDPLAS input file that all the variables have allowed values

    end subroutine read_rdplas_namelist

    subroutine call_init_mdm(pot_t, pot_t_lf, tessere)

        integer(i4b)  , optional, intent(in) :: tessere

        real(dbl) , dimension(tessere), intent(in) :: pot_t
        real(dbl) , dimension(tessere), intent(in) :: pot_t_lf

        call init_mdm(pot_t = pot_t, potf_t = pot_t_lf)

    end subroutine call_init_mdm

    subroutine call_prop_mdm(i, pot_t, pot_t_lf, tessere)

        integer(i4b)  , optional, intent(in) :: tessere

        real(dbl) , dimension(tessere), intent(in) :: pot_t
        real(dbl) , dimension(tessere), intent(in) :: pot_t_lf

        call init_prop(i, pot_t = pot_t, potf_t = pot_t_lf)

    end subroutine call_init_mdm


    subroutine get_q_nanop_q_t(q_t)

        real(dbl) , dimension(n_ci), intent(out) :: q_t


        q_t = qr_t + qx_t

    end subroutine call_init_mdm