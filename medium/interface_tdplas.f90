module interface_tdplas
      use tdplas
      implicit none




      public init_global_tdplas, call_prop_charges, call_init_charges  

      contains


            subroutine init_global_tdplas(dt, vts, n_tessere, n_states)

                real*8     , intent(in)  :: dt                          ! time step
                integer    , intent(in)  :: n_states, n_tessere         ! number of CIS states
                real*8     , intent(in)  :: vts(:,:,:)

                call  quantum_init(dt, vts, n_tessere, n_states)
                call  readio_and_init_tdplas_for_ocpy
            end subroutine


            subroutine call_init_charges(ci, field_vector, n_states)
                complex*16, intent(in) :: ci(n_states)
                real*8, dimension(3), intent(in) :: field_vector
                integer :: n_states
                call init_charges_int(ci, field_vector, quantum_n_ci, pedra_surf_n_tessere)
            end subroutine


            subroutine init_charges_int(ci, field_vector, n_states, n_tessere)
                complex*16, intent(in) :: ci(n_states)
                real*8, dimension(3), intent(in) :: field_vector
                real*8, allocatable  :: V_reactionf(:)
                real*8, allocatable :: V_localf(:)
                integer, intent(in) :: n_tessere, n_states


                allocate(V_reactionf(n_tessere))
                allocate(V_localf(n_tessere))

                call deallocate_BEM_public
                call deallocate_potential
                call finalize_prop

                call prepare_potentials(ci, field_vector, V_reactionf, V_localf, quantum_n_ci, n_tessere)

                call do_BEM_prop
                call init_potential(V_reactionf, V_localf)
                call init_charges
                call init_vv_propagator
                return
            end subroutine





            subroutine call_prop_charges(ci, field_vector, q_t, n_states, n_tessere)
                complex*16, intent(in) :: ci(n_states)
                real*8, dimension(3), intent(in) :: field_vector
                integer :: n_states, n_tessere
                real*8, intent(inout) :: q_t(n_tessere)
                call prop_charges(ci, field_vector, q_t, quantum_n_ci, pedra_surf_n_tessere)
            end subroutine



           subroutine prop_charges(ci, field_vector, q_t, n_states, n_tessere)
               complex*16, intent(in) :: ci(n_states)
               real*8, dimension(3), intent(in) :: field_vector
               real*8, intent(out) :: q_t(n_tessere)
               real*8, allocatable :: V_reactionf(:)
               real*8, allocatable :: V_localf(:)
               integer, intent(in) :: n_states, n_tessere



               allocate(V_reactionf(n_tessere))
               allocate(V_localf(n_tessere))
               call prepare_potentials(ci, field_vector, V_reactionf, V_localf,quantum_n_ci, n_tessere)
               call set_potential(V_reactionf, V_localf)
               call prop_chr
               call get_corrected_propagated_charges(q_t)

           end subroutine


            subroutine prepare_potentials(ci, field_vector, V_reactionf, V_localf, n_states, n_tessere)
                complex*16, dimension(n_states), intent(in) :: ci 
                real*8, dimension(3), intent(in) :: field_vector
                real*8, dimension(n_tessere), intent(inout) :: V_reactionf(:)
                real*8, dimension(n_tessere), intent(inout) :: V_localf(:)
                integer,intent(in)  :: n_states, n_tessere 
                integer :: i

                call doubleSum(ci, conjg(ci), quantum_vts, n_states, n_tessere, V_reactionf)


                do i=1, pedra_surf_n_tessere
                    V_localf(i) = -1.*(field_vector(1)*pedra_surf_tessere(i)%x + &
                                  field_vector(2)*pedra_surf_tessere(i)%y + &
                                  field_vector(3)*pedra_surf_tessere(i)%z)
                enddo
            end subroutine
!
!
            subroutine doubleSum(ket, bra, M , states, tessere,prod)
                integer :: its, j, k
                integer :: states
                integer :: tessere
                real*8, dimension(tessere,states,states),intent(in) :: M
                complex*16, dimension(states),intent(in) :: ket
                complex*16, dimension(states),intent(in) :: bra
                complex*16, dimension(tessere) :: ptmp
                real*8, dimension(tessere),intent(out) :: prod
                complex*16,dimension(tessere,states) :: ctmp
                ctmp=(0.d0, 0d0)
                do k=1,states
                    do j=1,states
                        do its=1,tessere
                            ctmp(its,k)=ctmp(its,k)+ M(its,k,j)*ket(j)
                        enddo
                    enddo
                enddo
                do its=1,tessere
                    ptmp(its)=dot_product(bra,ctmp(its,:))
                enddo
                prod = real(ptmp)
            end subroutine




        end module

