subroutine doubleSum(ket, bra, M , states, tessere,prod)
    integer :: its, j, k
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere,states,states),intent(in) :: M
    complex*16, dimension(states),intent(in) :: ket
    complex*16, dimension(states),intent(in) :: bra
    complex*16, dimension(tessere),intent(out) :: prod
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
        prod(its)=dot_product(bra,ctmp(its,:))
    enddo
!    write(*,*) prod(30)
end subroutine


subroutine single_summation_tessere(ket, M , states, tessere,prod)
    integer :: its, j, k
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere,states,states),intent(in) :: M
    complex*16, dimension(tessere),intent(in) :: ket
    complex*16, dimension(states, states),intent(out) :: prod
    prod=(0.d0, 0d0)
    do k=1,states
         do j=1,states
            do its=1,tessere
               prod(k,j)=prod(k,j)+ M(its,k,j)*ket(its)
            enddo
         enddo
    enddo
end subroutine



subroutine propagate_q_frozen(psi_prev, qijn, states, tessere, out_prop_charges)
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere, states, states),intent(in) :: qijn
    complex*16, dimension(states),intent(in) :: psi_prev
    complex*16, dimension(tessere),intent(out) :: out_prop_charges

    call doublesum(psi_prev, psi_prev, qijn, states, tessere, out_prop_charges)
end subroutine



subroutine eulero_pcm(psi_prev, qijn_t, Vijn, states, tessere, out_matmul_sum)
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere),intent(in) :: qijn_t
    complex*16, dimension(tessere, states, states),intent(in) :: Vijn
    complex*16, dimension(states),intent(in) :: psi_prev
    complex*16, dimension(states),intent(out) :: out_matmul_sum

    complex*16, dimension(states,states) :: out_single_summation_tessere

    call single_summation_tessere(qijn_t, Vijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = matmul(psi_prev, out_single_summation_tessere)



end subroutine


subroutine bwd_pcm(chi_prev, wf_fwd, qijn_t, Vijn, qijn, states, tessere, out_matmul_sum)
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere),intent(in) :: qijn_t
    complex*16, dimension(tessere, states, states),intent(in) :: qijn
    complex*16, dimension(tessere, states, states),intent(in) :: Vijn
    complex*16, dimension(states),intent(in) :: chi_prev
    complex*16, dimension(states),intent(in) :: wf_fwd
    complex*16, dimension(states),intent(out) :: out_matmul_sum


    complex*16, dimension(tessere) :: out_doublesum
    complex*16, dimension(states,states) :: out_single_summation_tessere


    call single_summation_tessere(qijn_t, Vijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = matmul(chi_prev, out_single_summation_tessere)

    call doublesum(chi_prev, wf_fwd, Vijn, states, tessere, out_doublesum)
    call single_summation_tessere(out_doublesum, qijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = out_matmul_sum + matmul(wf_fwd, out_single_summation_tessere)

    call doublesum(wf_fwd, chi_prev, qijn, states, tessere, out_doublesum)
    call single_summation_tessere(out_doublesum, Vijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = out_matmul_sum - matmul(wf_fwd, out_single_summation_tessere)
end subroutine




subroutine bwd_pcm_vecchia(chi_prev, wf_fwd, qijn, Vijn, subtract, states, tessere, out_matmul_sum)
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere, states, states),intent(in) :: qijn
    complex*16, dimension(tessere, states, states),intent(in) :: Vijn
    complex*16, dimension(states),intent(in) :: chi_prev
    complex*16, dimension(states),intent(in) :: wf_fwd
    complex*16, dimension(states),intent(out) :: out_matmul_sum
    complex*16, dimension(states, states),intent(in) :: subtract

    complex*16, dimension(tessere) :: out_doublesum
    complex*16, dimension(states,states) :: out_single_summation_tessere


    call doublesum(wf_fwd, wf_fwd, qijn, states, tessere, out_doublesum)
    call single_summation_tessere(out_doublesum, Vijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = matmul(chi_prev, out_single_summation_tessere - subtract)

    call doublesum(chi_prev, wf_fwd, Vijn, states, tessere, out_doublesum)
    call single_summation_tessere(out_doublesum, qijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = out_matmul_sum + matmul(wf_fwd, out_single_summation_tessere)

    call doublesum(wf_fwd, chi_prev, qijn, states, tessere, out_doublesum)
    call single_summation_tessere(out_doublesum, Vijn, states, tessere, out_single_summation_tessere)
    out_matmul_sum = out_matmul_sum - matmul(wf_fwd, out_single_summation_tessere)
end subroutine



subroutine fwd_pcm(psi_prev, qijn, Vijn, subtract, states, tessere, out_matmul_sum)
    integer :: states
    integer :: tessere
    complex*16, dimension(tessere, states, states),intent(in) :: qijn
    complex*16, dimension(tessere, states, states),intent(in) :: Vijn
    complex*16, dimension(states),intent(in) :: psi_prev
    complex*16, dimension(states),intent(out) :: out_matmul_sum
    complex*16, dimension(states, states),intent(in) :: subtract

    complex*16, dimension(tessere) :: out_doublesum
    complex*16, dimension(states,states) :: out_single_summation_tessere
    complex*16, dimension(states,states) :: tmp 

    call doublesum(psi_prev, psi_prev, qijn, states, tessere, out_doublesum)
    call single_summation_tessere(out_doublesum, Vijn, states, tessere, out_single_summation_tessere)
    tmp=out_single_summation_tessere - subtract
    out_matmul_sum = matmul(psi_prev, tmp)
end subroutine

