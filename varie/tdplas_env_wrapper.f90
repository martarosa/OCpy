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


subroutine single_summation_tessere(ket, M , states, tessere, prod)
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


subroutine fwd_pcm(psi_prev, Vijn, field_r_tessere, subtract, states, tessere, out_h_pcm)
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

    call doublesum(psi_prev, psi_prev, Vijn, states, tessere, out_Vnt)
    call faiEvolverelecaricheeridammele(out_Vnt, field_r_tessere, out_qt)
    call single_summation_tessere(out_qt, Vijn, states, tessere, out_Vijt)
    tmp=out_Vijt - subtract
    out_h_pcm = matmul(psi_prev, tmp)
end subroutine

