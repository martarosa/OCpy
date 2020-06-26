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
    complex*16, dimension(states,states) :: h_mdm_0


        h_mdm_0(1,1) = -3.1673748114875971E-004
        h_mdm_0(1,2) =   3.0074473421421086E-011
        h_mdm_0(1,3) =   -1.0057815130298987E-006
        h_mdm_0(1,4) =   2.5508172787542155E-006
        h_mdm_0(1,5) =   -2.8537359261022263E-007
        h_mdm_0(1,6) =   4.1166749930074150E-009
        h_mdm_0(1,7) =   4.0117437758327982E-005
        h_mdm_0(1,8) =   -4.3202735421925415E-006
        h_mdm_0(1,9) =   -1.6701464801351182E-006
        h_mdm_0(1,10) =   -8.2003814269903722E-005
        h_mdm_0(1,11) =  1.5317996883426760E-010
        h_mdm_0(2,1) =   3.0074473421421086E-011
        h_mdm_0(2,2) =   -1.3311771205228905E-004
        h_mdm_0(2,3) =   6.2790977160431575E-006
        h_mdm_0(2,4) =   7.5953328830808310E-009
        h_mdm_0(2,5) =   -1.2255569350941540E-006
        h_mdm_0(2,6) =   -4.2508041689286444E-002
        h_mdm_0(2,7) =   -1.0984901080683236E-007
        h_mdm_0(2,8) =   2.9296178671772056E-006
        h_mdm_0(2,9) =   6.2018299331079573E-012
        h_mdm_0(2,10)=   1.5075409334921121E-006
        h_mdm_0(2,11) =  -4.2835511815177703E-002
        h_mdm_0(3,1) =   -1.0057815130298987E-006
        h_mdm_0(3,2) =   6.2790977160431575E-006
        h_mdm_0(3,3) =   4.2749102209352172E-002
        h_mdm_0(3,4) =   1.8022228660542358E-006
        h_mdm_0(3,5) =   3.0684699279072555E-006
        h_mdm_0(3,6) =   6.1695842688281749E-009
        h_mdm_0(3,7) =   -3.4778395949477732E-005
        h_mdm_0(3,8) =   -4.5288651566854887E-006
        h_mdm_0(3,9) =   3.4364297984744141E-008
        h_mdm_0(3,10)=   3.8965649510141062E-006
        h_mdm_0(3,11) =  6.6728223218692782E-009
        h_mdm_0(4,1) =    2.5508172787542155E-006
        h_mdm_0(4,2) =    7.5953328830808310E-009
        h_mdm_0(4,3) =    1.8022228660542358E-006
        h_mdm_0(4,4) =    4.2746636476031809E-002
        h_mdm_0(4,5) =   -6.4806105189146222E-011
        h_mdm_0(4,6) =    2.9285599714200690E-006
        h_mdm_0(4,7) =    4.4257581021519440E-006
        h_mdm_0(4,8) =    3.5134993978056966E-005
        h_mdm_0(4,9) =    3.8969873724595548E-006
        h_mdm_0(4,10)=   -2.7880320713678917E-008
        h_mdm_0(4,11) =  -3.8327014026415203E-006
        h_mdm_0(5,1) =   -2.8537359261022263E-007
        h_mdm_0(5,2) =   -1.2255569350941540E-006
        h_mdm_0(5,3) =    3.0684699279072555E-006
        h_mdm_0(5,4) =   -6.4806105189146222E-011
        h_mdm_0(5,5) =    4.2920801560623133E-002
        h_mdm_0(5,6) =   -7.8075702034202114E-006
        h_mdm_0(5,7) =    1.1204077574058755E-005
        h_mdm_0(5,8) =   -1.0646618136421990E-007
        h_mdm_0(5,9) =   -1.5057360861352863E-006
        h_mdm_0(5,10)=    1.0562817681824461E-006
        h_mdm_0(5,11) =  -1.6127599685994748E-007
        h_mdm_0(6,1) =    4.1166749930074150E-009
        h_mdm_0(6,2) =   -4.2508041689286444E-002
        h_mdm_0(6,3) =    6.1695842688281749E-009
        h_mdm_0(6,4) =    2.9285599714200690E-006
        h_mdm_0(6,5) =   -7.8075702034202114E-006
        h_mdm_0(6,6) =    4.2748457260188534E-002
        h_mdm_0(6,7) =   -4.2506977672953070E-002
        h_mdm_0(6,8) =   -1.1302941532438900E-005
        h_mdm_0(6,9) =    1.0345040835858801E-006
        h_mdm_0(6,10)=   -4.2603861457066454E-002
        h_mdm_0(6,11) =  -4.2832552442104321E-002
        h_mdm_0(7,1) =    4.0117437758327982E-005
        h_mdm_0(7,2) =   -1.0984901080683236E-007
        h_mdm_0(7,3) =   -3.4778395949477732E-005
        h_mdm_0(7,4) =    4.4257581021519440E-006
        h_mdm_0(7,5) =    1.1204077574058755E-005
        h_mdm_0(7,6) =   -4.2506977672953070E-002
        h_mdm_0(7,7) =    4.2754451374714814E-002
        h_mdm_0(7,8) =    4.3899519361477078E-006
        h_mdm_0(7,9) =    3.6087318358728139E-009
        h_mdm_0(7,10)=   -1.1905343756191178E-005
        h_mdm_0(7,11) =  -1.5174954011805351E-006
        h_mdm_0(8,1) =   -4.3202735421925415E-006
        h_mdm_0(8,2) =    2.9296178671772056E-006
        h_mdm_0(8,3) =   -4.5288651566854887E-006
        h_mdm_0(8,4) =    3.5134993978056966E-005
        h_mdm_0(8,5) =   -1.0646618136421990E-007
        h_mdm_0(8,6) =   -1.1302941532438900E-005
        h_mdm_0(8,7) =    4.3899519361477078E-006
        h_mdm_0(8,8) =   -9.7820330654336304E-005
        h_mdm_0(8,9) =   -4.2844456949788112E-002
        h_mdm_0(8,10)=   -4.6272003324547460E-006
        h_mdm_0(8,11) =   3.8769343070200773E-007
        h_mdm_0(9,1) =   -1.6701464801351182E-006
        h_mdm_0(9,2) =    6.2018299331079573E-012
        h_mdm_0(9,3) =    3.4364297984744141E-008
        h_mdm_0(9,4) =    3.8969873724595548E-006
        h_mdm_0(9,5) =   -1.5057360861352863E-006
        h_mdm_0(9,6) =    1.0345040835858801E-006
        h_mdm_0(9,7) =    3.6087318358728139E-009
        h_mdm_0(9,8) =   -4.2844456949788112E-002
        h_mdm_0(9,9) =    4.2746632718318849E-002
        h_mdm_0(9,10)=   -4.6605521042411026E-006
        h_mdm_0(9,11) =  -4.2634225441137374E-002
        h_mdm_0(10,1) =   -8.2003814269903722E-005
        h_mdm_0(10,2) =    1.5075409334921121E-006
        h_mdm_0(10,3) =    3.8965649510141062E-006
        h_mdm_0(10,4) =   -2.7880320713678917E-008
        h_mdm_0(10,5) =    1.0562817681824461E-006
        h_mdm_0(10,6) =   -4.2603861457066454E-002
        h_mdm_0(10,7) =   -1.1905343756191178E-005
        h_mdm_0(10,8) =   -4.6272003324547460E-006
        h_mdm_0(10,9) =   -4.6605521042411026E-006
        h_mdm_0(10,10) =    4.2746636629112253E-002
        h_mdm_0(10,11) =  -8.7185162045226782E-011
        h_mdm_0(11,1) =    1.5317996883426760E-010
        h_mdm_0(11,2) =   -4.2835511815177703E-002
        h_mdm_0(11,3) =    6.6728223218692782E-009
        h_mdm_0(11,4) =   -3.8327014026415203E-006
        h_mdm_0(11,5) =   -1.6127599685994748E-007
        h_mdm_0(11,6) =   -4.2832552442104321E-002
        h_mdm_0(11,7) =   -1.5174954011805351E-006
        h_mdm_0(11,8) =    3.8769343070200773E-007
        h_mdm_0(11,9) =   -4.2634225441137374E-002
        h_mdm_0(11,10)=   -8.7185162045226782E-011
        h_mdm_0(11,11) =   5.3498186566630895E-006



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

