program Main
    implicit none
    save

    integer :: seed, i, j, k
    integer, parameter :: states =2
    integer, parameter :: tessere = 3
    complex*16, dimension(tessere,states,states) :: qijn, Vijn 
    complex*16, dimension(tessere) :: wf_fwd, chi_prev
    real :: start, finish
    complex*16, dimension(states,states) :: subtract 

    complex, dimension(states) :: output
    seed = 1


    do i=1,tessere
       do j=1,states
          do k=1,states
              qijn(i,j,k) = cmplx ( i+j+k+0.1,i+j+k+0.1)
              Vijn(i,j,k) = cmplx ( i+j+k+0.1,i+j+k+0.1)
          end do
       enddo
    end do



    do i=1,states
        wf_fwd(i) = i
        write(*,*) wf_fwd(i)
        chi_prev = -i
    enddo



    do i=1,tessere
       do j=1,states
            subtract(i,j) = i+j
       enddo
    end do

    call doubleSum(chi_prev, wf_fwd, qijn, states, tessere, output) 

end program




