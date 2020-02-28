subroutine matmult(M1, M2, d1, prod)
    integer :: d1
    complex, dimension(d1, d1),intent(in) :: M1
    complex, dimension(d1,d1),intent(in) :: M2
    complex, dimension(d1,d1), intent(out) :: prod
    prod = matmul(M1,M2)
end subroutine

subroutine dot(v, M2, d1, prod)
    integer :: d1
    complex, dimension(d1),intent(in) :: v
    complex, dimension(d1,d1),intent(in) :: M2
    complex, dimension(d1), intent(out) :: prod
    prod = matmul(v,M2)
end subroutine

subroutine dotv(v1, v2, d1, prod)
    integer :: d1
    real, dimension(d1),intent(in) :: v1
    real, dimension(d1),intent(in) :: v2
    real :: prod
    prod = dot_product(v1,v2)
end subroutine

