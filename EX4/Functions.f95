module Functions

    contains
        subroutine LoopMult(A,B,C,mm,nn)
            ! input variables
            integer, intent(in)                        ::  mm, nn
            integer                        ::  ii,jj,kk
            double precision,intent(in)                 ::  A(:,:),B(:,:)
            double precision,intent(out), dimension(:,:), allocatable    ::  C
            allocate(C(mm,nn))
            C=0.

            do ii=1,mm
                do jj=1,nn
                    do kk=1,size(A,dim=2)
                        C(ii,jj)=C(ii,jj)+A(ii,kk)*B(kk,jj)
                    end do      
                end do
            end do
        end subroutine

end module