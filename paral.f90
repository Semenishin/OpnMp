module Homework
    use omp_lib
    implicit none
    contains
        subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
        implicit none
        real(8),  dimension(:,:) :: A
        integer(4), intent(out) :: x1, y1, x2, y2
        integer(4) :: n, L, R, Up, Down, m, tmp,otr1,otr2,num_thread,i,amount_of_threads
        real(8), allocatable :: current_column(:),max_sum_arr(:), B(:,:)
        real(8) :: current_sum, max_sum,A_min,count1,count2
        logical :: transpos,otr
!         real(8), allocatable, dimension(:) :: max_sum   
        integer(4), allocatable, dimension(:):: X_1, X_2, Y_1, Y_2,Coords
        
        m = size(A, dim=1) 
        n = size(A, dim=2) 
        transpos = .FALSE.

        if (m < n) then 
            transpos = .TRUE.   
            B = transpose(A)
            m = size(B, dim=1) 
            n = size(B, dim=2) 
        else
            B = A     
        endif

        allocate(current_column(m))
        max_sum=-huge(0d0)
        x1=1
        y1=1
        x2=1
        y2=1

!$OMP PARALLEL
    !$OMP single
            amount_of_threads=omp_get_num_threads()
             allocate(max_sum_arr(amount_of_threads))
             allocate(X_1(amount_of_threads))
             allocate(X_2(amount_of_threads))
             allocate(Y_1(amount_of_threads))
             allocate(Y_2(amount_of_threads))
             max_sum_arr=max_sum
             X_1=1
             X_2=1
             Y_1=1
             Y_2=1
    !$OMP end single
    !$OMP   DO  SCHEDULE(dynamic) private(L,R,current_column,current_sum,Up,Down,x1,x2,y1,y2,num_thread,max_sum)
        do L=1, n
            current_column = B(:, L) 
            do R=L,n
                if (R > L) then 
                    current_column = current_column + B(:, R)
                endif
                call FindMaxInArray(current_column, current_sum, Up, Down) 
                if (current_sum > max_sum) then
                    max_sum = current_sum
                    num_thread=omp_get_thread_num()+1
                    max_sum_arr(num_thread)=max_sum
                    X_1(num_thread)=Up
                    X_2(num_thread)=Down
                    Y_1(num_thread)=L
                    Y_2(num_thread)=R
                endif
              end do
        end do     
      
    !$OMP END DO  
!$OMP END PARALLEL

        deallocate(current_column)
        max_sum=maxval(max_sum_arr)
        
        do i=1,amount_of_threads
            if(max_sum_arr(i)==max_sum) then
                 exit
            endif
        enddo
      
        x1=X_1(i)
        x2=X_2(i)
        y1=Y_1(i)
        y2=Y_2(i)

        deallocate(max_sum_arr)
        deallocate(X_1)
        deallocate(X_2)
        deallocate(Y_1)
        deallocate(Y_2)
     
        if (transpos) then  
            tmp = x1
            x1 = y1
            y1 = tmp
    
            tmp = y2
            y2 = x2
            x2 = tmp
            endif

        end subroutine


        subroutine FindMaxInArray(a, Sum, Up, Down)
            real(8), intent(in), dimension(:) :: a
            integer(4), intent(out) :: Up, Down
            real(8), intent(out) :: Sum
            real(8) :: cur_sum
            integer(4) :: minus_pos, i

            Sum = a(1)
            Up = 1
            Down = 1
            cur_sum = 0
            minus_pos = 0
            do i=1, size(a)
                cur_sum = cur_sum + a(i)
                 if (cur_sum > Sum) then
                     Sum = cur_sum
                     Up = minus_pos + 1
                     Down = i
                endif
         
            if (cur_sum < 0) then
                cur_sum = 0
                minus_pos = i
            endif

            enddo
     end subroutine FindMaxInArray
end module Homework







