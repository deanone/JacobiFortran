program jacobi
implicit none

	! Section for the declaration of variables used by the program
	! All declarations should be at the top of the program
	
	! Input
	real, dimension(6,6) :: A !coefficient matrix: a two dimensional array
	real, dimension(6) :: b !right-hand side vector: an one dimensional array
	integer :: max_iter = 50 !max number of iterations 
	real, dimension(6) :: x_init !initial vector
	real :: tol = 1.0e-5 !tolerance

	! Output
	real, dimension(6) :: x

	! Αuxiliary variables
	integer:: i, j, k
	real :: sum
	integer :: n
	real :: err
	
	n = size(b) ! or n = size(A, 1) or n = size(A, 2)

	print *, 'Size of my problem: ', n
	print *, 'Maximum number of iterations: ', max_iter
	print *, 'Tolerance: ', tol

	! Initialize x_init
	do i = 1, n
		x_init(i) = 0
	end do

	! Print x_init
	print *, 'x_init:'
	do i = 1, n
		write(*, fmt="(f3.0)", advance="no") x_init(i)
	end do
	write(*, '(/)')

	! Initialize x
	x = x_init

	! Print x
	print *, 'x:'
	do i = 1, n
		write(*, fmt="(f3.0)", advance="no") x(i)
	end do
	write(*, '(/)')

	!Assign values to A
	A = transpose(reshape((/ 4,-1,0,-1,0,0,-1,4,-1,0,-1,0,0,-1,4,0,0,-1,-1,0,0,4,-1,0,0,-1,0,-1,4,-1,0,0,-1,0,-1,4 /), shape(A)))

	!Print A
	print *, 'Coefficient matrix A:'
	do i = 1, n
		do j = 1, n
			write(*, fmt="(f3.0)", advance="no") A(i,j)
		end do
		write(*, '(/)')
	end do 

	!Assign values to b
	b = (/ 0,5,0,6,-2,6 /)

	!Print b
	print *, 'Right-hand side vector b:'
	do i = 1, n
		write(*, fmt="(f3.0)", advance="no") b(i)
	end do
	write(*, '(/)')

	!Main Jacobi loop
	print *, 'error during Jacobi iterations'
	do k = 1, max_iter
		do i = 1, n
			sum = 0
			do j = 1, (i - 1)
				sum = sum + A(i, j) * x_init(j)
			end do

			do j = (i + 1), n
				sum = sum + A(i, j) * x_init(j)
			end do 

			x(i) = (b(i) - sum) / A(i, i)
		end do

		!do i = 1, n
		!	write(*, fmt="(f20.14)", advance="no") x(i)
		!end do
		!write(*, '(/)')

		err = NORM2(x - x_init)
		write(*, fmt="(f20.14)", advance="no") err
		write(*, '(/)')

		if (err < tol) then
			exit
		end if

		x_init = x
	end do

	print *, "Iterations: ", k
	print *, 'x:'
	do i = 1, n
		write(*, fmt="(f20.14)", advance="no") x(i)
	end do
	write(*, '(/)')

end program jacobi