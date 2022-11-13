!Module containing functions and parameters used for calculations
! in Hormones_v2.f90

module functions

  use parameter_values, only : RP

  implicit none

  !Values used for calculations in main program
  real(kind=RP), parameter, public  :: pi = atan(1._RP) * 4 !Pi



  contains


  !Function for autocorelation
  function autocorrelatedE(lastE, alpha, StochScale, ValueOfE_min, ValueOfE_max, E_max) result(phiVE)
    !Things going in
    real(kind=RP), intent(in) :: lastE
    real(kind=RP), intent(in) :: alpha      !Autocorrelation coefficient - process has mean 0 and sd 1
    real(kind=RP), intent(in) :: StochScale !SD of autocorrelated process at which extreme of E-categories-axis is positioned.
    real(kind=RP), intent(in) :: ValueOfE_min, ValueOfE_max
    integer, intent(in)       :: E_max

    !Things going out
    real(kind=RP), dimension(1:3) :: phiVE
    
    !Variables within the function
    real(kind=RP) :: lastphi, phi
    real(kind=RP) :: ValueOfE
    real(kind=RP) :: E_real
    
    ValueOfE = -1000._RP
    
    !This is only true if StochScale != 0 and if ValueOfE_max != ValueOfE_min
    lastphi = (real(E_max-1,RP)*(ValueOfE_min-1._RP)+(lastE-1._RP)*(ValueOfE_max-ValueOfE_min))/(StochScale*real(E_max-1,RP)) 
      
    !Stochscale is variance in category units, number of sd that correspond to category=1 (and max)
    do while ((ValueOfE < ValueOfE_min) .OR. (ValueOfE > ValueOfE_max))
      !See Ripa and Lundberg 1996 for scaling of variance.
      phi = alpha*lastphi + sqrt(1._RP-alpha**2._RP)*rnorm()  !Phi-process is autocorrelated with mean 0 and sd 1
      ValueOfE = (StochScale * phi) + 1._RP            !Change the standard deviation and variance 
      
      !Saving phi as lastphi in case the loop has to run again
      lastphi = phi
      
      !If AutoCorr/alpha==1 exit loop to avoid infinite loop
      if(alpha==1._RP) exit
    end do
      
    E_real = 1._RP + ((ValueOfE-ValueOfE_min)/(ValueOfE_max-ValueOfE_min))*real(E_max-1,RP)
    
    !Saving values in array for returning
    phiVE(1) = phi
    phiVE(2) = ValueOfE
    phiVE(3) = E_real

  end function autocorrelatedE



  !Function for cumulative calculations of a 1D array containing real numbers
  ! and returning these values as an array with one more position than the raw array
  ! where the first position 0.
  function cumucalc(RawArray) result(CumuArray)
    real(kind=RP), dimension(:), intent(in)      :: RawArray !Array coming in
    real(kind=RP), dimension(1:size(RawArray)+1) :: CumuArray  !Array going out

    !Variables in function
    integer :: num !Integer used for loop

    !Adding 0 to the first place in the array
    CumuArray(1) = 0._RP

    !For every position in RawArray calculate the cumulative values
    do num = 1, size(RawArray)
        CumuArray(num+1) = CumuArray(num) + RawArray(num)
    end do

  end function cumucalc



  !Function for picking random values in a 1D array containing cumulative probabilities
  ! and returning the position of the nearest value in the array by ROUNDING DOWN
  function rcumu(CumuArray) result(ArrayPosition)
    real(kind=RP), dimension(:), intent(in) :: CumuArray !Array coming in
    integer                                 :: ArrayPosition !Position going out

    !Variables in function
    real(kind=RP) randomnr !Random number
    integer i !Counter in loop

    !Calculating a random number between the min and max values in the array
    randomnr = rreal(CumuArray(1),CumuArray(size(CumuArray)))

    !Goes through every number in the array and saves the position of the
    ! nearest value in the array by rounding down
    do i = 1, size(CumuArray)-1
      if (CumuArray(i)<=randomnr .AND. randomnr <= CumuArray(i+1)) ArrayPosition = i
    end do

  end function rcumu




  !FUNCTIONS FOR GENERATING RANDOM NUMBERS
  
  
  !Subroutine for starting the random number generators from a set seed
  ! useful for debugging
  subroutine fixed_seed(fix_seed)
    integer :: n 
    integer, intent(in) :: fix_seed
    integer, dimension(:), allocatable :: seed
    
    call random_seed(size = n)
    allocate(seed(n))
    seed = fix_seed
    call random_seed(PUT = seed)
    deallocate(seed)
    
  end subroutine fixed_seed



  !Functions that generates a normally distributed real random number
  !Slightly modified from RNORM_VAL_R8 in BASE_RANDOM.f90 in HEDTOOLS
  function rnorm() result(fn_val)
    !*******************************************************************************
    ! RNORM
    ! PURPOSE: Returns a normally distributed pseudo-random number with zero mean
    !          and unit variance. This version uses double precision, as a high
    !          precition from random numbers is not needed.
    ! NOTES:   Adapted from the following Fortran 77 code
    !          ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
    !          THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
    !          VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
    ! WARNING: This is a simple conversion of single precision 32 bit function
    !          statistical properties of the generated distribution are to be
    !          checked yet.
    ! The function rnorm() returns a normally distributed pseudo-random
    ! number with zero mean and unit variance.
    ! The algorithm uses the ratio of uniforms method of A.J. Kinderman
    ! and J.F. Monahan augmented with quadratic bounding curves.
    ! CODE SOURCE: http://www.netlib.org/ (random.f90)
    !*******************************************************************************

    REAL(8) :: fn_val !Random number

    !Local variables used for calculations
    REAL(8) :: s=0.449871_8, t=-0.386595_8,              &
               a=0.19600_8, b=0.25472_8,                 &
               r1=0.27597_8, r2=0.27846_8, u, v, x, y, q
    REAL(8) :: half = 0.5_8 !Moved here from header of the original module

    !Generate P = (u,v) uniform in rectangle enclosing acceptance region

    DO
      CALL random_number(u)
      CALL random_number(v)
      v = 1.7156_8 * (v - half)

    !     Evaluate the quadratic form
      x = u - s
      y = ABS(v) - t
      q = x**2 + y*(a*y - b*x)

    !Accept P if inside inner ellipse
      IF (q < r1) EXIT
    !Reject P if outside outer ellipse
      IF (q > r2) CYCLE
    !Reject P if outside acceptance region
      IF (v**2 < -4.0_8*LOG(u)*u**2) EXIT
    END DO

    !Return ratio of P's coordinates as the normal deviate
    fn_val = v/u
    RETURN

  end function rnorm



  !Functions that generates a random real number between from a and b
  ! (a <= r < b)
  function rreal(a,b) result(randreal)
    real(kind=RP), intent(in) :: a, b
    real(kind=RP)             :: randreal

    !Calling a random number between 0 and 1
    call random_number(randreal)

    !Scaling the random number to the user defined interval
    randreal = (randreal*(b - a)) + a

  end function rreal



  !Functions that generates a random integer number between a and b
  !Modified from RAND_I in BASE_RANDOM.f90 in HEDTOOLS
  function rint(a, b) result (randint)
    ! Trivial random integer (a <= r <= b)

    implicit none
    integer             :: randint
    integer, intent(in) :: a, b
    real(kind=RP)       :: randreal

    !Getting a random real number between 0 and 1
    call random_number(randreal)

    ! Use int or floor? floor returns more homogeneous numbers if negative
    ! integers are in the range. With positive arguments int and floor eqivalent.
    ! Note: int(-0.2)= 0, floor(-0.2)=-1
    !       int(+0.2)= 0, floor(+0.2)= 0 , so with int there is a small local
    !       raise in the density function, not good.
    randint = a + floor(randreal * (b - a + 1))

  end function rint



end module functions
