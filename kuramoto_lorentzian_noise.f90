! kuramoto_lorentzian_noise.f90
!
! Bryan Daniels
! 9-25-2004
! 1-04-2005 added noise
! 2-09-2005 finalized noise
!
! Numerical simulation of the Kuramoto model.
! N globally coupled oscillators.
! Gaussian noise term added each timestep.
!

implicit none

integer N
parameter (N=5000)  ! number of oscillators
double precision, dimension(N):: Theta, dTheta_dtau, omega, Theta_out, eta
double precision tau, Delta_tau, K, Delta_K, r, phi, gamma, pi, beta_squared, Delta_beta_squared
integer timesteps, i, j, t, num_K, l
K = 1.0d0          ! initial coupling
Delta_K = 0.2d0    ! K step size
Delta_tau = .01d0  ! time step size
timesteps = 5e5    ! number of time steps until we calculate r
num_K = 10         ! number of K values to test
gamma = 0.5d0      ! defines the width of the distribution g(omega)
beta_squared = 0.5d0  ! strength of noise
Delta_beta_squared = 0.5d0
pi = 4.d0*datan(1.d0)

open(unit=40, file='lorentzian_1000_gamma=0.5_beta_squared=varying.txt', status='unknown')
open(unit=50, file='Omega_1000_gamma=0.5_beta_squared=varying.txt', status='unknown')
open(unit=60, file='psi_1000_gamma=0.5_beta_squared=varying.txt', status='unknown')
open(unit=70, file='r_1000_gamma=0.5_beta_squared=varying.txt', status='unknown')

! set natural frequencies to lorentzian distribution
call lorentzian(gamma, N, omega)
close(40)

! loop over beta_squared values
do 950 l = 1,4
  write(70,*) "Beta^2 = ", beta_squared

! loop over K values
do 900 i = 1,num_K
  write(50,*) "K = ", K
  write(50,*) "Beta^2 = ", beta_squared
  write(60,*) "K = ", K
  write(60,*) "Beta^2 = ", beta_squared

  tau = 0.0d0

  ! initialize phases randomly
  call random_seed
  do 111 t=1,N
call random_number(Theta(t))
Theta(t) = Theta(t) * 2.d0*pi
111 continue

  do 800 t = 1,timesteps
    call random_array(eta,N,beta_squared/Delta_tau)
    call derivs(tau, Theta, dTheta_dtau, N, K, eta, omega)
    call rk4(Theta, dTheta_dtau, N, tau, Delta_tau, Theta_out, K, beta_squared, eta, omega)
    Theta = Theta_out

    tau = tau + Delta_tau

800       continue

  call find_order_param(Theta, r, phi, N)
  write(*,*) K, r
  write(70,*) K, r

  do 500 j = 1,N
    write(50,*) omega(j), "     ", dTheta_dtau(j)
    write(60,*) omega(j), "     ", MOD(Theta(j),2*pi)

500  continue
  K = K + Delta_K

900  continue
  K = 1.0d0
  beta_squared = beta_squared + Delta_beta_squared

950  continue

stop
end

! lorentzian
!
! produces an array of random values for the natural frequencies omega(i)
! uses lorentzian distribution (rejection method); produces values
!  from -10*gamma to 10*gamma
! takes a value for gamma (defines width), and N, the size of the array
! returns the array gamma with random values
subroutine lorentzian(gamma, N, omega)

implicit none
  integer N
double precision, dimension(N)::i_thermal, omega
double precision pi, gamma, random1, random2, p, p_max
integer d
call random_seed

pi = 4.d0*datan(1.d0)

p_max = 1/(pi*gamma)

do 75 d=1,N

50 call random_number(random1)
random1 = random1*20.d0*gamma - 10.d0*gamma
p = gamma / (pi * (gamma*gamma + random1*random1))
call random_number(random2)
random2 = random2*p_max
if (random2 > p) go to 50
omega(d) = random1

75 continue

return

end

! rk4
!
! rk4 uses the fourth-order runge-kutta method to advance the
!   solution over an interval h
! returns the advanced value yout
! uses subroutine derivs to obtain values for the derivatives
! (from Numerical Recipies)
! **with noise**
subroutine rk4(y,dydx,n,x,h,yout,K,beta_squared,eta,omega)

implicit none
  integer n
double precision, dimension(n) :: y, dydx, yout, yt, dyt, dym, omega, eta
double precision h, hh, h6, x, xh, K, beta_squared
integer i
external derivs

hh=h*0.5
h6=h/6.d0
xh=x+hh
do 11 i=1,n
yt(i)=y(i)+hh*dydx(i)
11  continue

  call derivs(xh,yt,dyt,n,K,eta,omega)
do 12 i=1,n
yt(i)=y(i)+hh*dyt(i)
12  continue

  ! use same noise array as last
call derivs(xh,yt,dym,n,K,eta,omega)
do 13 i=1,n
yt(i)=y(i)+h*dym(i)
dym(i)=dyt(i)+dym(i)
13  continue

  call derivs(x+h,yt,dyt,n,K,eta,omega)
do 14 i=1,n
yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
14  continue

return
end

! derivs
!
! calculates the time derivative of Theta using the Kuramoto model
! **with noise**
! (using equation 1-3 in notes plus noise)
! returns the array dTheta_dtau
subroutine derivs(tau, Theta, dTheta_dtau, N, K, eta, omega)
  implicit none
  integer N
  double precision, dimension(N)::Theta, dTheta_dtau, omega, eta
  double precision tau, K, r, phi
  integer i

  call find_order_param(Theta, r, phi, N)

  do 100 i=1,N

    dTheta_dtau(i) = omega(i) + K*r*dsin(phi-Theta(i)) + eta(i)

100  continue

  return
  end

! find_order_param
!
! computes the complex order parameter,
! returned as the variables r and phi,
! where r is the magnitude and phi is the angle
subroutine find_order_param(theta, r, phi, N)

implicit none
  integer N
double precision, dimension(N)::theta
double precision r, phi, real_sum, imag_sum
integer j

real_sum = 0.d0
imag_sum = 0.d0

do 200  j=1,N
real_sum = real_sum + dcos(theta(j))
imag_sum = imag_sum + dsin(theta(j))
200 continue
real_sum = real_sum/N
imag_sum = imag_sum/N

r = dsqrt((real_sum)**2 + (imag_sum)**2)
  phi = dacos(real_sum/r)

  return
end

! random_array(x)
!
! returns an array of N random numbers
! gaussian distribution, mean zero, width specified as parameter
subroutine random_array(array, N, width)
  implicit none
  double precision, dimension(N):: array
  double precision width
  integer N, i
  do 300 i=1,N
    call normal_random_num(array(i))
    array(i) = dsqrt(width) * array(i)

300  continue

  return
  end

! normal_random_num(x)
!
! returns a random value from a gaussian probability distribution
! with mean 0 and width 1.
! uses gaussian  distribution (rejection method); produces values
!  from -10 to 10
subroutine normal_random_num(x)

  implicit none
  double precision x, pi, p, p_max, random1, random2

  pi = 4.d0*datan(1.d0)
  p_max = 1/sqrt(2*pi)

  50 call random_number(random1)
random1 = random1*20.d0 - 10.d0
p = exp(-(random1*random1)/2) / sqrt(2*pi)
call random_number(random2)
random2 = random2*p_max
if (random2 > p) go to 50
x = random1
return
end
