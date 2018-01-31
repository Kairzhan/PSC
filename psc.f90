!
!   psc - FVM method for solving scalar convection in 2D uniform grid
!
!   Created by Kairzhan Karzhaubayev on 3/17/12.
!   Copyright 2012 Kairzhan Karzhaubayev. All rights reserved.
!
program main
implicit none
real*8 :: xmin,xmax,ymin,ymax
real*8, allocatable, dimension(:) :: x,xc,y,yc
real*8, allocatable, dimension(:,:) :: me,mn,ae,aw,an,as,ap,q
real*8, allocatable, dimension(:,:) :: new_phi, phi
real*8 :: dx,dy
integer :: n,m
integer :: i,j,l
real*8 :: dens, diff
real*8 :: aed,awd,and,asd,aec,awc,anc,asc
real*8 :: tmp,res,eps
integer :: scheme
logical :: debug

eps=0.001
dens=1
diff=0.1
scheme=2 ! 1 - CDS ; 2 - UDS
debug=.false.

xmin=0
xmax=1
n=80
ymin=0
ymax=1
m=80

! ---------------------------------------------------------------------------------------------------
!
! Setting up a grid in X direction
!
!print *,"enter xmin,xmax,n"
!read (*,*) xmin,xmax,n

allocate(x(0:n+2),xc(n+2))

dx=(xmax-xmin)/n

x(0)=xmin
x(1)=xmin
do i=2,n
    x(i)=x(i-1)+dx
enddo
x(n+1)=xmax
x(n+2)=xmax

do i=1,n+2
    xc(i)=0.5*(x(i)+x(i-1))
enddo

! ---------------------------------------------------------------------------------------------------
!
! Setting up a grid in Y direction
!
!print *,"enter ymin,ymax,m"
!read (*,*) ymin,ymax,m

allocate(y(0:m+2),yc(m+2))

dy=(ymax-ymin)/m

y(0)=ymin
y(1)=ymin
do j=2,m
    y(j)=y(j-1)+dy
enddo
y(m+1)=ymax
y(m+2)=ymax

do j=1,m+2
    yc(j)=0.5*(y(j)+y(j-1))
enddo

! ---------------------------------------------------------------------------------------------------
!
! Calculating mass flux
!
allocate(me(n+1,m),mn(n,m+1))
do i=1,n+1
    do j=1,m
        me(i,j)=dens*(x(i))*(dy)
    enddo
enddo
do i=1,n
    do j=1,m+1
        mn(i,j)=-dens*(y(j))*(dx)
    enddo
enddo

! ---------------------------------------------------------------------------------------------------
!
! Coefficients
!
allocate(ae(0:n+1,0:m+1),aw(0:n+1,0:m+1),an(0:n+1,0:m+1),as(0:n+1,0:m+1),ap(n,m),q(n,m))
ae=0
aw=0
an=0
as=0
ap=0
q=0
do i=1,n
    do j=1,m
        ! Central difference for diffusion
        aed=-diff*dy/(xc(i+2)-xc(i+1))
        awd=-diff*dy/(xc(i+1)-xc(i))
        and=-diff*dx/(yc(j+2)-yc(j+1))
        asd=-diff*dx/(yc(j+1)-yc(j))

        ! CDS or UDS ?
        if (scheme==1) then ! 
            aec=me(i+1,j)/2.
            awc=-me(i,j)/2.
            anc=mn(i,j+1)/2.
            asc=-mn(i,j)/2.
        elseif (scheme==2) then ! UDS scheme
            aec=min(me(i+1,j),0.)
            awc=-max(me(i,j),0.)
            anc=min(mn(i,j+1),0.)
            asc=-max(mn(i,j),0.)
        else
            print *,"unknown conv. scheme"
            stop(1)
        endif

        ae(i,j)=aed+aec
        aw(i,j)=awc+awd
        an(i,j)=anc+and
        as(i,j)=asc+asd
        ap(i,j)=-(ae(i,j)+aw(i,j)+an(i,j)+as(i,j))
        q(i,j)=0
    enddo
enddo

! ---------------------------------------------------------------------------------------------------
!
! Boundary conditions
!
do j=1,m
    ! West
    q(1,j)=q(1,j)-aw(1,j)*( &
            (1.-((yc(j)+yc(j+1))/2.-ymin)/(ymax-ymin)) &
    )
    aw(1,j)=0

    ! East
    ap(n,j)=ap(n,j)+ae(n,j)
    ae(n,j)=0
enddo

do i=1,n
    ! South
    ap(i,1)=ap(i,1)+as(i,1)
    as(i,1)=0

    ! North
    q(i,m)=q(i,m)-an(i,m)*0.0
    an(i,m)=0
enddo

if (debug) then
    write(*,'(6(A7))') 'AE','AW','AN','AS','AP','Q'
    do i=1,n
    do j=1,m
    write(*,'(6(F7.4,A1))') ae(i,j),' ',aw(i,j),' ',an(i,j),' ',as(i,j),' ',ap(i,j),' ',q(i,j),' '
    enddo
    enddo
endif

! ---------------------------------------------------------------------------------------------------
!
! Simple Gauss-Zeidel iteration
!
allocate(phi(0:n+1,0:m+1),new_phi(0:n+1,0:m+1))
phi=0
new_phi=0
do l=1,2000
    do i=1,n
        do j=1,m
            tmp=ae(i,j)*phi(i+1,j) &
                +aw(i,j)*new_phi(i-1,j) &
                +an(i,j)*phi(i,j+1) &
                +as(i,j)*new_phi(i,j-1)
            new_phi(i,j)=(q(i,j)-tmp)/ap(i,j)
        enddo
    enddo
    res=0
    do i=1,n
        do j=1,m
            res=res+abs(phi(i,j)-new_phi(i,j))
        enddo
    enddo
    if (res<eps) goto 5
    do i=1,n
        do j=1,m
            phi(i,j)=new_phi(i,j)
        enddo
    enddo
enddo
4 print *,"***** not converged in",l,"iterations"
goto 6
5 print *,"converged in",l,"iterations"
6 continue

! ---------------------------------------------------------------------------------------------------
!
! Output to file
!
open(unit=111,file="output.dat",status="replace")
write(111,*) 'TITLE = "Stagnation conv"'
write(111,*) 'Variables="X","Y","Phi", "U", "V", "UV","m_e","m_n"'
write(111,*) 'Zone I=',n, 'J=', m, 'F=POINT'
do i=1,n
    do j=1,m
        write(111,*) i,j,phi(i,j),(x(i)+x(i+1))/2.,-(y(j)+y(j+1))/2.,-(x(i)+x(i+1))/2.*(y(j)+y(j+1))/2., (me(i,j)),(mn(i,j))
    enddo
enddo
close(111)

open(unit=112,file="output.txt",status="replace")
do i=1,n
    write(112,10) (phi(i,j),j=1,m)
enddo
10 FORMAT(100F6.3)
close(112)

end program main
