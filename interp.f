c Interpolate complex function y(x)
c 20140919 Remove the use of MINLOC statement, use bisection to fine the
c          location in the x array closest to x0. Speed is faster by
c          >6 times.
      MODULE Interp
      CONTAINS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Natural spline for 2nd derivative of complex function y(x)
      SUBROUTINE SplineC(x,y,n,yp1,ypn,y2)
      INTEGER::j,n
      REAL::sig
      COMPLEX::pxy,ryx,yp1,ypn,qn,un
c      REAL,ALLOCATABLE::x(:)
c      COMPLEX,ALLOCATABLE::y(:),y2(:),u(:)
c      ALLOCATE(u(n))
      REAL::x(n)
      COMPLEX::y(n),y2(n),u(n)
      y2(1)=(-0.5,0.0)
      u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      DO j=2,n-1
        sig=(x(j)-x(j-1))/(x(j+1)-x(j-1))
        pxy=sig*y2(j-1)+2.0
        y2(j)=(sig-1.0)/pxy
        ryx=(y(j+1)-y(j))/(x(j+1)-x(j))
     &     -(y(j)-y(j-1))/(x(j)-x(j-1))
        u(j)=(6.0*ryx/(x(j+1)-x(j-1))-sig*u(j-1))/pxy
      END DO
      qn=(0.5,0.0)
      un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
      DO j=n-1,1,-1
        y2(j)=y2(j)*y2(j+1)+u(j)
      END DO
      RETURN
      END SUBROUTINE SplineC

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Interpolate complex function y(x) for y(x0)
      SUBROUTINE SplIntC(x,y,n,y2,x0,y0)
      INTEGER::k0,klo,khi,n
      REAL::x0,h,A,B
      COMPLEX::y0
c      REAL,ALLOCATABLE::x(:)
c      COMPLEX,ALLOCATABLE::y(:),y2(:)
      REAL::x(n)
      COMPLEX::y(n),y2(n)
c      k0=SUM(MINLOC(ABS(x-x0)))
c      IF (x0>x(k0).AND.k0<n) THEN
c        klo=k0
c        khi=k0+1
c      ELSE IF (x0<x(k0).AND.k0>1) THEN
c        khi=k0
c        klo=k0-1
c      ELSE ! x0 is on top of x(k0)
c        y0=y(k0)
c        RETURN
c      END IF
      klo=1
      khi=n
      DO WHILE (khi-klo>1)
        k0=(khi+klo)/2
        IF (x(k0)>x0) THEN
          khi=k0
        ELSE
          klo=k0
        END IF
      END DO
      h=x(khi)-x(klo)
      A=(x(khi)-x0)/h
      B=(x0-x(klo))/h
      y0=A*y(klo)+B*y(khi)
     &  +((A**3-A)*y2(klo)+(B**3-B)*y2(khi))*h*h/6.0
      RETURN
      END SUBROUTINE SplIntC

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Natural spline for 2nd derivative of real function y(x)
      SUBROUTINE SplineR(x,y,n,yp1,ypn,y2)
      INTEGER::j,n
      REAL::sig,pxy,ryx,yp1,ypn,qn,un
c      REAL,ALLOCATABLE::x(:),y(:),y2(:),u(:)
c      ALLOCATE(u(n))
      REAL::x(n),y(n),y2(n),u(n)
      y2(1)=(-0.5,0.0)
      u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      DO j=2,n-1
        sig=(x(j)-x(j-1))/(x(j+1)-x(j-1))
        pxy=sig*y2(j-1)+2.0
        y2(j)=(sig-1.0)/pxy
        ryx=(y(j+1)-y(j))/(x(j+1)-x(j))
     &     -(y(j)-y(j-1))/(x(j)-x(j-1))
        u(j)=(6.0*ryx/(x(j+1)-x(j-1))-sig*u(j-1))/pxy
      END DO
      qn=(0.5,0.0)
      un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
      DO j=n-1,1,-1
        y2(j)=y2(j)*y2(j+1)+u(j)
      END DO
      RETURN
      END SUBROUTINE SplineR

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Interpolate real function y(x) for y0(x0)
      SUBROUTINE SplIntR(x,y,n,y2,x0,y0)
      INTEGER::k0,klo,khi,n
      REAL::x0,h,A,B,y0
c      REAL,ALLOCATABLE::x(:),y(:),y2(:)
      REAL::x(n),y(n),y2(n)
c      k0=SUM(MINLOC(ABS(x-x0)))
c      IF (x0>x(k0).AND.k0<n) THEN
c        klo=k0
c        khi=k0+1
c      ELSE IF (x0<x(k0).AND.k0>1) THEN
c        khi=k0
c        klo=k0-1
c      ELSE ! x0 is on top of x(k0)
c        y0=y(k0)
c        RETURN
c      END IF
      klo=1
      khi=n
      DO WHILE (khi-klo>1)
        k0=(khi+klo)/2
        IF (x(k0)>x0) THEN
          khi=k0
        ELSE
          klo=k0
        END IF
      END DO
      h=x(khi)-x(klo)
      A=(x(khi)-x0)/h
      B=(x0-x(klo))/h
      y0=A*y(klo)+B*y(khi)
     &  +((A**3-A)*y2(klo)+(B**3-B)*y2(khi))*h*h/6.0
      RETURN
      END SUBROUTINE SplIntR
      

      END MODULE Interp
