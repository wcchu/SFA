c SFA model for dipole oscillation in linear polarized field of
c arbitrary envelope
c Wei-Chun Chu / Max Planck Institute for the Science of Light 2015
c
c The theory/algorithm are based on:
c Lewenstein et al, Phys Rev A 49, 2117 (1994)
c
c The hydrogenlike dipole matrix elements are in convention defined by:
c Jin et al, Phys Rev A 79, 053413 (2009)
c
c Atomic parameters:
c   Ip = binding energy
c Field parameters:
c   I0 = peak pulse intensity
c   tlen = FWHM pulse duration
c   wlen = central wavelength
c   CEP = carrier envelope phase
c   t0 = peak time of pulse
c   alternatively the electric field can be numeric input
c Numeric approximation parameters:
c   taumin = lower limit of tau
c   taumax = upper limit of tau
c   eps = epsilon
c Dimensional parameters:
c   tmin = start of total physical time
c   tmax = end of total physical time
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MODULE Global
      USE Const
      USE Interp
      IMPLICIT NONE
      COMPLEX,PARAMETER::i=(0.0,1.0)
      INTEGER::
     & Nt,Lt,itau1,itau2,ofld,Ntin
      REAL::
     & Ip,E0,td,w0,tmin,tmax,dt,t0,CEP,eps,taumin,taumax,B
      REAL,ALLOCATABLE::
     & t(:),E(:),A(:),acc(:),x(:),Pst(:,:),Sst(:,:)
      COMPLEX,ALLOCATABLE::
     & cg(:),H1(:,:),H2(:,:)

      CONTAINS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Input
      REAL::Ip_eV,I0_Wcm2,wlen_m,tmin_s,tmax_s,t0_s,tlen_s,
     & taumin_s,taumax_s,dt_s
      NAMELIST/ATM/Ip_eV
      NAMELIST/FLD/ofld,I0_Wcm2,wlen_m,CEP,tlen_s,Ntin,t0_s
      NAMELIST/RTN/taumin_s,taumax_s,eps
      NAMELIST/DMN/dt_s,tmin_s,tmax_s,Lt
      OPEN(15,FILE='sfaqdip.cfg',STATUS='OLD')
      READ(15,NML=ATM)
      READ(15,NML=FLD)
      READ(15,NML=RTN)
      READ(15,NML=DMN)
      CLOSE(15)
c all units are converted from SI to atomic units
      Ip=Ip_eV/atu_eV
      B=SQRT(2.0*Ip)
      E0=SQRT(I0_Wcm2/Icyc)
      td=tlen_s/atu_s/1.665
      w0=(atu_lam_nm*1E-9)/wlen_m
      dt=dt_s/atu_s
      tmin=tmin_s/atu_s
      tmax=tmax_s/atu_s
      t0=t0_s/atu_s
      taumin=taumin_s/atu_s
      taumax=taumax_s/atu_s
      RETURN
      END SUBROUTINE Input

c Define the dimensions of calculation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Dimensions
      INTEGER::it
c time
      Nt=NINT((tmax-tmin)/dt)+1
      ALLOCATE(t(nt),E(nt),A(nt),x(nt),acc(nt),cg(Nt))
      DO it=1,nt
        t(it)=tmin+dt*REAL(it-1)
      END DO
c tau range
      itau1=NINT(taumin/dt)
      itau2=NINT(taumax/dt)
      ALLOCATE(Pst(Nt,itau1:itau2),Sst(Nt,itau1:itau2),
     & H1(Nt,itau1:itau2),H2(Nt,itau1:itau2))
      RETURN
      END SUBROUTINE Dimensions

c External electric field E(t) and vector potential A(t)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Field
      INTEGER::it
      REAL::tt,tmp,t1,E1
      REAL,ALLOCATABLE::tin(:),Ein(:),Ein2(:)
      SELECT CASE (ofld)
      CASE (1)
        E=0.0
        DO it=1,Nt
          tt=t(it)-t0
          E(it)=E0*EXP(-0.5*(tt/td)**2)*COS(w0*tt+CEP)
        END DO
      CASE (2)
        OPEN(17,FILE='field.inp',STATUS='OLD')
        ALLOCATE(tin(Ntin),Ein(Ntin),Ein2(Ntin))
        DO it=1,Ntin
          READ(17,*) tin(it),Ein(it)
        END DO
        CLOSE(17)
        tin=tin/atu_s
        Ein=Ein/atu_Vm
c interpolation from E_input to E_grid
c        write(*,*) 'Interpolating input field'
        CALL SplineR(tin,Ein,Ntin,0.0,0.0,Ein2)
        DO it=1,Nt
          t1=t(it)
          CALL SplIntR(tin,Ein,Ntin,Ein2,t1,E1)
          E(it)=E1
        END DO
      END SELECT
      A=0.0
      DO it=2,nt
        A(it)=A(it-1)-dt*E(it)
      END DO
c print
 11   FORMAT(ES13.5,10ES11.3)
 21   FORMAT(A13,10A11)
      OPEN(20,FILE='E.dat',STATUS='UNKNOWN')
      WRITE(20,21) '#t(s)','E(V/m)','A(V/m*s)'
      DO it=1,nt
        WRITE(20,11) t(it)*atu_s,E(it)*atu_Vm,A(it)*atu_Vm*atu_s
      END DO
      CLOSE(20)
      RETURN
      END SUBROUTINE Field

c Stationary momentum pst(t,tau)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      SUBROUTINE Stationary
      INTEGER::it,jt,it2
      REAL::a1,a2,a3,a4
      WRITE(*,*) 'preparing Pst & Sst'
      Pst=0.0
      Sst=0.0
!$OMP PARALLEL PRIVATE(it,jt,it2,a1,a2,a3,a4)
!$OMP DO SCHEDULE(Dynamic)
      DO it=1,Nt
        IF (MOD(it,Lt)==1) WRITE(*,'(ES16.7E3)') t(it)*atu_s
        DO jt=itau1,itau2
          IF (jt==0) THEN
            Pst(it,jt)=-A(it)
            Sst(it,jt)=0.0
          ELSE
            a1=-1.0/(dt*REAL(jt))
            a2=0.0
            DO it2=it-jt,it
              IF (it2==it-jt.OR.it2==it) THEN
                a2=a2+dt*0.5*A(it2)
              ELSE IF (it2>it-jt.AND.it2<it) THEN
                a2=a2+dt*A(it2)
              END IF
            END DO
            Pst(it,jt)=a1*a2
            a4=0.0
            DO it2=it-jt,it
              a3=0.5*(Pst(it,jt)+A(it2))**2+Ip
              IF (it2==it-jt.OR.it2==it) THEN
                a4=a4+dt*0.5*a3
              ELSE IF (it2>it-jt.AND.it2<it) THEN
                a4=a4+dt*a3
              END IF
            END DO
            Sst(it,jt)=a4
          END IF
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      WRITE(*,*) 'stationary momentum, action done'
      RETURN
      END SUBROUTINE Stationary

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ReScaFunc
      INTEGER::it,jt
      COMPLEX::dip1,dip2,H
      COMPLEX,ALLOCATABLE::fac(:)
      ALLOCATE(fac(itau1:itau2))
      DO jt=itau1,itau2 ! factor after integral over p
c        fac(jt)=(pi/(eps+i*0.5*dt*REAL(jt)))**0.5                       ! 1D space
        fac(jt)=(pi/(eps+i*0.5*dt*REAL(jt)))**1.5                       ! 3D space
      END DO
!$OMP PARALLEL PRIVATE(it,jt,dip1,dip2,H)
!$OMP DO SCHEDULE(Dynamic)
      DO it=1,Nt
        DO jt=itau1,itau2
          IF (it-jt<1) CYCLE
          dip1=dip(Pst(it,jt)+A(it))
          dip2=dip(Pst(it,jt)+A(it-jt))
          H=fac(jt)*dip2*E(it-jt)*EXP(-i*Sst(it,jt))
          H1(it,jt)=H*dip1 ! for depletion
          H2(it,jt)=H*CONJG(dip1) ! for dipole moment
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      WRITE(*,*) 'rescattering function done'
      RETURN
      END SUBROUTINE ReScaFunc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION dip(p)
      IMPLICIT NONE
      REAL::p
      COMPLEX::dip
c      dip=-i*SQRT(8.0*B**3/pi)*p/(p*p+B*B)**2                           ! 1D space, delta pot. with plane wave
      dip=-i*(2.0**3.5)*(B**2.5)/pi*p/(p*p+B*B)**3                      ! 3D space, H-like pot. with plane wave
      RETURN
      END FUNCTION dip

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Evolution
      INTEGER::it,jt
      COMPLEX::c2,c2a,c3
      COMPLEX,ALLOCATABLE::fac(:),c1(:)
      ALLOCATE(fac(itau1:itau2),c1(Nt))
      cg=(1.0,0.0)
!$OMP PARALLEL PRIVATE(it,jt)
!$OMP DO SCHEDULE(Dynamic)
      DO it=1,Nt
        IF (MOD(it,Lt)==1) WRITE(*,'(ES16.7E3)') t(it)*atu_s
        c1(it)=(0.0,0.0)
        c2=(0.0,0.0)
        c3=(0.0,0.0)
        DO jt=itau1,itau2
          IF (it-jt<1) CYCLE
          IF (jt==itau1.OR.jt==itau2) THEN
            c1(it)=c1(it)+dt*0.5*H1(it,jt)*cg(it-jt)
            c2=c2+dt*0.5*H1(it,jt)
            c3=c3+dt*0.5*H2(it,jt)*cg(it-jt)
          ELSE IF (jt>itau1.AND.jt<itau2) THEN
            c1(it)=c1(it)+dt*H1(it,jt)*cg(it-jt)
            c2=c2+dt*H1(it,jt)
            c3=c3+dt*H2(it,jt)*cg(it-jt)
          END IF
        END DO
c full calculation of cg
        IF (it>1) THEN
          cg(it)=cg(it-1)+dt*0.5*(E(it-1)*c1(it-1)+E(it)*c1(it))
        END IF
c approximation of cg
c        Gamma(it)=E(it)*c2
c        IF (it>1) THEN
c          c2a=c2a+dt*0.5*(Gamma(it-1)+Gamma(it))
c          cg(it)=EXP(c2a)
c        END IF
c dipole moment
        x(it)=2.0*REAL(i*CONJG(cg(it))*c3)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      acc=0.0
      DO it=2,nt-1
        acc(it)=(x(it+1)+x(it-1)-2.0*x(it))/dt/dt
      END DO
c print
      OPEN(60,FILE='ne.dat',STATUS='UNKNOWN')
      OPEN(62,FILE='cg.dat',STATUS='UNKNOWN')
      OPEN(70,FILE='x.dat',STATUS='UNKNOWN')
      OPEN(72,FILE='a.dat',STATUS='UNKNOWN')
      DO it=1,Nt
        WRITE(60,10) t(it)*atu_s,1.0-ABS(cg(it))**2
        WRITE(62,10) t(it)*atu_s,cg(it)
        WRITE(70,10) t(it)*atu_s,x(it)*atu_Cm
        WRITE(72,10) t(it)*atu_s,acc(it)*atu_Cm/atu_s/atu_s
      END DO
      CLOSE(60)
      CLOSE(62)
      CLOSE(70)
      CLOSE(72)
 10   FORMAT(10ES16.7E3)
      RETURN
      END SUBROUTINE Evolution


      END MODULE Global

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      PROGRAM Main
      USE Global
      CALL Input
      CALL Dimensions
      CALL Field
      CALL Stationary
      CALL ReScaFunc
      CALL Evolution
      STOP
      END PROGRAM Main


