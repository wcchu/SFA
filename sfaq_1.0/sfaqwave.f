cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MODULE Global
      USE Const
      USE Interp
      IMPLICIT NONE
      COMPLEX,PARAMETER::i=(0.0,1.0)
      INTEGER::
     & it,Nt,iv,Nv,kt,jt,itau1,itau2,Lt,ir,Nr,ofld,Ntin
      REAL::
     & Ip,E0,td,w0,tmin,tmax,dt,t0,CEP,eps,taumin,taumax,B,
     & dv,vmax,dr,rmax
      REAL,ALLOCATABLE::
     & t(:),E(:),A(:),v(:),r(:)
      COMPLEX,ALLOCATABLE::
     & exvr(:,:)

      CONTAINS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Input
      REAL::Ip_eV,I0_Wcm2,wlen_m,tmin_s,tmax_s,t0_s,tlen_s,
     & taumin_s,taumax_s,dt_s,vmax_ms,rmax_m
      NAMELIST/ATM/Ip_eV
      NAMELIST/FLD/ofld,I0_Wcm2,wlen_m,CEP,tlen_s,Ntin,t0_s
      NAMELIST/RTN/taumin_s,taumax_s,eps
      NAMELIST/DMN/dt_s,tmin_s,tmax_s,Lt,
     & Nv,vmax_ms,rmax_m,Nr
      OPEN(15,FILE='sfaqwave.cfg',STATUS='OLD')
      READ(15,NML=ATM)
      READ(15,NML=FLD)
      READ(15,NML=RTN)
      READ(15,NML=DMN)
      CLOSE(15)
      Ip=Ip_eV/atu_eV
      B=SQRT(2.0*Ip)
      E0=SQRT(I0_Wcm2/Icyc)
      td=tlen_s/atu_fs/1.665
      w0=(atu_lam_nm*1E-9)/wlen_m
      dt=dt_s/atu_s
      tmin=tmin_s/atu_s
      tmax=tmax_s/atu_s
      t0=t0_s/atu_s
      taumin=taumin_s/atu_s
      taumax=taumax_s/atu_s
      vmax=vmax_ms/(atu_m/atu_s)
      rmax=rmax_m/atu_m
      RETURN
      END SUBROUTINE Input

c Define the dimensions of calculation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Dimensions
c time
      Nt=INT((tmax-tmin)/dt)+1
      ALLOCATE(t(Nt),E(Nt),A(Nt))
      DO it=1,Nt
        t(it)=tmin+dt*REAL(it-1)
      END DO
c tau range
      itau1=NINT(taumin/dt)
      itau2=NINT(taumax/dt)
c velocity
      dv=2.0*vmax/REAL(Nv)
      Nv=Nv+1
      ALLOCATE(v(Nv))
      DO iv=1,Nv
        v(iv)=-vmax+dv*REAL(iv-1)
      END DO
c position
      dr=2.0*rmax/REAL(Nr)
      Nr=Nr+1
      ALLOCATE(r(Nr))
      DO ir=1,Nr
        r(ir)=-rmax+dr*REAL(ir-1)
      END DO
c exvr(v,r)
      ALLOCATE(exvr(Nv,Nr))
      FORALL (iv=1:Nv,ir=1:Nr)
        exvr(iv,ir)=cFT*EXP(i*v(iv)*r(ir))
      END FORALL
      RETURN
      END SUBROUTINE Dimensions

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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION dip(p)
      REAL::p
      COMPLEX::dip
c      dip=-i*SQRT(8.0*B**3/pi)*p/(p*p+B*B)**2                           ! 1D space, delta pot. with plane wave
      dip=-i*(2.0**3.5)*(B**2.5)/pi*p/(p*p+B*B)**3                      ! 3D space, H-like pot. with plane wave
      RETURN
      END FUNCTION dip

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gnd(pos)
      REAL::gnd,pos
c      gnd=SQRT(B)*EXP(-B*ABS(pos)) ! delta-func
      gnd=SQRT(B**3/pi)*EXP(-B*ABS(pos)) ! hydrogen
      RETURN
      END FUNCTION gnd

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION Act(ivel,itim,itau)
      INTEGER::ivel,itim,itau,it2
      REAL::inner,cc,Act
      cc=0.0
      DO it2=itim-itau,itim
        IF (it2<1) CYCLE
        inner=(v(ivel)-A(itim)+A(it2))**2/2.0+Ip
        IF (it2==itim-itau.OR.it2==itim) THEN
          cc=cc+0.5*dt*inner
        ELSE IF (it2>itim-itau.AND.it2<itim) THEN
          cc=cc+dt*inner
        END IF
      END DO
      Act=cc
      RETURN
      END FUNCTION Act

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Evolve
      INTEGER::OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,tid,nthreads
      INTEGER::it2,Nt2
      COMPLEX::aa
      COMPLEX,ALLOCATABLE::cv(:,:),c1(:),inner1(:),cr(:,:)
      Nt2=(Nt-1)/Lt+1
      ALLOCATE(cv(Nt2,Nv),c1(Nv),inner1(Nv),cr(Nt2,Nr))
      OPEN(90,FILE='thread.dat',STATUS='UNKNOWN')
!$OMP PARALLEL PRIVATE(it,it2,jt,kt,iv,ir,inner1,c1,aa,tid,nthreads)
!$OMP DO SCHEDULE(Dynamic)
c kt: launch, it: return, jt: drift duration
      DO it=1,Nt,Lt ! wide time mesh
        WRITE(*,10) t(it)*atu_s
c        tid=OMP_GET_THREAD_NUM()
c        WRITE(90,'(F12.5,I12)') t(it)*atu_fs,tid
c        IF (tid==0) THEN
c          nthreads=OMP_GET_NUM_THREADS()
c          write(90,*) 'number of threads',nthreads
c        END IF
        it2=(it-1)/Lt+1
        c1=(0.0,0.0)
        DO jt=itau1,itau2
          kt=it-jt
          IF (kt<1) CYCLE
          DO iv=1,Nv
            inner1(iv)=E(kt)*dip(v(iv)-A(it)+A(kt))*EXP(i*Act(iv,it,jt))
          END DO
          IF (jt==itau1.OR.jt==itau2) THEN
            c1=c1+dt*0.5*inner1
          ELSE IF (jt>itau1.AND.jt<itau2) THEN
            c1=c1+dt*inner1
          END IF
        END DO
        cv(it2,:)=i*EXP(i*Ip*t(it))*c1(:)
        DO ir=1,Nr
          aa=(0.0,0.0)
          DO iv=1,Nv
            aa=aa+dv*exvr(iv,ir)*cv(it2,iv)
          END DO
c          cr(it2,ir)=aa-EXP(i*Ip*t(it))*gnd(r(ir))                     ! free-electron + ground state
          cr(it2,ir)=aa                                                 ! only free-electron part
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      CLOSE(90)
      OPEN(40,FILE='cv.dat',STATUS='UNKNOWN')
      OPEN(50,FILE='cr.dat',STATUS='UNKNOWN')
      DO it=1,Nt,Lt
        it2=(it-1)/Lt+1
        DO iv=1,Nv
          WRITE(40,10) t(it)*atu_s,v(iv)*atu_m/atu_s,
     &                 ABS(cv(it2,iv))**2,cv(it2,iv)
        END DO
        WRITE(40,*)
        DO ir=1,Nr
          WRITE(50,10) t(it)*atu_s,r(ir)*atu_m,
     &                 ABS(cr(it2,ir))**2,cr(it2,ir)
        END DO
        WRITE(50,*)
      END DO
      CLOSE(40)
      CLOSE(50)
 5    FORMAT(A16,F16.7)
 10   FORMAT(10ES16.7E3)
 21   FORMAT(A16,10A16)
 22   FORMAT(2A16,10A16)
      RETURN
      END SUBROUTINE Evolve


      END MODULE Global

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      PROGRAM Main
      USE Global
      CALL Input
      CALL Dimensions
      CALL Field
      CALL Evolve
      STOP
      END PROGRAM Main


