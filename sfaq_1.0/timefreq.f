cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM timefreq
      USE Const
      IMPLICIT NONE
      REAL(8),PARAMETER::rzero=0.0
      COMPLEX(8),PARAMETER::i=(0.0,1.0),czero=(0.0,0.0)
      INTEGER::it,Nt,iw,Nw,is,Ns,itmin,itmax
      REAL::wmin,wmax,dw,smin,smax,ds,height,s0,t1,t2,tramp1,tramp2
      COMPLEX::upper,lower
      REAL,ALLOCATABLE::t(:),y(:),s(:),w(:),win(:)
      COMPLEX,ALLOCATABLE::yw(:,:)
      NAMELIST/XPARA/Nt,Nw,wmin,wmax,smin,smax,Ns,s0,tramp1,tramp2
c numerics parameters
      OPEN(22,FILE='timefreq.cfg',STATUS='OLD')
      READ(22,NML=XPARA)
      CLOSE(22)
c input function y(t)
      ALLOCATE(t(Nt),y(Nt))
      DO it=1,Nt
        READ(*,*) t(it),y(it)
      END DO
c sI to atomic
      wmin=wmin/atu_eV
      wmax=wmax/atu_eV
      smin=smin/atu_s
      smax=smax/atu_s
      s0=s0/atu_s
      tramp1=tramp1/atu_s
      tramp2=tramp2/atu_s
      t=t/atu_s
c time window
      ALLOCATE(win(Nt))
      t1=t(1)+tramp1
      t2=t(Nt)-tramp2
      DO it=1,Nt
        IF (t(it)<t1) THEN
          win(it)=EXP(-((t1-t(it))/tramp1*2.0)**4)
        ELSE IF (t(it)>t2) THEN
          win(it)=EXP(-((t(it)-t2)/tramp2*2.0)**4)
        ELSE IF (t(it)>=t1.AND.t(it)<=t2) THEN
          win(it)=1.0
        END IF
      END DO
      y=y*win
c XFROG
      dw=(wmax-wmin)/REAL(Nw)
      Nw=Nw+1
      ALLOCATE(w(Nw))
      DO iw=1,Nw
        w(iw)=wmin+dw*REAL(iw-1)
      END DO
      ds=(smax-smin)/REAL(Ns)
      Ns=Ns+1
      ALLOCATE(s(Ns))
      DO is=1,Ns
        s(is)=smin+ds*REAL(is-1)
      END DO
      ALLOCATE(yw(Ns,Nw))
!$OMP PARALLEL PRIVATE(is,itmin,itmax,iw,it,lower,upper,height)
!$OMP DO SCHEDULE(Dynamic)
      DO is=1,Ns
        WRITE(*,5) 'delay(s)',s(is)*atu_s
        itmin=SUM(MINLOC(ABS(t-(s(is)-4.0*s0)))) ! the it for t(it) that is closest to s-4*s0
        itmax=SUM(MINLOC(ABS(t-(s(is)+4.0*s0)))) ! the it for t(it) that is closest to s+4*s0
        DO iw=1,Nw
          yw(is,iw)=czero
          DO it=itmin,itmax
            lower=EXP(-((t(it)-s(is))/s0)**2-i*t(it)*w(iw))*y(it)
            upper=EXP(-((t(it+1)-s(is))/s0)**2-i*t(it+1)*w(iw))*y(it+1)
            height=(t(it+1)-t(it))/2.0
            yw(is,iw)=yw(is,iw)+height*(lower+upper)
          END DO
          yw(is,iw)=yw(is,iw)*cFT
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      OPEN(35,FILE='xfr.dat',STATUS='UNKNOWN')
      DO is=1,Ns
        DO iw=1,Nw
          WRITE(35,10) s(is)*atu_s,w(iw)*atu_eV,
     &                 ABS(yw(is,iw))**2,yw(is,iw)
        END DO
        WRITE(35,*)
      END DO
      CLOSE(35)
      STOP
 5    FORMAT(A16,100ES16.7E3)
 10   FORMAT(100ES16.7E3)
      END PROGRAM timefreq
