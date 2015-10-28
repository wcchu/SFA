cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM spectr
      USE Const
      IMPLICIT NONE
      COMPLEX,PARAMETER::i=(0.0,1.0)
      INTEGER::it,Nt,iw,Nw,icol,Ncol
      REAL::tramp1,tramp2,t1,t2,wmin,wmax,dw,height,dum
      COMPLEX::lower,upper
      REAL,ALLOCATABLE::t(:),f(:,:),w(:),win(:)
      COMPLEX,ALLOCATABLE::fw(:,:)
      NAMELIST/PAR/Nt,tramp1,tramp2,wmin,wmax,Nw,Ncol
c configuration
      OPEN(15,FILE='spectr.cfg',STATUS='OLD')
      READ(15,PAR)
      CLOSE(15)
c input f(t)
      ALLOCATE(t(Nt),f(Ncol,Nt))
      DO it=1,Nt
        READ(*,*) t(it),(f(icol,it),icol=1,Ncol)
      END DO
c sI to atomic
      tramp1=tramp1/atu_s
      tramp2=tramp2/atu_s
      wmin=wmin/atu_eV
      wmax=wmax/atu_eV
      t=t/atu_s
c      f=f/atu_Cm*atu_s*atu_s
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
      FORALL (icol=1:Ncol,it=1:Nt)
        f(icol,it)=f(icol,it)*win(it)
      END FORALL
      OPEN(21,FILE='win.dat',STATUS='UNKNOWN')
      DO it=1,Nt
        WRITE(21,10) t(it)*atu_fs,win(it),(f(icol,it),icol=1,Ncol)
      END DO
      CLOSE(21)
c fourier transform
      dw=(wmax-wmin)/REAL(Nw)
      Nw=Nw+1
      ALLOCATE(w(Nw),fw(Ncol,Nw))
!$OMP PARALLEL PRIVATE(iw,it,lower,upper,height)
!$OMP DO SCHEDULE(Dynamic)
      DO iw=1,Nw
        w(iw)=wmin+dw*REAL(iw-1)
        fw(:,iw)=(0.0,0.0)
        DO it=1,Nt-1
          height=(t(it+1)-t(it))/2.0
          DO icol=1,Ncol
            lower=EXP(-i*t(it)*w(iw))*f(icol,it)
            upper=EXP(-i*t(it+1)*w(iw))*f(icol,it+1)
            fw(icol,iw)=fw(icol,iw)+height*(lower+upper)
          END DO
        END DO
        fw(:,iw)=fw(:,iw)*cFT
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      OPEN(40,FILE='amp.dat',STATUS='UNKNOWN')
      OPEN(45,FILE='int.dat',STATUS='UNKNOWN')
      DO iw=1,Nw
        WRITE(40,10) w(iw)*atu_eV,(fw(icol,iw),icol=1,Ncol)
        WRITE(45,10) w(iw)*atu_eV,(ABS(fw(icol,iw))**2,icol=1,Ncol)
      END DO
      CLOSE(40)
      CLOSE(45)
      STOP
 10   FORMAT(100ES16.7E3)
      END PROGRAM spectr
