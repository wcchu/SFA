      MODULE Const
      IMPLICIT NONE
      REAL,PARAMETER::
     & atu_m=5.2917721E-11,                                                ! 1 atu length in [m]
     & atu_s=2.4188843E-17,                                                ! 1 atu time in [s]
     & atu_fs=0.024188843,
     & atu_Vm=5.1422065E11,                                                ! 1 atu field strength in [V/m]
     & atu_Cm=8.4783533E-30,                                               ! 1 atu dipole in [Coul*m]
     & atu_J=4.3597442E-18,                                                ! 1 atu energy in [J]
     & atu_lam_nm=45.563354,                                                  ! 1 atu photon energy to wavelength in [nm]
     & cFT=0.3989422804,                                                ! 1/sqrt(2*pi)
c     & clight=2.9979246E8,                                              ! speed of light in [m/s]
     & clight=137.036,                                                  ! speed of light in atu
c     & eps0=8.8541878E-12,                                              ! free-space permittivity in [Coul/V/m]
     & eps0=0.07957747154,                                              ! free-space permittivity in atu
     & atu_Hz=6.579683920729E15,                                                  ! 1 atu photon energy in [1/s]
     & atu_PHz=6.579683920729,
     & atu_eV=27.21138505,
c     & Rgas=1.38065E-28,                                                ! R for ideal gas [m^3*bar/K]
     & Rgas=931.709,                                                    ! R for ideal gas [a0^3*bar/K]
     & pi=3.14159265359,
     & Icyc=3.51E16,                                                  ! 1 atu cycle-ave intensity in [W/cm2]
     & Ifun=6.44E15,                                                  ! fundamental units combined
     & reff=0.56553644,                                                ! reff = I_eff / I_center
     & rave=0.26951411,                                                ! rave = I_ave / I_center
     & u11=2.4048256,
     & a11=8.6584065
      END MODULE Const
