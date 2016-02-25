c     
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program dynamic - run molecular or stochastic dynamics  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     'dynamic' computes a molecular or stochastic dynamics 
c     trajectory in one of the standard statistical mechanical 
c     ensembles using any of several possible integration methods
c
c
      program dynamic_amoeba_waterBox_THz
      use sizes
      use atoms
      use atomid
      use bath
c      use bond
      use boxes
      use bound
      use charge
      use dipole
      use files
      use group
      use inform
      use iounit
      use keys
      use mdstuf
      use moldyn
      use moment
      use mpole
      use polar
      use potent
      use rgddyn
      use solute
      use stodyn
      use units
      use usage
      implicit none

      integer i,istep,istep0,nstep
      integer mode,next
      real*8 dt,dtdump
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string

      INTEGER NMAX,TMAX,TMAX2,LMAX
      PARAMETER (NMAX=5000,TMAX=5001,TMAX2=5001)

      INTEGER NMAX2
      PARAMETER (NMAX2=5000)
 
      INTEGER ISTOREV,DTCALCV,ISKIPV,ISTORED,DTCALCD,ISKIPD,NUM
      integer J,K,L,M,Np,Nc,Nw,II,JJ,II2
      integer ISTEP2,TT,DTC,ISTOREM,DTCALCM,ISKIPM
      integer wlist(NMAX2)
      real*8 xcm,ycm,zcm,norm,dx,dy,dz,Rg,T
      real*8 XC(NMAX),YC(NMAX),ZC(NMAX)


      INTEGER COUNTM(TMAX),COUNTV(TMAX),COUNTD(TMAX2)
      INTEGER countmass
      INTEGER RESTARTTIMESTEP,RESTARTFLAG,SAVE_LENGTH                ! restart variables
      INTEGER SAVE_DATA_FREQ,SAVE_REST_FREQ,t_stp                    ! restart variables
      INTEGER O_IDX,H_IDX,LYR_IDX,ATM_IDX                            ! restart variables
      real*8 XO(NMAX,TMAX),VXO(NMAX,TMAX),VXH(NMAX,TMAX)
      real*8 YO(NMAX,TMAX),VYO(NMAX,TMAX),VYH(NMAX,TMAX)
      real*8 ZO(NMAX,TMAX),VZO(NMAX,TMAX),VZH(NMAX,TMAX)
      real*8 MXC,MYC,MZC,deltax,deltay,deltaz

ccc      real*8 MXI(NMAX,TMAX),MYI(NMAX,TMAX),MZI(NMAX,TMAX)
      real*8 MXI(TMAX),MYI(TMAX),MZI(TMAX)

ccc      real*8 MXIC(NMAX),MYIC(NMAX),MZIC(NMAX)

      real*8 MXCT,MYCT,MZCT

      real*8 DR2,MSDO(TMAX)
      real*8 VACFO(TMAX),VACFH(TMAX)
      real*8 VX(NMAX),VY(NMAX),VZ(NMAX),eksum,ekin,temp
      real*8 drsq
      real*8 DCF(TMAX2)

      real*8 tempMSDO,tempVACFO,tempVACFH                       ! temporary local variables for omp
      real*8 tempCOUNTM,tempCOUNTV                              ! temporary local variables for omp
      real*8 privateVACFO,privateVACFH                          ! temporary local variables for omp

      real*8 summass,drsqmin

      integer KK,NGRP2

      real*8 startTime, finishTime, netStartTime, netFinishtime         ! to get timings
      real*8 wall,cpu       ! To get timings

c     Wall time
      INTEGER nb_ticks_initial,nb_ticks_final,nb_ticks_max,nb_ticks_sec
      INTEGER nb_ticks
      REAL elapsed_time  ! real time in seconds

      CHARACTER FILE1*99,FILE2*99,FILE3*99,FILE4*99,FILE5*99,FILE6*99
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                       c
c     These files save the values from the tape, for use to restart the simulation      c
c                                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CHARACTER FILE7*99,FILE8*99,FILE9*99  ! These files save the X, Y and Z coordinates for the Oxygens of the waters.
      CHARACTER FILE10*99,FILE11*99,FILE12*99 ! These files save the X, Y and Z velocities for the Oxygens of the waters.
      CHARACTER FILE13*99,FILE14*99,FILE15*99 ! These files save the X, Y and Z velocities for the Hydrogens of the waters.
      CHARACTER FILE16*99,FILE17*99,FILE18*99 ! These files save the X, Y and Z dipole moments for the 13 layers.
      CHARACTER FILE19*99,FILE24*99           ! This file writes out the file which contains the index upto which the simulation has progressed, in order to restart from that point.
      CHARACTER FILE20*99,FILE21*99,FILE22*99 ! These files write out the COUNTM, COUNTV and COUNTD values so that they can be multiplied into the values of MSD, VACF_O, VACF_H and DCF read in from the restart files.
      CHARACTER FILE23*99 ! This file writes out the summass array.
ccc     call cpu_time(startTime)
ccc     call cpu_time(netStartTime)
ccc     CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
ccc     CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
ccc      call settime
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BEEMAN'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
      end do
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,160)
  160          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if
c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,330)
  330    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,340)
  340    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,350)
  350    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,360)
  360    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,370)
  370    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,380)
  380    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
c
c     set modified group variables for amoeba
c
      NGRP2=N/3
c
c     read time intervals for calculations
c     
      call nextarg (string,exist)
      if (exist)  read (string,*)  ISTOREM !length of MSD    
      call nextarg (string,exist)
      if (exist)  read (string,*)  DTCALCM !time-step for MSD
      call nextarg (string,exist)
      if (exist)  read (string,*)  ISKIPM !interval to calculate MSD
      call nextarg (string,exist)
      if (exist)  read (string,*)  ISTOREV !length of VACF
      call nextarg (string,exist)
      if (exist)  read (string,*)  DTCALCV !time-step for VACF
      call nextarg (string,exist)
      if (exist)  read (string,*)  ISKIPV !interval to calculate VACF
      call nextarg (string,exist)
      if (exist)  read (string,*)  ISTORED !length of DCF
      call nextarg (string,exist)
      if (exist)  read (string,*)  DTCALCD !time-step for DCF
      call nextarg (string,exist)
      if (exist)  read (string,*)  ISKIPD !interval to calculate DCF
      call nextarg (string,exist)
      if (exist)  read (string,*)  SAVE_DATA_FREQ ! Freq to write the dcf and other production arrays to file
      call nextarg (string,exist)
      if (exist)  read (string,*)  SAVE_REST_FREQ  ! Freq to write out coordinates and velocities for restart files.
      ISKIPM=DTCALCM*INT(ISKIPM/DTCALCM)! make multiple of DTCALCM
      ISKIPV=DTCALCV*INT(ISKIPV/DTCALCV)! make multiple of DTCALCV
      ISKIPD=DTCALCD*INT(ISKIPD/DTCALCD)! make multiple of DTCALCD
c
c     read file names for writing output
c         
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE1 !name of dipole moment file
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE2 !name of MSD file
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE3 !name of O vacf file
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE4 !name of H vacf file
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE5 !name of dipole correlation file
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE6 !name of Rg file
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE7 !store file for X positions of Oxygens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE8 !store file for Y positions of Oxygens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE9 !store file for Z positions of Oxygens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE10 !store file for X velocities of Oxygens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE11 !store file for Y velocities of Oxygens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE12 !store file for Z velocities of Oxygens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE13 !store file for X velocities of Hydrogens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE14 !store file for Y velocities of Hydrogens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE15 !store file for Z velocities of Hydrogens of waters
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE16 !store file for X dipole moments for layers
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE17 !store file for Y dipole moments for layers
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE18 !store file for Z dipole moments for layers
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE19 !restart file which stores the index of the time step to which the simulation has progressed
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE20 !File with the count values for the msd array.
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE21 !File with the count values for the vacf array.
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE22 !File with the count values for the dcf array.
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE23 !File with the summass vector.
      call nextarg (string,exist)
      if (exist)  read (string,*)  FILE24 !Test file to check if the right thing is being read in.
c
c     open output files
c     
c     INQUIRE(FILE=FILE1,EXIST=EXIST)
c        IF (EXIST) THEN
c           OPEN (UNIT=201,FILE=FILE1,STATUS='OLD',
c    &           ACTION='WRITE',POSITION='APPEND')
c        ELSE
c           OPEN (UNIT=201,FILE=FILE1,STATUS='UNKNOWN')
c        ENDIF
c     INQUIRE(FILE=FILE6,EXIST=EXIST)
c        IF (EXIST) THEN
c           OPEN (UNIT=206,FILE=FILE6,STATUS='OLD',
c    &           ACTION='WRITE',POSITION='APPEND')
c        ELSE
c           OPEN (UNIT=206,FILE=FILE6,STATUS='UNKNOWN')
c        ENDIF
c
c     initialize time-correlation functions
c
      DO DTC=1,ISTOREM+1
         MSDO(DTC)=0d0
         COUNTM(DTC)=0d0
      ENDDO
      DO DTC=1,ISTOREV+1
         VACFO(DTC)=0d0
         VACFH(DTC)=0d0
         COUNTV(DTC)=0d0
      ENDDO
      DO DTC=1,ISTORED+1 
         DCF(DTC)=0d0
         COUNTD(DTC)=0d0
      ENDDO
c
c     count water atoms 
c
      Nw=0
      do I=1,N
         Nw=Nw+1
         wlist(Nw)=I
      enddo
c
c     X/Y/Z = boxed water coord, XC/YC/ZC = unboxed water coords
c
      DO I=1,Nw,3
         II=wlist(I) ! Oxygen index
         DO J=0,2
            JJ=wlist(I)+J ! all water indices
            XC(JJ)=X(JJ)-DNINT((X(JJ)-X(II))/XBOX)*XBOX
            YC(JJ)=Y(JJ)-DNINT((Y(JJ)-Y(II))/YBOX)*YBOX
            ZC(JJ)=Z(JJ)-DNINT((Z(JJ)-Z(II))/ZBOX)*ZBOX
         ENDDO
      ENDDO
c
c     calc dipole moment
c
      MXC=0.0d0
      MYC=0.0d0
      MZC=0.0d0
      DO II=1,npole
         K=ipole(II)
         MXC=MXC+X(K)*rpole(1,ii)+rpole(2,ii)+uind(1,ii)
         MYC=MYC+Y(K)*rpole(1,ii)+rpole(3,ii)+uind(2,ii)
         MZC=MZC+Z(K)*rpole(1,ii)+rpole(4,ii)+uind(3,ii)
      ENDDO
c     
c     save initial coords to tape
c
      TT=ISTOREM
      DO I=1,Nw,3
         K=wlist(I) ! Oxygen index
         XO(I,TT)=XC(K)
         YO(I,TT)=YC(K)
         ZO(I,TT)=ZC(K)
c      add velocities to tape         
         II=(I+2)/3
         J=Nw/3+II ! I (J) - Hydrogen 1 (2)   
         L=K+1 ! Hydrogen 1
         M=K+2 ! Hydrogen 2
         VXO(II,TT)=V(1,K)
         VYO(II,TT)=V(2,K)
         VZO(II,TT)=V(3,K)
         VXH(II,TT)=V(1,L)
         VYH(II,TT)=V(2,L)
         VZH(II,TT)=V(3,L)
         VXH(J,TT)=V(1,M)
         VYH(J,TT)=V(2,M)
         VZH(J,TT)=V(3,M)
      ENDDO
c      add dipole moments to tape
ccc      DO I=1,N
ccc         MXI(I,TT)=MXIC(I)
ccc         MYI(I,TT)=MYIC(I)
ccc         MZI(I,TT)=MZIC(I)
ccc      ENDDO
      MXI(TT)=MXC
      MYI(TT)=MYC
      MZI(TT)=MZC
      summass=0d0
      countmass=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     If simulation is being restarted, the ISTEP value is read         c
c     in from the restarttimestep.txt file to start from that point     c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         RESTARTFLAG=0
         RESTARTTIMESTEP=0
c        INQUIRE(FILE=FILE19,EXIST=EXIST)
c        IF (EXIST) THEN
c           OPEN (UNIT=219,FILE=FILE19, STATUS='OLD', ACTION='READ')
c              READ(219,*) RESTARTTIMESTEP
c              RESTARTFLAG=1
c           CLOSE(219)
c        ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                 c
c     All the arrays that have been written out are read back into memory         c
c     These include: MSD, H_vacf, O_vacf, M, M2, dcf, dcf2,                       c
c     X0, Y0, Z0, VXO, VYO, VZO, VXH, VYH, VZH, MXI, MYI, MZI, MXI2, MYI2, MZI2   c
c                                                                                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This reads in the MSD, O_VACF, H_vacf and dcf files and 
c       updates the MSD, O_vacf, H_vacf and dcf arrays, which are zero
c       since the simulation has been restarted.
c
c       INQUIRE(FILE=FILE2,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=202,FILE=FILE2, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=203,FILE=FILE3, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=204,FILE=FILE4, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=205,FILE=FILE5, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=220,FILE=FILE20, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=221,FILE=FILE21, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=222,FILE=FILE22, STATUS='OLD', ACTION='READ')
c             DO I=1,ISTOREV
c                READ(202,*) t_stp,MSDO(I)
c                READ(203,*) t_stp,VACFO(I)
c                READ(204,*) t_stp,VACFH(I)
c                READ(205,*) t_stp,DCF(I)
c                READ(220,*) t_stp,COUNTM(I)
c                READ(221,*) t_stp,COUNTV(I)
c                READ(222,*) t_stp,COUNTD(I)
c             ENDDO
c          CLOSE(202)
c       ENDIF
c
c       This file reads in the summass array.
c
c       INQUIRE(FILE=FILE23,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=223,FILE=FILE23, STATUS='OLD', ACTION='READ')
c             READ(223,*) summass
c          CLOSE(223)
c       ENDIF
c
c       This multiplies the values of the MSD, VACF_O, VACF_H and DCF read in
c       from the files by the appropriate COUNTM, COUNTV and COUNTD values to
c       get the actual values for those quantities
c
c       INQUIRE(FILE=FILE19,EXIST=EXIST)
c       IF (EXIST) THEN
c          DO I = 1, ISTOREV
c             MSDO(I)=MSDO(I)*COUNTM(I)
c             VACFO(I)=(VACFO(I)*COUNTV(I)*
c    &         (VACFO(1)/COUNTV(1)))
c             VACFH(I)=(VACFH(I)*COUNTV(I)*
c    &         (VACFH(1)/COUNTV(1)))
c             DCF(I)=DCF(I)*COUNTD(I)*summass
c          ENDDO
c       ENDIF

c
c       This reads in the Oxygen X, Y, Z coordinates and velocities file,
c       written out from the previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE7,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=207,FILE=FILE7, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=208,FILE=FILE8, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=209,FILE=FILE9, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=210,FILE=FILE10, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=211,FILE=FILE11, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=212,FILE=FILE12, STATUS='OLD', ACTION='READ')
c             IF(RESTARTTIMESTEP.LT.ISTOREV) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREV
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(207,*) (XO(O_IDX,DTC),O_IDX=1,Nw)
c                READ(208,*) (YO(O_IDX,DTC),O_IDX=1,Nw)
c                READ(209,*) (ZO(O_IDX,DTC),O_IDX=1,Nw)
c                READ(210,*) (VXO(O_IDX,DTC),O_IDX=1,Nw/3)
c                READ(211,*) (VYO(O_IDX,DTC),O_IDX=1,Nw/3)
c                READ(212,*) (VZO(O_IDX,DTC),O_IDX=1,Nw/3)
c             ENDDO
c          CLOSE(207)
c          CLOSE(208)
c          CLOSE(209)
c          CLOSE(210)
c          CLOSE(211)
c          CLOSE(212)
c       ENDIF
c
c       This reads in the Hydrogen X, Y, Z velocities file, written out from the 
c       previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE13,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=213,FILE=FILE13, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=214,FILE=FILE14, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=215,FILE=FILE15, STATUS='OLD', ACTION='READ')
c             IF (RESTARTTIMESTEP.LT.ISTOREV) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREV
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(213,*) (VXH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c                READ(214,*) (VYH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c                READ(215,*) (VZH(H_IDX,DTC),H_IDX=1,2*Nw/3)
c             ENDDO
c          CLOSE(213)
c          CLOSE(214)
c          CLOSE(215)
c       ENDIF
c
c       This reads in the X, Y, Z dipole moments file, written out from the 
c       previous run as the simulation has been restarted.
c
c       INQUIRE(FILE=FILE16,EXIST=EXIST)
c       IF (EXIST) THEN
c          OPEN (UNIT=216,FILE=FILE16, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=217,FILE=FILE17, STATUS='OLD', ACTION='READ')
c          OPEN (UNIT=218,FILE=FILE18, STATUS='OLD', ACTION='READ')
c             IF (RESTARTTIMESTEP.LT.ISTOREM) THEN
c                SAVE_LENGTH = RESTARTTIMESTEP
c             ELSE 
c                SAVE_LENGTH = ISTOREM
c             END IF
c             DO DTC=1,SAVE_LENGTH
c                READ(216,*) (MXI(DTC))
c                READ(217,*) (MYI(DTC))
c                READ(218,*) (MZI(DTC))
c             ENDDO
c          CLOSE(216)
c          CLOSE(217)
c          CLOSE(218)
c       ENDIF
ccc       call cpu_time(finishTime)
ccc       write(70,*) 'Time for initialization =',finishTime-startTime
ccc       CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc       nb_ticks = nb_ticks_final - nb_ticks_initial
ccc       IF (nb_ticks_final < nb_ticks_initial) 
ccc    &           nb_ticks = nb_ticks + nb_ticks_max
ccc       elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc       write(80,*) 'Wall time for initialization =',elapsed_time
c
c
c     This is where the simulation actually starts
c     integrate equations of motion to take a time step
c
      do istep0 = 1, nstep
ccc        call cpu_time(startTime)
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
         if (integrate .eq. 'VERLET') then
            call verlet (istep0,dt)
c            call verlet_4site_alexCode(istep0,dt)
         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep (istep0,dt)
         else if (integrate .eq. 'BUSSI') then
            call bussi (istep0,dt)
         else if (integrate .eq. 'NOSE-HOOVER') then
            call nose (istep0,dt)
         else if (integrate .eq. 'GHMC') then
            call ghmcstep (istep0,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep0,dt)
         else if (integrate .eq. 'RESPA') then
            call respa (istep0,dt)
         else
c            call beeman_4site(istep0,dt)
             call beeman(istep0,dt)
         end if
ccc        call cpu_time(finishTime)
ccc        write(71,*) 'Time for integration step',istep0,
ccc    &               '=',finishTime-startTime
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc        nb_ticks = nb_ticks_final - nb_ticks_initial
ccc        IF (nb_ticks_final < nb_ticks_initial) 
ccc    &            nb_ticks = nb_ticks + nb_ticks_max
ccc        elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc        write(81,*) 'Wall time for integration step',istep0,
ccc    &               ' =',elapsed_time
c     
c      Take coordinates out of the box for diffusion calculations
c
         DO I=1,N
            XC(I)=X(I)-DNINT((X(I)-XC(I))/XBOX)*XBOX
            YC(I)=Y(I)-DNINT((Y(I)-YC(I))/YBOX)*YBOX
            ZC(I)=Z(I)-DNINT((Z(I)-ZC(I))/ZBOX)*ZBOX
         ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     ISTEP is recaliberated if the simulation is being restarted       c
c     ISTEP=current_simulation_ISTEP_value+ISTEP value read in          c
c     from where previous simulation died.                              c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ISTEP=ISTEP0+RESTARTTIMESTEP
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c     The new modified value of ISTEP is written out to use for         c
c     subsequent restarts.                                              c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        IF(MOD(ISTEP,SAVE_REST_FREQ).EQ.0) THEN
c           OPEN (UNIT=219,FILE=FILE19, STATUS='REPLACE')
c              WRITE(219,'(I8)') ISTEP
c           CLOSE(219)
c        ENDIF
c
c        calculate dipole moment
c
ccc        call cpu_time(startTime)
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
         IF(MOD(ISTEP,DTCALCD).EQ.0) THEN
            MXC=0.0d0
            MYC=0.0d0
            MZC=0.0d0
            DO II=1,npole
                K=ipole(II)
                MXC=MXC+X(K)*rpole(1,ii)+rpole(2,ii)+uind(1,ii)
                MYC=MYC+Y(K)*rpole(1,ii)+rpole(3,ii)+uind(2,ii)
                MZC=MZC+Z(K)*rpole(1,ii)+rpole(4,ii)+uind(3,ii)
            ENDDO
c           WRITE(201,'(4ES16.8)') DT*DBLE(ISTEP),MXC,MYC,MZC
         ENDIF
ccc        call cpu_time(finishTime)
ccc        write(74,*) 'Time for dipole moment at step',istep0,
ccc    &               '=',finishTime-startTime
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc        nb_ticks = nb_ticks_final - nb_ticks_initial
ccc        IF (nb_ticks_final < nb_ticks_initial) 
ccc    &               nb_ticks = nb_ticks + nb_ticks_max
ccc        elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc        write(84,*) 'Wall time for dipole moment at step',istep0,
ccc    &               '=',elapsed_time
c          
c        calc msd & vacf and dcf
c        
         DO I=1,Nw
            K=wlist(I)    
            VX(K)=V(1,K)
            VY(K)=V(2,K)
            VZ(K)=V(3,K)
         ENDDO
         IF (MOD(ISTEP,ISKIPM).EQ.0.AND.ISTEP.GE.ISTOREM*DTCALCM) THEN
ccc           call cpu_time(startTime)
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
            DO ISTEP2=ISTEP-ISTOREM*DTCALCM,ISTEP,DTCALCM
               TT=MOD(ISTEP2/DTCALCM-1,ISTOREM)+1        
               IF (TT.EQ.0) TT=ISTOREM
               DTC=(ISTEP-ISTEP2)/DTCALCM
               tempMSDO=MSDO(DTC+1)
               tempVACFO=VACFO(DTC+1)
               tempVACFH=VACFH(DTC+1)
               tempCOUNTM=COUNTM(DTC+1)
               tempCOUNTV=COUNTV(DTC+1)
!$OMP PARALLEL default(private) shared(wlist,
!$OMP& TT,NW,XO,YO,ZO,XC,YC,ZC,VXO,VYO,VZO,VXH,VYH,VZH,VX,VY,VZ,ISTEP,
!$OMP& ISTEP2)
!$OMP& shared(tempMSDO,tempVACFO,tempVACFH,tempCOUNTM,tempCOUNTV)
!$OMP DO reduction(+:tempMSDO,tempVACFO,tempVACFH,tempCOUNTM,tempCOUNTV)
               DO I=1,Nw,3
                  K=wlist(I) ! Oxygen
                  IF (ISTEP2.LT.ISTEP) THEN
                     DX=XO(I,TT)-XC(K)
                     DY=YO(I,TT)-YC(K)
                     DZ=ZO(I,TT)-ZC(K)
c                    write(*,*) I,TT,XO(I,TT)
                  ELSE
                     DX=0d0
                     DY=0d0
                     DZ=0d0
                  ENDIF
                  DR2=DX*DX+DY*DY+DZ*DZ
                  tempMSDO=tempMSDO+DR2
                  tempCOUNTM=tempCOUNTM+1
c
c             Calculate vacf
c
                  II=(I+2)/3
                  J=Nw/3+II ! I (J) - Hydrogen 1 (2)           
                  L=K+1 ! Hydrogen 1
                  M=K+2 ! Hydrogen 2
                  IF (ISTEP2.LT.ISTEP) THEN
                     privateVACFO=
     +               VXO(II,TT)*VX(K)+VYO(II,TT)*VY(K)+VZO(II,TT)*VZ(K)
                     privateVACFH=
     +               VXH(II,TT)*VX(L)+VYH(II,TT)*VY(L)+VZH(II,TT)*VZ(L)
     +               +VXH(J,TT)*VX(M)+VYH(J,TT)*VY(M)+VZH(J,TT)*VZ(M)
                  ELSE
                     privateVACFO=
     +               VX(K)*VX(K)+VY(K)*VY(K)+VZ(K)*VZ(K)
                     privateVACFH=
     +               VX(L)*VX(L)+VY(L)*VY(L)+VZ(L)*VZ(L)
     +               +VX(M)*VX(M)+VY(M)*VY(M)+VZ(M)*VZ(M)
                  ENDIF
                  tempVACFO=tempVACFO+privateVACFO
                  tempVACFH=tempVACFH+privateVACFH
                  tempCOUNTV=tempCOUNTV+1
               ENDDO
!$OMP END DO
!$OMP END PARALLEL
               MSDO(DTC+1)=tempMSDO
               VACFO(DTC+1)=tempVACFO
               VACFH(DTC+1)=tempVACFH
               COUNTM(DTC+1)=tempCOUNTM
               COUNTV(DTC+1)=tempCOUNTV
c
c             Calculate DCF
c
ccc               DO II=1,N
ccc                  DCF(DTC+1)=DCF(DTC+1)+
ccc     &               MXI(II,TT)*MXC+MYI(II,TT)*MYC+MZI(II,TT)*MZC
ccc               ENDDO
               IF (ISTEP2.LT.ISTEP) THEN
                  DCF(DTC+1)=DCF(DTC+1)+MXI(TT)*MXC+MYI(TT)*MYC+
     &                       MZI(TT)*MZC
               ELSE
                  DCF(DTC+1)=DCF(DTC+1)+MXC*MXC+MYC*MYC+
     &                       MZC*MZC
               ENDIF
               COUNTD(DTC+1)=COUNTD(DTC+1)+1
               summass = 0d0
               do i=1,N
                  summass=summass+mass(i)
               enddo    
               countmass=countmass+1
            ENDDO
ccc           call cpu_time(finishTime)
ccc           write(75,*) 'Time for calc msd, vacf, dcf at step',
ccc    &                  istep0,'=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(85,*) 'Wall time for calc msd, vacf, dcf',istep0,
ccc    &                  '=',elapsed_time
            IF (MOD(ISTEP,SAVE_DATA_FREQ).EQ.0) THEN
ccc           call cpu_time(startTime)
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c           OPEN(UNIT=202,FILE=FILE2,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=203,FILE=FILE3,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=204,FILE=FILE4,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=205,FILE=FILE5,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=220,FILE=FILE20,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=221,FILE=FILE21,STATUS='REPLACE',ACTION='WRITE')
c           OPEN(UNIT=222,FILE=FILE22,STATUS='REPLACE',ACTION='WRITE')
c           DO DTC=0,ISTOREM
c              WRITE(202,'(2ES16.8)') DT*DBLE(DTC*DTCALCM),
c    &           MSDO(DTC+1)/DBLE(COUNTM(DTC+1))
c              WRITE(203,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             ((VACFO(DTC+1)/(DBLE(COUNTV(DTC+1))))
c    &             /(VACFO(1)/DBLE(COUNTV(1))))
c              WRITE(204,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             ((VACFH(DTC+1)/(DBLE(COUNTV(DTC+1))))
c    &             /(VACFH(1)/DBLE(COUNTV(1))))
               WRITE(205,'(2ES16.8)') DT*DBLE(DTC*DTCALCD),
     &             DCF(DTC+1)/DBLE(COUNTD(DTC+1))/summass
c              WRITE(220,'(2ES16.8)') DT*DBLE(DTC*DTCALCM),
c    &             DBLE(COUNTM(DTC+1))
c              WRITE(221,'(2ES16.8)') DT*DBLE(DTC*DTCALCV),
c    &             DBLE(COUNTV(DTC+1))
c              WRITE(222,'(2ES16.8)') DT*DBLE(DTC*DTCALCD),
c    &             DBLE(COUNTD(DTC+1))
c           ENDDO
c           CLOSE(202)
c           CLOSE(203)
c           CLOSE(204)
c           CLOSE(205)
c           CLOSE(220)
c           CLOSE(221)
c           CLOSE(222)
c           OPEN(UNIT=223,FILE=FILE23,STATUS='REPLACE',ACTION='WRITE')
c               WRITE(223,'(1ES16.8)') summass
c           CLOSE(223)
ccc           call cpu_time(finishTime)
ccc           write(76,*) 'Time for write msd, vacf, dcf at step',
ccc    &                  istep0,'=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(86,*) 'Wall time for write msd, vacf, dcf',istep0,
ccc    &                  '=',elapsed_time
            ENDIF
         ENDIF
c        
c        add current positions to tape
c        
         IF(MOD(ISTEP,DTCALCM).EQ.0) THEN 
            TT=MOD(ISTEP/DTCALCM-1,ISTOREM)+1
            DO I=1,Nw,3
               K=wlist(I) ! Oxygen index
               XO(I,TT)=XC(K)
               YO(I,TT)=YC(K)
               ZO(I,TT)=ZC(K)
c              write(*,*) I,TT
c              write(*,*) XO(I,TT),YO(I,TT),ZO(I,TT)
               II=(I+2)/3
c         Add velocities to tape               
               J=Nw/3+II ! I (J) - Hydrogen 1 (2) 
               L=K+1 ! Hydrogen 1
               M=K+2 ! Hydrogen 2
               VXO(II,TT)=VX(K)
               VYO(II,TT)=VY(K)
               VZO(II,TT)=VZ(K)
               VXH(II,TT)=VX(L)
               VYH(II,TT)=VY(L)
               VZH(II,TT)=VZ(L)
               VXH(J,TT)=VX(M)
               VYH(J,TT)=VY(M)
               VZH(J,TT)=VZ(M)
            ENDDO
c         Add dipole moments to tape            
ccc            DO I=1,N
ccc               MXI(I,TT)=MXIC(I)
ccc               MYI(I,TT)=MYIC(I)
ccc               MZI(I,TT)=MZIC(I)
ccc            ENDDO
            MXI(TT)=MXC
            MYI(TT)=MYC
            MZI(TT)=MZC
         ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                                               c
c     Write out tape: All the tape values are written out, so that they can be read in          c
c     to use for restarting.                                                                    c
c                                                                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c        Store files for the Oxygen positions
c     
         IF(MOD(ISTEP,SAVE_REST_FREQ).EQ.0) THEN
ccc        call cpu_time(startTime)
ccc        CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
c           OPEN (UNIT=207,FILE=FILE7,STATUS='REPLACE') ! store file for X positions of Oxygens of waters
c           OPEN (UNIT=208,FILE=FILE8,STATUS='REPLACE') ! store file for Y positions of Oxygens of waters
c           OPEN (UNIT=209,FILE=FILE9,STATUS='REPLACE') ! store file for Z positions of Oxygens of waters
c           OPEN (UNIT=210,FILE=FILE10,STATUS='REPLACE') !store file for X velocities of Oxygens of waters
c           OPEN (UNIT=211,FILE=FILE11,STATUS='REPLACE') !store file for Y velocities of Oxygens of waters
c           OPEN (UNIT=212,FILE=FILE12,STATUS='REPLACE') !store file for Z velocities of Oxygens of waters
c           OPEN (UNIT=213,FILE=FILE13,STATUS='REPLACE') !store file for X velocities of Hydrogens of waters
c           OPEN (UNIT=214,FILE=FILE14,STATUS='REPLACE') !store file for Y velocities of Hydrogens of waters
c           OPEN (UNIT=215,FILE=FILE15,STATUS='REPLACE') !store file for Z velocities of Hydrogens of waters
c           OPEN (UNIT=216,FILE=FILE16,STATUS='REPLACE') !store file for X dipole moment
c           OPEN (UNIT=217,FILE=FILE17,STATUS='REPLACE') !store file for y dipole moment
c           OPEN (UNIT=218,FILE=FILE18,STATUS='REPLACE') !store file for Z dipole moment
c              IF (ISTEP.LT.ISTOREV*DTCALCV) THEN
c                 SAVE_LENGTH = ISTEP/DTCALCV
c              ELSE 
c                 SAVE_LENGTH = ISTOREV
c              END IF
c              DO DTC=1,SAVE_LENGTH
c                DO O_IDX=1,Nw
c                  WRITE(207,'(ES16.8,X)',advance='no') XO(O_IDX,DTC)
c                  WRITE(208,'(ES16.8,X)',advance='no') YO(O_IDX,DTC)
c                  WRITE(209,'(ES16.8,X)',advance='no') ZO(O_IDX,DTC)
c                ENDDO
c                DO O_IDX=1,Nw/3
c                  WRITE(210,'(ES16.8,X)',advance='no') VXO(O_IDX,DTC)
c                  WRITE(211,'(ES16.8,X)',advance='no') VYO(O_IDX,DTC)
c                  WRITE(212,'(ES16.8,X)',advance='no') VZO(O_IDX,DTC)
c                ENDDO
c                DO H_IDX=1,2*Nw/3
c                  WRITE(213,'(ES16.8,X)',advance='no') VXH(H_IDX,DTC)
c                  WRITE(214,'(ES16.8,X)',advance='no') VYH(H_IDX,DTC)
c                  WRITE(215,'(ES16.8,X)',advance='no') VZH(H_IDX,DTC)
c                ENDDO
c                WRITE(216,'(ES16.8)') MXI(DTC)
c                WRITE(217,'(ES16.8)') MYI(DTC)
c                WRITE(218,'(ES16.8)') MZI(DTC)
c                WRITE(207,'(/)',advance='no')
c                WRITE(208,'(/)',advance='no')
c                WRITE(209,'(/)',advance='no')
c                WRITE(210,'(/)',advance='no')
c                WRITE(211,'(/)',advance='no')
c                WRITE(212,'(/)',advance='no')
c                WRITE(213,'(/)',advance='no')
c                WRITE(214,'(/)',advance='no')
c                WRITE(215,'(/)',advance='no')
c              ENDDO
c           CLOSE(207)
c           CLOSE(208)
c           CLOSE(209)
c           CLOSE(210)
c           CLOSE(211)
c           CLOSE(212)
c           CLOSE(213)
c           CLOSE(214)
c           CLOSE(215)
c           CLOSE(216)
c           CLOSE(217)
c           CLOSE(218)
ccc           call cpu_time(finishTime)
ccc           write(77,*) 'Time for savefiles at step',istep0,
ccc    &                  '=',finishTime-startTime
ccc           CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
ccc           nb_ticks = nb_ticks_final - nb_ticks_initial
ccc           IF (nb_ticks_final < nb_ticks_initial) 
ccc    &                  nb_ticks = nb_ticks + nb_ticks_max
ccc           elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
ccc           write(87,*) 'Wall time for savefiles at step =',istep0,
ccc    &                  '=',elapsed_time
         ENDIF
      enddo
c     CLOSE(201) 
c     CLOSE(206)
ccc     call cpu_time(netFinishTime)
ccc     write(78,*) 'Time for whole program =',netFinishTime-netStartTime
ccc     call gettime(wall,cpu)            ! To get timings
ccc     print *, "Walltime is ", wall     ! To get timings
ccc     print *, "CPU time is ", cpu      ! To get timings
c
c     perform any final tasks before program exit
c
      write(*,*) "I made it to the end of dynamic."
      call final
      end
