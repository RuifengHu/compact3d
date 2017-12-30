c     FOURTH COMPACT DIFFERENCE DISCRETIZATION OF  
c     THE INCOMPRESSIBLE NAVIER-STOKES EQUATIONS
c     FOR THREE DIMENSIONAL PERIODIC TURBULENT CHANNEL FLOW
c
c     ORIGINALLY WRITTEN BY FLOW CONTROL LAB AT KAIST, KOREA
c     MODIFIED BY XIDIAN UNIVERSITY & LANZHOU UNIVERSITY
c
c     09/19/17 iwmles module for rough wall by rfhu
c     09/23/17 iwmles module for smooth wall by rfhu
                                                    
      PROGRAM MAIN
      INCLUDE 'COMMON1.FI'
      include 'omp_lib.h'  
      
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WMEANO/WMO(M2,3),WRMSO(M2,3)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)     
      COMMON/PARA/RE
      COMMON/FLOWRATE/FLOW1,FLOW3
c      COMMON/SIZE/ALX,ALY,ALZ,VOL
      
      COMMON/SPEC/XSPCTR(M1,M2,3),ZSPCTR(M3,M2,3)
      COMMON/SPECO/XSPCTRO(M1,M2,3),ZSPCTRO(M3,M2,3)
      COMMON/SPECTAU/SPECTAUW(M1,M3,2),SPECTAUB(M1,M3,2)
      COMMON/SPECTAUO/SPECTAUWO(M1,M3,2),SPECTAUBO(M1,M3,2)
      COMMON/CORR/RXX(M1,M2,3),RZZ(M3,M2,3)
      COMMON/CORRO/RXXO(M1,M2,3),RZZO(M3,M2,3)
      COMMON/CORR2/RXY(M1,M2,3)
      COMMON/CORR2O/RXYO(M1,M2,3)
      COMMON/QUAD/Q(M2,4),QP(M2,4)     
      COMMON/QUADO/QO(M2,4),QPO(M2,4)   

c-----for turbulent kinetic budgets 
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/PVMEANO/PVO(M2)
      COMMON/VHIGHO/VSKEWO(M2,3),VFLATO(M2,3),U2VO(M2),VW2O(M2)
      COMMON/PSTRO/PDUDXO(M2),PDVDYO(M2),PDWDZO(M2)
      COMMON/DISSUO/DUDX2O(M2),DUDY2O(M2),DUDZ2O(M2)
      COMMON/DISSVO/DVDX2O(M2),DVDY2O(M2),DVDZ2O(M2)
      COMMON/DISSWO/DWDX2O(M2),DWDY2O(M2),DWDZ2O(M2)      
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/WALLTAO/TAOW
      COMMON/WALLTAOAVE/TAOWO 
      COMMON/PARTICLE/DENRATIO,DIMP,PITIME,PDT     
      
      COMMON/ROUGHNESS/Y0
c      COMMON/SIZESCALE/HHH
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL IFWALL
      REAL HH(0:M1,0:M2,0:M3,3)  ! RIGHT-HAND SIDE CONVECTIVE TERM      
      REAL SS(0:M1,0:M2,0:M3,3)   ! RIGHT-HAND SIDE SUBGRID VISCOS TERM

      REAL UBULK,TBULK,TPLUS
      
!      nthds=num_parthds()
      nthds = omp_get_max_threads()
      
      tt0=rtc()      
      
      HHH=1000.0      
      Y0=0.0/HHH
      
      OPEN(31,FILE='WALLSS.plt',STATUS='UNKNOWN')
      OPEN(32,FILE='FLOW.plt',STATUS='UNKNOWN')
      OPEN(33,FILE='CFL.plt',STATUS='UNKNOWN')
      OPEN(34,FILE='RSH.plt',STATUS='UNKNOWN')
      OPEN(315,FILE='Retau.plt',STATUS='UNKNOWN')
      
      CALL SETUP  ! read inputs, mesh, index, etc.
      
      IF(NREAD.EQ.0) Call INIUP(U,P,HH,SS,PRESG,PRESG3)
      IF(NREAD.EQ.1) Call READUP(U,P,HH,SS,PRESG,PRESG3) 
      
      TIME=0.0
      TBULK=0.0
      
c      WRITE(*,*) DIVMAX      
      
      IMORE=0                ! Index of written file
      IPLUS=0                ! Index of written file
      IFWALL=0.0
      NAVG=1                
      DTR=DT
      NTIME=1                
      Init_P=1
      
      CALL DIVCHECK(U,TIME,DIVMAX,NTIME)  ! check divergence
      
      CALL CHKMF(U,TIME)  ! calculate flow rate  
      
!      DO 10 NTIME=1,NTST
20    t0=rtc()
      t1=mclock()
      t2 = omp_get_wtime()
      
      CALL CFL(U,CFLM) ! calc maximum cfl number
      DT=DTR
      IF (CFLM*DT.GE.CFLMAX.AND.IDTOPT.EQ.1) DT=CFLMAX/CFLM
c      IF (CFLM*DTR.LE.CFLMAX.OR.IDTOPT.NE.1)  DT=DTR
c      DT=DTR      
      IF(NTIME.LE.100.AND.NREAD.EQ.0) DT=0.001    
      
      TIME=TIME+DT      
      
      IF(IAVG.EQ.1.AND.NTIME.GT.50000.AND.MOD(NTIME,10).EQ.0) THEN  !.AND.NTIME.GE.10000
      NAVG=NAVG+1
      CALL ENERGY(U,P,NAVG)   ! x-z plane average of velocity, vorticity, pressure, etc       
      CALL SGSTRESS(U,SGSVIS)     !  x-z plane average of sgs stress 
      CALL TAVER(NAVG)    ! time-average of velocity, etc      
      ENDIF    
      
      IF(NTIME.GT.90000.AND.MOD(NTIME,200).EQ.0) THEN
      CALL WRITEUP(U,P,PRESG,PRESG3,NTIME) 
      ENDIF

      IF(WMODEL.EQ.0.) CALL WALLSS(U,P,TIME)    !  x-z plane average of wall shear stress
      IF(WMODEL.NE.0.) TAOWO=TAOW            
      
      ! IF(NTIME.Ge.500)CALL GET_PAR(U,NTIME) ! calculate the particle  
      
      CALL GETUP(U,P,HH,SS,NTIME,PRESG,PRESG3,SGSVIS)     
      
      CALL WALLSS_WRITE(U,P,PRESG,TIME)    
      
      CALL DIVCHECK(U,TIME,DIVMAX,NTIME) 
      
      CALL CHKMF(U,TIME)      
      
      CALL CFL(U,CFLM)      

      IF(WMODEL.EQ.0)THEN
      IF(IAVG.EQ.1.AND.NTIME.LE.50000) THEN
      UTAU=SQRT(WSM(1))
      ELSE
      UTAU=SQRT(WSMO(1))
      ENDIF
      ELSE
      UTAU=SQRT(TAOWO)
      ENDIF 

      UBULK=FLOW1/ALZ/ALY
      TBULK=TIME*UBULK/ALY
      TPLUS=TIME*UTAU
      
      
      IF(MOD(NTIME,100).EQ.0) THEN
      WRITE(*,100) NTIME,PRESG,DT,TBULK,TPLUS,DIVMAX,CFLM*DT,NAVG,
     >              UTAU,UTAU*RE
c      t0=rtc()-t0
c      t1=mclock()-t1
      t3 = omp_get_wtime()
      write(6,*) '-----------------------------------'
      write(6,*) '-----------------------------------'
      write(6,*) 'Elapsed time',(t3-t2)*100
c      write(6,*) 'CPU     time',t1*1.e-2
      END IF     
      
      
      WRITE(33,*) TIME,DT,CFLM*DT
      NTIME=NTIME+1
 
      IF (NTIME.GT.NTST) GOTO 10
      GOTO 20


 10   CONTINUE

      CALL WRITEUP(U,P,PRESG,PRESG3,NTIME) 
      CALL WRITEAVG(U,P,TIME,NAVG,IPLUS,NTIME)      
      CALL CWRITE(U,P,PRESG,PRESG3,NTIME,TIME)
      
C      IF (IAVG.EQ.1) CALL WRITEAVG(U,P,TIME,NAVG,IPLUS,NTIME)
C      CALL PROFILE(U,P,NTIME)

 100  FORMAT('STEP=',I6,2X,'DP=',E10.3,2X,'DT=',E12.5,2X,'TBULK=',E12.5,
     >        2X,'TPLUS=',E12.5,2X,'DIVMAX=',E10.3,2X,'CFL=',E10.3,2X,
     >        'NAVG=',I6.2,2X,'UTAU=',E10.3,2X,'RETAU=',E10.3)

      STOP
      END

c***************** SETUP ***********************     
      SUBROUTINE SETUP
      INCLUDE 'PARAM.H'
      
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/VPERIN/VPER
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/NEARWAL/NEARWALL
      !!!!!!!!!!!!!!!!!!
      
      Common/PAR1/ IFPAR,Nstep_P,Nprint_P,NTIMEP 
      COMMON/PARTICLE/DENRATIO,DIMP,PITIME,PDT,NUMPI
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*20 fileini,filegrd,fileout,fileavg
      CHARACTER*70 DUMMY
      integer NEARWALL
      
      OPEN(1,FILE='parame.ter',STATUS='OLD')
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,302) DUMMY,FILES
      WRITE(*,302) DUMMY,FILES
      READ (1,302) DUMMY,WMODEL
      WRITE(*,302) DUMMY,WMODEL
      READ (1,301) DUMMY,NEARWALL
      WRITE(*,301) DUMMY,NEARWALL
      READ (1,301) DUMMY,N1
      WRITE(*,301) DUMMY,N1
      READ (1,301) DUMMY,N2
      WRITE(*,301) DUMMY,N2
      READ (1,301) DUMMY,N3
      WRITE(*,301) DUMMY,N3
      READ (1,302) DUMMY,RE
      WRITE(*,302) DUMMY,RE
      READ (1,302) DUMMY,ALX
      WRITE(*,302) DUMMY,ALX
      READ (1,302) DUMMY,ALZ
      WRITE(*,302) DUMMY,ALZ
      READ (1,301) DUMMY,INCODE
      WRITE(*,301) DUMMY,INCODE
      READ (1,301) DUMMY,NTST
      WRITE(*,301) DUMMY,NTST
      READ (1,302) DUMMY,VPER
      WRITE(*,302) DUMMY,VPER
      READ (1,302) DUMMY,DT
      WRITE(*,302) DUMMY,DT
      READ (1,301) DUMMY,IDTOPT
      WRITE(*,301) DUMMY,IDTOPT
      READ (1,302) DUMMY,CFLMAX
      WRITE(*,302) DUMMY,CFLMAX
      READ (1,301) DUMMY,NWRITE
      WRITE(*,301) DUMMY,NWRITE
      READ (1,301) DUMMY,NREAD
      WRITE(*,301) DUMMY,NREAD
      READ (1,301) DUMMY,IAVG
      WRITE(*,301) DUMMY,IAVG
      READ (1,301) DUMMY,NPRN
      WRITE(*,301) DUMMY,NPRN
      READ (1,301) DUMMY,INSF
      WRITE(*,301) DUMMY,INSF
      READ (1,301) DUMMY,NINS
      WRITE(*,301) DUMMY,NINS
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,303) DUMMY,fileini
      WRITE(*,303) DUMMY,fileini
      READ (1,303) DUMMY,filegrd
      WRITE(*,303) DUMMY,filegrd
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      READ (1,300) DUMMY
      WRITE(*,300) DUMMY
      

  300 FORMAT(A65)    
  301 FORMAT(A45,I15)
  302 FORMAT(A45,E15.7)
  303 FORMAT(A45,A20)

C------------------------------------
C     PHYSICAL LENGTH
      IF(INCODE.EQ.1) THEN
      PI=ACOS(-1.0)
      ALX=4.0*PI     
      ALZ=2.0*PI
      ALY=1.0     
      ENDIF
C------------------------------------
      
      N1M=N1-1
      N2M=N2-1
      N3M=N3-1

      CALL MESH
      CALL INDICES 
      CALL INIWAVE
     
      RETURN
      END

c***************** MESH ***********************     
      SUBROUTINE MESH
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/MESH4/DYM(M2),DYC(M2),DYP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLCO/HPW(M2),HMW(M2),HCW(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/ROUGHNESS/Y0
      
      CHARACTER*20 fileini,filegrd,fileout,fileavg

c      CREATE THE UNIFORM GRID IN X2 DIRECTION      
c!$omp parallel do
      DO I=1,N1M
      XG(I)=DBLE(I-1)*ALX/DBLE(M1-1)
      ENDDO
      XG(0)=0.0
      XG(M1)=ALX
      
c!$omp parallel do      
      DO K=1,N3M
      ZG(K)=DBLE(K-1)*ALZ/DBLE(M3-1)
      ENDDO
      ZG(0)=0.0
      ZG(M3)=ALZ 
      
CCCCC    
      
      CKS=2.5      
      IF (WMODEL.EQ.0.0) THEN
      Y(0)=0.0
      YC(0)=0.0
      YG(1)=0.0
      ELSE
      Y(0)=Y0
      YC(0)=Y0
      YG(1)=Y0
      ENDIF      

      DO J=1,N2     
      IF (WMODEL.EQ.0.0) THEN 
      IF (ALY.EQ.2.) THEN
      Y(J)=1-TANH(CKS*(1-ALY*(J-1)/(N2-1)))/TANH(CKS)
	ELSE
	Y(J)=1-TANH(CKS*(1-2.*ALY*(J-1)/(2*N2-1)))/TANH(CKS)
	ENDIF
      ELSE
      Y(J)=(ALY-Y0)*(J-1)/(N2-1)+Y0
      ENDIF
      YG(J)=Y(J)
      ENDDO      
      ALY = Y(N2)

c!$omp parallel do       
      DO I=1,N1M
      XC(I)=0.5*(XG(I+1)+XG(I))
      ENDDO
      XC(0)=0.0
      XC(M1)=XG(M1)  
   
c!$omp parallel do
      DO K=1,N3M
      ZC(K)=0.5*(ZG(K+1)+ZG(K))
      ENDDO
      ZC(0)=0.0
      ZC(M3)=ZG(M3)

c!$omp parallel do      
      DO J=1,N2M
      YC(J)=0.5*(Y(J+1)+Y(J))
      ENDDO
      YC(N2)=Y(N2)
      
      VOL=ALX*ALY*ALZ 
      DX1=DBLE(N1M)/ALX
      DX3=DBLE(N3M)/ALZ
      DX1Q=DX1**2.0
      DX3Q=DX3**2.0
      
      DY(1)=Y(2)
c!$omp parallel do
      DO 20 J=2,N2M
      DY(J)=Y(J+1)-Y(J)
      H(J)=0.5*(DY(J)+DY(J-1))
 20   CONTINUE
      H(1)=0.5*DY(1)
      H(N2)=0.5*DY(N2M)

c!$omp parallel do      
      DO 30 J=2,N2M-1
      HP(J)=1.0/H(J+1)/DY(J)
      HC(J)=(H(J+1)+H(J))/H(J+1)/H(J)/DY(J)
      HM(J)=1.0/H(J)/DY(J)
 30   CONTINUE 
 
!     AKSELVOLL&MOIN(1995)'S SPECIAL TREATMENT AT THE WALL.
!     SEE THE PAGE 33-34 OF REPORT NO. TF-63 (1995)
      HP(1)=4.0*H(1)/H(2)/(H(1)+H(2))/DY(1)
      HC(1)=4.0/H(2)/DY(1)
      HM(1)=4.0/(H(1)+H(2))/DY(1)
      HP(N2M)=4.0/(H(N2M)+H(N2))/DY(N2M)
      HC(N2M)=4.0/H(N2M)/DY(N2M)
      HM(N2M)=4.0*H(N2)/H(N2M)/(H(N2M)+H(N2))/DY(N2M)    
      

      IF(WMODEL.NE.0.)THEN
      HP(1)=1.0/H(2)/DY(1)
      HC(1)=1.0/H(2)/DY(1)
      HM(1)=0.0
      HP(N2M)=0.0
      HC(N2M)=1.0/H(N2M)/DY(N2M)
      HM(N2M)=1.0/H(N2M)/DY(N2M)
      ENDIF    
     
      DO 40 J=2,N2M
      DYP(J)=1.0/H(J)/DY(J)
      DYC(J)=1.0/H(J)*(1.0/DY(J)+1.0/DY(J-1))
      DYM(J)=1.0/H(J)/DY(J-1)
 40   CONTINUE
 
  
      OPEN(203,FILE='GRID.PLT',STATUS='UNKNOWN')
      DO J=0,N2
      IF(J.GT.0.AND.J.LT.N2)THEN
      DETAY=Y(J+1)-Y(J)
      ELSE 
      DETAY=0.0
      ENDIF
      WRITE(203,*)j,Y(J),DETAY
      END DO
      CLOSE(203)     
      
      RETURN
      END

c***************** INDICES ***********************     
      SUBROUTINE INDICES
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINDX3/FUP(M2),FDOWN(M2)

c!$omp parallel do
      DO 10 IC=1,N1M
      IPA(IC)=IC+1   ! I+1
 10   IMA(IC)=IC-1   ! I-1
      IPA(N1M)=1
      IMA(1)=N1M

c!$omp parallel do
      DO 20 KC=1,N3M
      KPA(KC)=KC+1
 20   KMA(KC)=KC-1
      KPA(N3M)=1
      KMA(1)=N3M

c!$omp parallel do
      DO 30 JC=1,N2M
      JPA(JC)=JC+1
      JMU(JC)=JC-1
 30   JMV(JC)=JC-1
      JPA(N2M)=N2M
      JMU(1)=1
      JMV(2)=2

c!$omp parallel do      
      DO 40 JC=1,N2M
      FJPA(JC)=JPA(JC)-JC   ! At boundary, these values are "0".
      FJMU(JC)=JC-JMU(JC)   ! Otherwise, these are "1".
 40   FJMV(JC)=JC-JMV(JC)

c!$omp parallel do                            
      DO 50 JC=1,N2M/2
      FDOWN(JC)=1
 50   FUP  (JC)=0

c!$omp parallel do 
      DO 60 JC=N2M/2+1,N2M
      FDOWN(JC)=0
 60   FUP  (JC)=1
                            
      RETURN
      END

C  ****************************** INIUP **********************
C     THIS ROUTINE IS TO GIVE INITIAL FLOW FIELDS 
C     ON A PARABOLIC PROFILES
C     THIS SUBROUTINE IS MODIFIED BY JICHOI 99/07/27
C     TO MAKE CONSTANT FLOW RATE AT EACH PLANE
C     U IN YZ-PLANE, W IN XY-PLANE

      SUBROUTINE INIUP(U,P,HH,SS,PRESG,PRESG3)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/PARA/RE
      COMMON/VPERIN/VPER
      COMMON/SIZE/ALX,ALY,ALZ,VOL

      REAL U(0:M1,0:M2,0:M3,3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL RDNUM(3)
      INTEGER ISEED(4)
      REAL NXZ

      ISEED(1)=1001
      ISEED(2)=2001
      ISEED(3)=3001
      ISEED(4)=4001

C     Random Number Generation in MKL
C     Math Kernel Library Manual pp.5-182 ~5-183
C     call dlarnv(idist, iseed, n, x)
C     idist INTEGER.iseed On exit, the seed is updated.
C          =1:uniform(0,1)
C          =2:uniform(-1,1)
C          =3:normal(0,1)
C     iseed INTEGER, Array, DIMENSION(4)
C           On entry, the seed of the random number generator;     
C           the array elements
C           must be between 0 and 4095, and iseed(4) must be odd.
C     n integer. The number of random numbers to be generated.
C     x DOUBLE PRECISION for dlarnv, Array, DIMENSION(n)
C           The generated random numbers.
C     initialization U
!      U=0.
!      PI=ACOS(-1.) 
!      DO 10 K=1,N3M
!      DO 10 J=1,N2M
!      DO 10 I=1,N1M
!      CALL DLARNV(2,ISEED,3,RDNUM)
!      U(I,J,K,1)=VPER*RDNUM(1)
!      U(I,J,K,2)=VPER*RDNUM(2)
!      U(I,J,K,3)=VPER*RDNUM(3)
!10    CONTINUE

      call random_seed ()           
      U=0.
      PI=ACOS(-1.) 
!$omp parallel do private(ru11,ru21,ru31,ru1,ru2,ru3)
      DO 10 I=1,N1M
      DO 10 J=1,N2M
      DO 10 K=1,N3M
      call random_number (ru11)
      ru1=-1.0+2.0*ru11
      U(I,J,K,1)=ru1 !*VPER
      call random_number (ru21)
      ru2=-1.0+2.0*ru21
      U(I,J,K,2)=ru2 !*VPER
      call random_number (ru31)
      ru3=-1.0+2.0*ru31
      U(I,J,K,3)=ru3 !*VPER 
10    CONTINUE



C     ELIMINATE MEAN QUANTITIES OF RANDOM FLUCTUATIONS
!$omp parallel do private(V1M)
      DO 21 I=1,N1M
      V1M=0.
      DO 22 K=1,N3M
      DO 22 J=1,N2M
      V1M=V1M+U(I,J,K,1)*DY(J)/DX3
22    CONTINUE
      V1M=V1M/ALY/ALZ
      DO 23 K=1,N3M
      DO 23 J=1,N2M
      U(I,J,K,1)=U(I,J,K,1)-V1M
23    CONTINUE
21    CONTINUE

!$omp parallel do private(V2M)
      DO 31 J=2,N2M
      V2M=0.
      DO 32 K=1,N3M
      DO 32 I=1,N1M
      V2M=V2M+U(I,J,K,2)/DX1/DX3
32    CONTINUE
      V2M=V2M/ALX/ALZ
      DO 33 K=1,N3M
      DO 33 I=1,N1M
      U(I,J,K,2)=U(I,J,K,2)-V2M
33    CONTINUE
31    CONTINUE

!$omp parallel do private(V3M)
      DO 41 K=1,N3M
      V3M=0.
      DO 42 J=1,N2M
      DO 42 I=1,N1M
      V3M=V3M+U(I,J,K,3)*DY(J)/DX1
42    CONTINUE
      V3M=V3M/ALX/ALY
      DO 43 J=1,N2M
      DO 43 I=1,N1M
      U(I,J,K,3)=U(I,J,K,3)-V3M
43    CONTINUE
41    CONTINUE

C     IMPOSE LAMINAR VELOCITY PROFILES IN U VELOCITIESA

      RET=0.09*RE**0.88  ! FROM REC TO RE_TAU
      UTAU_INI=RET/RE

!$omp parallel do private (YH,JP,REALU)
      DO 30 K=1,N3M
      DO 30 J=1,N2M
      JP=J+1
      DO 30 I=1,N1M
      YH=0.5*(Y(J)+Y(JP))
      REALU=YH*(ALY-YH)*(ALY-1.0)+YH*(2*ALY-YH)*(2.0-ALY)  ! laminar flow, full channel 0828 (NON-DIMENSIONAL)
! !     REALU=   !laminar flow, half channel 0828
!      IF(J.LE.N2M/2)THEN
!      REALU=1.0-UTAU_INI*(-1.0/0.41*ALOG(YH)+0.2)  !log law
!      ELSE
!      REALU=1.0-UTAU_INI*(-1.0/0.41*ALOG((ALY-YH))+0.2) !log law
!      ENDIF      

      !U(I,J,K,1)=U(I,J,K,1)+REALU

      U(I,J,K,1)=U(I,J,K,1)*VPER*REALU+REALU
      U(I,J,K,2)=U(I,J,K,2)*VPER   
      U(I,J,K,3)=U(I,J,K,3)*VPER
30    CONTINUE


C     CHECK FLOW RATE
      FLOW1=0.0
c      FLOW3=0.0
!$omp parallel do reduction(+:FLOW1)
      DO 40 K=1,N3M
      DO 40 J=1,N2M
      DO 40 I=1,N1M
      FLOW1=FLOW1+U(I,J,K,1)*DY(J)/DX1/DX3
c      FLOW3=FLOW3+U(I,J,K,3)*DY(J)/DX1/DX3
40    CONTINUE
      FLOW1=FLOW1/ALX/ALZ
c      FLOW3=FLOW3/ALX/ALZ/ALY

      FLOWR=(2*ALY)**3.0/6.0/2.0*(2.0-ALY)+(ALY)**3.0/6.0*(ALY-1.0) ! LAMINAR MASS FLOW RATE IN X1 DIRECTION
!$omp parallel do
      DO 45 K=1,N3M
      DO 45 J=1,N2M
      DO 45 I=1,N1M
      U(I,J,K,1)=FLOWR/FLOW1*U(I,J,K,1)
      U(I,J,K,3)=FLOWR/FLOW1*U(I,J,K,3)
   45 CONTINUE

C     IMPOSE ZERO VELOCITY AT BOUNDARY
!$omp parallel do
      DO 15 K=1,N3M  
      DO 15 I=1,N1M
      U(I,0,K,1)=0.0
      U(I,N2,K,1)=U(I,N2M,K,1)*(2.0-ALY)
      U(I,1,K,2)=0.0
      U(I,N2,K,2)=0.0
      U(I,0,K,3)=0.0
      U(I,N2,K,3)=U(I,N2M,K,3)*(2.0-ALY)
15    CONTINUE

C     IMPOSE ZERO-PRESSURE FLUCTUATIONS
!$omp parallel do
      DO 60 K=1,N3M
      DO 60 J=1,N2M
      DO 60 I=1,N1M
      P(I,J,K)=0.0
   60 CONTINUE
      
C     IMPOSE ZERO HH

      DO 65 NV=1,3
!$omp parallel do     
      DO 65 K=1,N3M
      DO 65 J=1,N2M
      DO 65 I=1,N1M
      HH(I,J,K,NV)=0.0
      SS(I,J,K,NV)=0.0
   65 CONTINUE
      
C     INITIAL MEAN PRESSURE GRADIENT AT LAMINAR FLOW FIELD
      PRESG=0.0  
      PRESG3=0.0      ! Z-direction pressure gradient 

      OPEN(201,FILE='INIU.PLT')
      WRITE(201,*)'VARIABLES="X","Y","Z","U","V","W"'
      WRITE(201,*)'ZONE T="',1,'" I=',N1M ,' J=',N2+1,' K=',N3M
      DO K=1,N3M
      DO J=0,N2
      DO I=1,N1M
      IF(J.NE.N2)THEN
      WRITE(201,333) (I-1)/DX1,0.5*(Y(J)+Y(J+1)),(K-1)/DX3,
     >            U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
      ENDIF
      IF(J.EQ.N2)THEN
      WRITE(201,333) (I-1)/DX1,Y(J),(K-1)/DX3,
     >            U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
333   FORMAT(6(E12.5,2X))
      CLOSE(201)
 
      
      
      RETURN
      END

C  ************************  READUP **********************
C     READ FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE READUP(U,P,HH,SS,PRESG,PRESG3)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*20 fileini,filegrd,fileout,fileavg

      REAL U(0:M1,0:M2,0:M3,3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)

      OPEN(3,FILE='CONTINUE',FORM='UNFORMATTED',STATUS='OLD')
      READ(3) PRESG1,PRESG3
      READ(3) (((U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
     >            ,K=1,N3M),J=0,N2),I=1,N1M)
      READ(3) (((P(I,J,K),K=1,N3M),J=1,N2M),I=1,N1M)
      CLOSE(3)

      RETURN
      END



C  ************************  WRITEUP **********************
C     WRITE FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE WRITEUP(U,P,PRESG,PRESG3,NTIME)
      INCLUDE 'PARAM.H'
      
	COMMON/MESH2/Y(0:M2)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
     
      COMMON/PARA/RE    
     
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/SIZESCALE/HHH
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/TAUXYZO/TAUWO(0:M1,0:M3,2),TAUBO(0:M1,0:M3,2)
      COMMON/WALLTAO/TAOW
      COMMON/WALLTAOAVE/TAOWO
      COMMON/WALLMODEL/WMODEL
      
      CHARACTER(6) FILENAME
      REAL U(0:M1,0:M2,0:M3,3),TAUW0(0:M1,0:M3)
      REAL P(M1,M2,M3)
      REAL UF,VF,WF,TERMP
      REAL UTAU,TAUDOWN
      WRITE(FILENAME,'(I6)') 100000+NTIME
      
      IF(WMODEL.EQ.0)THEN
      UTAU=SQRT(WSMO(1))
      ELSE
      UTAU=SQRT(TAOWO)
      ENDIF      

      Open (1,file='RESULT/FLOW-'//FILENAME//'.plt')	  
      WRITE(1,*)'VARIABLES="x","y","z","U","V","W","P" '
      Write(1,*)'ZONE T="',NTIME,'"I= ',N1M,'J=',N2M,'K=',N3M
      TERMP=1.0 !RE*1.5/100000.0/HHH
      DO K=1,N3M
      KP=KPA(K)
      DO J=1,N2M
      JP=JPA(J)
      DO I=1,N1M
      IP=IPA(I)
	  
      UF=0.5*(U(I,J,K,1)+U(IP,J,K,1))
      VF=0.5*(U(I,J,K,2)+U(I,JP,K,2))
      WF=0.5*(U(I,J,K,3)+U(I,J,KP,3))
      UF=UF*TERMP
      VF=VF*TERMP
      WF=WF*TERMP
      PP=P(I,J,K)*TERMP
	Write(1,'(7F20.10)')XC(I),YC(J),ZC(K),UF,VF,WF,PP
      Enddo
      ENDDO  
      ENDDO
      Close(1)
      
      TAUDOWN=0.
      DO K=1,N3M
      DO I=1,N1M   
      TAUDOWN=TAUDOWN+TAUW(I,K,1)
      ENDDO
      ENDDO
      TAUDOWN=TAUDOWN/DBLE(N1M*N3M)
      
      Open (2,file='RESULT/STRAIN-'//FILENAME//'.plt')	  
      WRITE(2,*)'VARIABLES="x","z","strain"'
      Write(2,*)'ZONE T="',NTIME,'"I= ',N1M,'K=',N3M
      DO K=1,N3M
      DO I=1,N1M
      TAUW0(I,K)=TAUW(I,K,1)-TAUDOWN
      WRITE(2,'(7F20.10)')  XC(I),ZC(K),TAUW0(I,K)/(UTAU**2)
      ENDDO
      ENDDO
      Close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END 
!------------------------------------------------------------------------------
      
C  ************************CWRITE*************************
C     WRITE FLOW FIELD AND BOUNDARY CONDITIONS
C     AND MEAN PRESSURE GRADIENT

      SUBROUTINE CWRITE(U,P,PRESG,PRESG3,NTIME,TIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/FILENAME/fileini,filegrd,fileout,fileavg
      CHARACTER*20 fileini,filegrd,fileout,fileavg

      REAL U(0:M1,0:M2,0:M3,3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)

      OPEN(3,FILE='CONTINUE')
      WRITE(3,*) PRESG,PRESG3
      WRITE(3,*) NTIME
      WRITE(3,*) TIME
      DO K=1,N3M
      DO J=0,N2
      DO I=1,N1M
      WRITE(3,*) U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)
      ENDDO
      ENDDO
      ENDDO
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      WRITE(3,*) P(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(3)

      RETURN
      END   


c***************** GETUP ************************************************     
      SUBROUTINE GETUP(U,P,HH,SS,NTIME,PRESG,PRESG3,SGSVIS)
      INCLUDE 'COMMON1.FI'
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL HH2(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL SS2(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL DP(M1,M2,M3)
      REAL SGSH(M2)
      
C INITIALIZE THE INTERMEDIATE VELOCITY AND PRESSURE

      DO 10 NV=1,3  
!$omp parallel do      
      DO 10 K=1,N3
      DO 10 J=1,N2
      DO 10 I=1,N1
 10   UH(I,J,K,NV)=0.0
      
!$omp parallel do
      DO 20 K=1,N3
      DO 20 J=1,N2
      DO 20 I=1,N1
 20   DP(I,J,K)=0.0
      
!$omp parallel do
      DO 30 K=1,N3
      DO 30 J=0,N2
      DO 30 I=1,N1
 30   SGSVIS(I,J,K)=0.0
      

      DO 40 NV=1,3  
!$omp parallel do      
      DO 40 K=1,N3
      DO 40 J=1,N2
      DO 40 I=1,N1
      HH2(I,J,K,NV)=HH(I,J,K,NV)
      SS2(I,J,K,NV)=SS(I,J,K,NV)
      HH(I,J,K,NV)=0.0
 40   SS(I,J,K,NV)=0.0 

      IF(FILES.EQ.1.0) CALL SGSDY(U,SGSVIS,NTIME)  ! DYNAMIC MODEL
      IF(FILES.EQ.2.0) CALL SGSCS(U,SGSVIS,NTIME)  ! SMAGRINSKY MODEL
      IF(FILES.EQ.3.0) CALL SGSSD(U,SGSVIS,NTIME)  ! SCALE-DEPENDENT MODEL
       
      IF(FILES.NE.0.0) CALL WALL_MODEL(U,P,PRESG,PRESG3,SGSVIS,NTIME)   ! WALL MODEL
C      IF(FILES.NE.0.0) CALL NEARWALLTREAD(U,SGSVIS,NTIME)  
       
      CALL BCOND(U,SGSVIS,NTIME) ! lower and upper boundary condition
      
C CALCULATE THE INTERMEDIATE VELOCITY
      CALL UHCALC(U,UH,P,HH,HH2,SS,SS2,PRESG,PRESG3,SGSVIS,NTIME)
      
C CALCULATE DP
      CALL DPCALC(U,UH,DP,NTIME)
      
C UPDATE THE N+1 TIME STEP VELOCITY AND PRESSURE
      CALL UPCALC(U,P,UH,DP,PRESG,PRESG3,NTIME)     
      
      
      IF(MOD(NTIME,100) .EQ.0.0 .AND.1.EQ.1)THEN 	
      OPEN(1,FILE='SGSVIS-H.plt')       
      DO J=1,N2M
      SGSH(J)=0.0
      DO K=1,N3M
      DO I=1,N1M
      SGSH(J)=SGSH(J)+SGSVIS(I,J,K)
      ENDDO
      ENDDO
      SGSH(J)=SGSH(J)/DBLE(N1M*N3M)
      WRITE(1,*) YC(J),SGSH(J)
      ENDDO
      CLOSE(1)      
      ENDIF

      RETURN
      END

c***************** BCOND **************************************     
      SUBROUTINE BCOND(U,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/BCON/UBC(0:M1,0:M3,3,2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      REAL U(0:M1,0:M2,0:M3,3)
      REAL SGSVIS(0:M1,0:M2,0:M3)

!$omp parallel do
      DO 10 K=1,N3M
      DO 10 I=1,N1M
      
      UBC(I,K,1,1)=0.0    ! lower wall : no-slip condition
      UBC(I,K,2,1)=0.0
      UBC(I,K,3,1)=0.0
      
      U(I,0,K,1)=UBC(I,K,1,1)
      U(I,0,K,2)=UBC(I,K,2,1)
      U(I,0,K,3)=UBC(I,K,3,1)
      
      UBC(I,K,1,2)=U(I,N2M,K,1)*(2.0-ALY)  ! upper wall condition (ALY=2: no-slip; ALY=1: free-slip)
      UBC(I,K,2,2)=0.0 
c      UBC(I,K,2,2)=U(I,N2M,K,2)*(2.0-ALY)
      UBC(I,K,3,2)=U(I,N2M,K,3)*(2.0-ALY)
      
      U(I,N2,K,1)=UBC(I,K,1,2)
      U(I,N2,K,2)=UBC(I,K,2,2)
      U(I,N2,K,3)=UBC(I,K,3,2)
      
      IF(FILES.NE.0.0) THEN
      SGSVIS(I,0,K)=0.0
      SGSVIS(I,N2,K)=SGSVIS(I,N2M,K)*(2.0-ALY)
      ENDIF
      
 10   CONTINUE
      RETURN
      END

c***************** UHCALC *************************************     
      SUBROUTINE UHCALC(U,UH,P,HH,HH2,SS,SS2,PRESG,PRESG3,SGSVIS,NTIME)
      INCLUDE 'COMMON1.FI'    
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL HH2(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL SS2(0:M1,0:M2,0:M3,3)

      REAL UH(0:M1,0:M2,0:M3,3)
      REAL RUH1(0:M1,0:M2,0:M3)
      REAL RUH2(0:M1,0:M2,0:M3)
      REAL RUH3(0:M1,0:M2,0:M3)
      SAVE RUH1,RUH2,RUH3
      
      ! calculate right-hand side of momemtum equations
      CALL RHS1(U,P,HH,HH2,SS,SS2,RUH1,PRESG,SGSVIS,NTIME)
      CALL RHS2(U,P,HH,HH2,SS,SS2,RUH2,SGSVIS,NTIME)
      CALL RHS3(U,P,HH,HH2,SS,SS2,RUH3,PRESG3,SGSVIS,NTIME)

      !  CALL SOURCE_PAR(U,RUH1,RUH2,RUH3,NTIME)
      
      ! calculate intermediate velocity
      CALL GETUH1(U,UH,RUH1,SGSVIS,NTIME)
      CALL GETUH2(U,UH,RUH2,SGSVIS,NTIME)
      CALL GETUH3(U,UH,RUH3,SGSVIS,NTIME)

      RETURN
      END

c***************** DPCALC ***********************     
      SUBROUTINE DPCALC(U,UH,DP,NTIME)
      INCLUDE 'PARAM.H'
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL DP(M1,M2,M3)
      REAL RDP(M1,M2,M3)

      CALL RHSDP(RDP,UH,U,NTIME)
      CALL GETDP(DP,RDP,NTIME)

      RETURN
      END

c***************** UPCALC ***********************     
      SUBROUTINE UPCALC(U,P,UH,DP,PRESG,PRESG3,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/BCON/UBC(0:M1,0:M3,3,2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/SIZE/ALX,ALY,ALZ,VOL

      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL DP(M1,M2,M3)
      REAL DMPRESG
      REAL DMPRESG3
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1),X1(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1),X3(M3-1)
      
C     U1 VELOCITY UPDATE
C     To keep a constant mass flow rate,
C     calculate the mean pressure gradient difference
C     between n+1/2 time step And n-1/2 time step 
      
      DMPRESG=0.0 
      FLOW1=0.0
!$omp parallel do private(IM) reduction(+:FLOW1,DMPRESG)
      DO 1 K=1,N3M
      DO 1 J=1,N2M
      DO 1 I=1,N1M
      IM=IMA(I)
      FLOW1=FLOW1+U(I,J,K,1)*DY(J)/DX3/DX1
      DMPRESG=DMPRESG
     >       +UH(I,J,K,1)*DY(J)/DX3/DX1
     >       -DT*(DP(I,J,K)-DP(IM,J,K))*DY(J)/DX3
  1   CONTINUE      
C      WRITE(*,*)DMPRESG,FLOW1
       
      DMPRESG=(DMPRESG-FLOW1)/VOL/DT   
     
cc!$omp parallel do private(IM)
c      DO 10 K=1,N3M
c      DO 10 J=1,N2M
c      DO 10 I=1,N1M
c      IM=IMA(I)
c      U(I,J,K,1)=UH(I,J,K,1)
c     >          -DT*(DP(I,J,K)-DP(IM,J,K))*DX1
c     >          -DT*DMPRESG
c      write(*,*) I,J,K,U(I,J,K,1)
c  10  CONTINUE
     
!$omp parallel do
      DO 35 I=1,N1M
      A1(I)=1.0
      B1(I)=22.0
      C1(I)=1.0
 35   CONTINUE 
     
!$omp parallel do private (D1,X1)
      DO 40 K=1,N3M
      DO 40 J=1,N2M     
      
      DO I=2,N1M
      D1(I)=24.*DX1*(-DP(I-1,J,K)+DP(I,J,K))
      ENDDO      
      D1(1)=24.*DX1*(-DP(N1M,J,K)+DP(1,J,K))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      U(I,J,K,1)=UH(I,J,K,1)-DT*X1(I)-DT*DMPRESG    
      ENDDO      
      
  40  CONTINUE  
       
!$omp parallel do 
      DO 11 K=1,N3M
      DO 11 I=1,N1M
      U(I,0,K,1) =UBC(I,K,1,1)
      U(I,N2,K,1)=UBC(I,K,1,2)
  11  CONTINUE
      
      
c     U2 VELOCITY UPDATE
!$omp parallel do 
      DO 20 K=1,N3M
      DO 20 J=2,N2M
      DO 20 I=1,N1M
      U(I,J,K,2)=UH(I,J,K,2)
     >          -DT*(DP(I,J,K)-DP(I,J-1,K))/H(J)
  20  CONTINUE
       
!$omp parallel do 
      DO 21 K=1,N3M
      DO 21 I=1,N1M
      U(I,1,K,2) =UBC(I,K,2,1)
      U(I,N2,K,2)=UBC(I,K,2,2)
  21  CONTINUE
      
      
C     U3 VELOCITY UPDATE
C     To make the mass flow rate in x3 direction zero
C     calculate the mean pressure gradient in x3 direction
C     this term may be zero at quasi-steady state

      DMPRESG3=0.0 
      FLOW3=0.0
!$omp parallel do private(KM) reduction(+:FLOW3,DMPRESG3)
      DO 3 K=1,N3M
      KM=KMA(K)
      DO 3 J=1,N2M
      DO 3 I=1,N1M
      FLOW3=FLOW3+U(I,J,K,3)*DY(J)/DX1/DX3
      DMPRESG3=DMPRESG3
     >       +UH(I,J,K,3)*DY(J)/DX3/DX1
     >       -DT*(DP(I,J,K)-DP(I,J,KM))*DY(J)/DX1
  3   CONTINUE
      DMPRESG3=(DMPRESG3-FLOW3)/VOL/DT 
      
C      WRITE(*,*) DMPRESG3

C      c!$omp parallel do private(KM)
C      DO 30 K=1,N3M
C      KM=KMA(K)
C      DO 30 J=1,N2M
C      DO 30 I=1,N1M
C      U(I,J,K,3)=UH(I,J,K,3)
C     >          -DT*(DP(I,J,K)-DP(I,J,KM))*DX3
C     >          -DT*DMPRESG3 
C  30  CONTINUE

!$omp parallel do
      DO 45 K=1,N3M
      A3(K)=1.0
      B3(K)=22.0
      C3(K)=1.0
 45   CONTINUE 
     
!$omp parallel do private (D3,X3)
      DO 50 I=1,N1M
      DO 50 J=1,N2M    
      
      DO K=2,N3M
      D3(K)=24.*DX3*(-DP(I,J,K-1)+DP(I,J,K))
      ENDDO      
      D3(1)=24.*DX3*(-DP(I,J,N3M)+DP(I,J,1))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      U(I,J,K,3)=UH(I,J,K,3)-DT*X3(K)-DT*DMPRESG3  
      ENDDO      
      
  50  CONTINUE 
       
!$omp parallel do
      DO 31 K=1,N3M
      DO 31 I=1,N1M
      U(I,0,K,3) =UBC(I,K,3,1)
      U(I,N2,K,3)=UBC(I,K,3,2)
  31  CONTINUE
      
C     PRESSURE UPDATE 
!$omp parallel do private (XP1,XP3)
      DO 60 K=1,N3M
      DO 60 J=1,N2M
      DO 60 I=1,N1M
C      XP1=(REAL(I-1)+0.5)/DX1
C      XP3=(REAL(K-1)+0.5)/DX3
      P(I,J,K)=P(I,J,K)+DP(I,J,K)
c      P(I,J,K)=P(I,J,K)+DP(I,J,K)+XP1*DMPRESG+XP3*DMPRESG3
c      IF (I.EQ.2 .AND. J.EQ.2 .AND. K.EQ.2 .AND.
c     >        MOD(NTIME,100).EQ.0) THEN
c      WRITE(*,*) U(I,J,K,1),U(I,J,K,2),U(I,J,K,3)    
c      ENDIF
  60  CONTINUE

      PRESG=PRESG+DMPRESG
      PRESG3=PRESG3+DMPRESG3
      
C      write(*,*) NTIME,PRESG,PRESG3
      
      RETURN
      END

C************************** SGSDS *****************************
C     THIS SUBROUTINE IS TO OBTAIN SUBGRID-SCALE VISCOSITY. 
C     VALUES ARE EVALUATED AT THE CELL CENTER.
C     #(I,K,1) = #_{11}
C     #(I,K,2) = #_{12}
C     #(I,K,3) = #_{13}
C     #(I,K,4) = #_{22}
C     #(I,K,5) = #_{23}

C     DYNAMIC SMAGORINSKY MODEL
      SUBROUTINE SGSDY(U,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/SGSMEAN/STRM(M2,6),SGSNUT(M2),SGSDS(M2),SMAGC(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)  
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/RANSVISO/RANSVIS,RANSVIS2

      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL COM(M1,M3,6),SR(M1,M3,6),UC(M1,M3,3)
      REAL NXZ,CS_DAMPING
      REAL SMAGC1(M2),SGS1(M2)
      
      REAL XX(5),YY(5),XA,XB
      NXZ=1.0/DBLE(N1M*N3M)
      RENINV=1./DBLE(RE)



      DO 1000 J=1,N2M
      
      DXYZ=1.0/DX1*DY(J)/DX3
 
      CALL STRAIN(U,SR,J) 
c      CALL STRAIN1(U,SR,J)
C     SECOND TERM OF M_IJ AND |S_IJ|
!$omp parallel do private(DELTA2,SUMSR,ABSSR,DELABS)     
      DO 10 K=1,N3M
      DO 10 I=1,N1M
      DELTA2=DXYZ**(2./3.)
      SUMSR=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >     +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >     +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSSR=SQRT(2.*SUMSR)
      SGSVIS(I,J,K)=ABSSR
      DELABS=DELTA2*ABSSR
      COM(I,K,1)=DELABS*SR(I,K,1)
      COM(I,K,2)=DELABS*SR(I,K,2)        
      COM(I,K,3)=DELABS*SR(I,K,3)
      COM(I,K,4)=DELABS*SR(I,K,4)
      COM(I,K,5)=DELABS*SR(I,K,5)
      COM(I,K,6)=DELABS*SR(I,K,6)
   10 CONTINUE
      CALL FILTER(COM,J,6)

C     FILTERED S_IJ
      CALL FILTER(SR,J,6)
                  
C     ABSOLUTE VALUES OF FILTERED S_IJ
C     M_IJ
!$omp parallel do private(SUMFSR,ABSFSR,FDXYZ,FDABS)   
      DO 20 K=1,N3M
      DO 20 I=1,N1M  
      SUMFSR=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >      +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >      +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSFSR=SQRT(2.*SUMFSR)
      FDXYZ=4.0*DXYZ**(2./3.)
      FDABS=FDXYZ*ABSFSR
      SR(I,K,1)=FDABS*SR(I,K,1)-COM(I,K,1)
      SR(I,K,2)=FDABS*SR(I,K,2)-COM(I,K,2)
      SR(I,K,3)=FDABS*SR(I,K,3)-COM(I,K,3)
      SR(I,K,4)=FDABS*SR(I,K,4)-COM(I,K,4)
      SR(I,K,5)=FDABS*SR(I,K,5)-COM(I,K,5)
      SR(I,K,6)=FDABS*SR(I,K,6)-COM(I,K,6)     
   20 CONTINUE

C     L_IJ
      JP=J+1
!$omp parallel do private(KP,IP)     
      DO 30 K=1,N3M
      KP=KPA(K)
      DO 30 I=1,N1M
      IP=IPA(I)
      UC(I,K,1)=(U(I,J,K,1)+U(IP,J ,K ,1))*0.5
      UC(I,K,2)=(U(I,J,K,2)+U(I ,JP,K ,2))*0.5
      UC(I,K,3)=(U(I,J,K,3)+U(I ,J ,KP,3))*0.5
      COM(I,K,1)=UC(I,K,1)*UC(I,K,1)   
      COM(I,K,2)=UC(I,K,1)*UC(I,K,2)   
      COM(I,K,3)=UC(I,K,1)*UC(I,K,3)   
      COM(I,K,4)=UC(I,K,2)*UC(I,K,2)   
      COM(I,K,5)=UC(I,K,2)*UC(I,K,3)   
      COM(I,K,6)=UC(I,K,3)*UC(I,K,3)   
   30 CONTINUE
      CALL FILTER(UC,J,3)
      CALL FILTER(COM,J,6)
                              
C     RESOLVED TURBULENT STRESS, L_IJ
C     OBTAIN DYNAMIC SMAGORINSKY CONSTANT, C AND SUBGRID-SCALE VISCOSITY.
C     C CAN BE A VERY BIG NEGATIVE VALUE, WHICH MEANS THE PRESENT VERSION
C     OF SGS MODEL WOULD HAVE A PROBLEM FOR 3D GEOMETRIES.
C     CLIPPING OPERATION IS ADDED IN 45-DO LOOP.
      SUMNUM=0.
      SUMDEN=0.
!$omp parallel do private(XL11,XL12,XL13,XL22,XL23,XL33,SQMM,SQLM)
!$omp& reduction(+:SUMDEN,SUMNUM)      
      DO 40 K=1,N3M
      DO 40 I=1,N1M
      XL11=COM(I,K,1)-UC(I,K,1)*UC(I,K,1)
      XL12=COM(I,K,2)-UC(I,K,1)*UC(I,K,2)
      XL13=COM(I,K,3)-UC(I,K,1)*UC(I,K,3)
      XL22=COM(I,K,4)-UC(I,K,2)*UC(I,K,2)
      XL23=COM(I,K,5)-UC(I,K,2)*UC(I,K,3)
      XL33=COM(I,K,6)-UC(I,K,3)*UC(I,K,3)
      SQMM=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >    +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >    +2.*SR(I,K,5)**2+SR(I,K,6)**2
      SQLM=2.*XL12*SR(I,K,2)+XL11*SR(I,K,1)
     >    +2.*XL13*SR(I,K,3)+XL22*SR(I,K,4)
     >    +2.*XL23*SR(I,K,5)+XL33*SR(I,K,6)
      SUMDEN=SUMDEN+SQMM
      SUMNUM=SUMNUM+SQLM
   40 CONTINUE
      SUMNUM=SUMNUM*NXZ
      SUMDEN=SUMDEN*NXZ
      SMAGC(J)=-0.5*SUMNUM/SUMDEN ! C=-0.5*(L_IJ*M_IJ)/(M_IJ*M_IJ)
!      SMAGC(J)=0.065**2          ! Constant Smagorinsky constant

!      CALL DAMPINGFUNC(CS_DAMPING,J) ! Constant Smagorinsky constant
      CS_DAMPING=1.0                ! Constant Smagorinsky constant
      
!$omp parallel do private(TURBV)       
      DO 45 K=1,N3M
      DO 45 I=1,N1M
      TURBV=SMAGC(J)*DXYZ**(2./3.)*SGSVIS(I,J,K)        
      SGSVIS(I,J,K)=AMAX1(TURBV,-RENINV) ! Clipping
C      SGSVIS(I,J,K)=TURBV                 ! No clipping
      SGSVIS(I,J,K)=SGSVIS(I,J,K)*CS_DAMPING !Constant Smagorinsky constant
   45 CONTINUE
      

!      IF (J.LT.N2M) GO TO 1000
1000  CONTINUE

       
      RETURN
      END
      
C************************** SGSCS *****************************
C     THIS SUBROUTINE IS TO OBTAIN SUBGRID-SCALE VISCOSITY. 
C     VALUES ARE EVALUATED AT THE CELL CENTER.
C     #(I,K,1) = #_{11}
C     #(I,K,2) = #_{12}
C     #(I,K,3) = #_{13}
C     #(I,K,4) = #_{22}
C     #(I,K,5) = #_{23}
C     #(I,K,6) = #_{33}
 
C     SMAGORINSKY MODEL
      SUBROUTINE SGSCS(U,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/SGSMEAN/STRM(M2,6),SGSNUT(M2),SGSDS(M2),SMAGC(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)  
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/RANSVISO/RANSVIS,RANSVIS2

      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL COM(M1,M3,6),SR(M1,M3,6),UC(M1,M3,3)
      REAL NXZ,CS_DAMPING,CS
      REAL SMAGC1(M2),SGS1(M2)
      
      REAL XX(5),YY(5),XA,XB
      
      NXZ=1.0/DBLE(N1M*N3M)
      RENINV=1./DBLE(RE)
      
      CS=0.0
      UTAU=SQRT(SQRT(WSM(1)**2+WSM(3)**2))
      

      DO 1000 J=1,N2M
      
      DXYZ=1.0/DX1*DY(J)/DX3
 
      CALL STRAIN(U,SR,J) 
c      CALL STRAIN1(U,SR,J)
C     |S_IJ|
!$omp parallel do private(SUMSR,ABSSR)       
      DO K=1,N3M
      DO I=1,N1M
      SUMSR=2.*SR(I,K,2)**2.+SR(I,K,1)**2.
     >     +2.*SR(I,K,3)**2.+SR(I,K,4)**2.
     >     +2.*SR(I,K,5)**2.+SR(I,K,6)**2.
      ABSSR=SQRT(2.*SUMSR)
      SGSVIS(I,J,K)=ABSSR
      ENDDO
      ENDDO

c      CS=0.3  ! why so large?
      CS=0.1      
      A_plus=26.0
      Y_plus=(Y(J)+DY(J)*0.5)*UTAU*RE
      CS    =CS*(1-Exp(-Y_plus/A_plus))
      SMAGC(J)=Cs**2.0
      
!$omp parallel do private(TURBV)       
      DO 45 K=1,N3M
      DO 45 I=1,N1M
      TURBV=SMAGC(J)*DXYZ**(2./3.)*SGSVIS(I,J,K)        
      SGSVIS(I,J,K)=AMAX1(TURBV,-RENINV)      ! Clipping
C      SGSVIS(I,J,K)=TURBV                    ! No clipping
c      SGSVIS(I,J,K)=SGSVIS(I,J,K)*CS_DAMPING  ! what is CS_DAMPING
   45 CONTINUE

1000  CONTINUE
     
       
      RETURN
      END  
      
      
C************************** SGSSD *****************************
C     THIS SUBROUTINE IS TO OBTAIN SUBGRID-SCALE VISCOSITY. 
C     VALUES ARE EVALUATED AT THE CELL CENTER.
C     #(I,K,1) = #_{11}
C     #(I,K,2) = #_{12}
C     #(I,K,3) = #_{13}
C     #(I,K,4) = #_{22}
C     #(I,K,5) = #_{23}
C     #(I,K,6) = #_{33}

C     SCALE-DEPENDENT DYNAMIC SMAGORINSKY MODEL
      SUBROUTINE SGSSD(U,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)

      COMMON/SGSMEAN/STRM(M2,6),SGSNUT(M2),SGSDS(M2),SMAGC(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)  
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)

      COMMON/RANSVISO/RANSVIS,RANSVIS2

      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL COM(M1,M3,6),SR(M1,M3,6)
      REAL COM1(M1,M3,6),SR1(M1,M3,6),UC(M1,M3,3)
      REAL NXZ,CS_DAMPING
      REAL A1,B1,C1,D1,E1
      REAL A2,B2,C2,D2,E2
      REAL PA(5),BETA,BETAZ(M2)

      NXZ=1.0/DBLE(N1M*N3M)
      RENINV=1./DBLE(RE)

!-------for parallel---------------------------------------------

    
      DO 1000 J=1,N2M
      
      DXYZ=1.0/DX1*DY(J)/DX3
 
!---------A1,B1,C1,D1,E1---------------------------------------------      

      CALL STRAIN(U,SR,J) 
c      CALL STRAIN1(U,SR,J)
C     SECOND TERM OF M_IJ AND |S_IJ|
!$omp parallel do private(DELTA2,SUMSR,ABSSR,DELABS) 
      DO 10 K=1,N3M
      DO 10 I=1,N1M
      DELTA2=DXYZ**(2./3.)
      SUMSR=2.*SR(I,K,2)**2+SR(I,K,1)**2  !total 9
     >     +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >     +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSSR=SQRT(2.*SUMSR)
      DELABS=DELTA2*ABSSR       !*0.5 test scale-dependent cs
      COM1(I,K,1)=DELABS*SR(I,K,1)
      COM1(I,K,2)=DELABS*SR(I,K,2)        
      COM1(I,K,3)=DELABS*SR(I,K,3)
      COM1(I,K,4)=DELABS*SR(I,K,4)
      COM1(I,K,5)=DELABS*SR(I,K,5)
      COM1(I,K,6)=DELABS*SR(I,K,6)
   10 CONTINUE
      
      CALL FILTER(COM1,J,6) !6 indicates number of variables

C     FILTERED S_IJ
      CALL FILTER(SR,J,6)
                  
C     ABSOLUTE VALUES OF FILTERED S_IJ
!$omp parallel do private(SUMFSR,ABSFSR,FDXYZ,FDABS) 
      DO 20 K=1,N3M
      DO 20 I=1,N1M  
      SUMFSR=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >      +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >      +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSFSR=SQRT(2.*SUMFSR)
      FDXYZ=DXYZ**(2./3.)
      FDABS=FDXYZ*ABSFSR
      SR1(I,K,1)=FDABS*SR(I,K,1)
      SR1(I,K,2)=FDABS*SR(I,K,2)
      SR1(I,K,3)=FDABS*SR(I,K,3)
      SR1(I,K,4)=FDABS*SR(I,K,4)
      SR1(I,K,5)=FDABS*SR(I,K,5)
      SR1(I,K,6)=FDABS*SR(I,K,6)  
   20 CONTINUE

C     U_I
C     U_I*U_J
      JP=J+1
!$omp parallel do private(KP,IP) 
      DO 30 K=1,N3M
      KP=KPA(K)
      DO 30 I=1,N1M
      IP=IPA(I)
      
      UC(I,K,1)=(U(I,J,K,1)+U(IP,J ,K ,1))*0.5
      UC(I,K,2)=(U(I,J,K,2)+U(I ,JP,K ,2))*0.5
      UC(I,K,3)=(U(I,J,K,3)+U(I ,J ,KP,3))*0.5
      COM(I,K,1)=UC(I,K,1)*UC(I,K,1)   
      COM(I,K,2)=UC(I,K,1)*UC(I,K,2)   
      COM(I,K,3)=UC(I,K,1)*UC(I,K,3)   
      COM(I,K,4)=UC(I,K,2)*UC(I,K,2)   
      COM(I,K,5)=UC(I,K,2)*UC(I,K,3)   
      COM(I,K,6)=UC(I,K,3)*UC(I,K,3)   
   30 CONTINUE
      
      CALL FILTER(UC,J,3)
      CALL FILTER(COM,J,6)
                                    
C     RESOLVED TURBULENT STRESS, L_IJ
C     OBTAIN DYNAMIC SMAGORINSKY CONSTANT, C AND SUBGRID-SCALE VISCOSITY.
C     C CAN BE A VERY BIG NEGATIVE VALUE, WHICH MEANS THE PRESENT VERSION
C     OF SGS MODEL WOULD HAVE A PROBLEM FOR 3D GEOMETRIES.
C     CLIPPING OPERATION IS ADDED IN 45-DO LOOP.

      A1=0.0
      B1=0.0
!$omp parallel do private(XL11,XL12,XL13,XL22,XL23,XL33) 
!$omp& reduction(+:A1,B1)      
      DO 40 K=1,N3M
      DO 40 I=1,N1M

      XL11=COM(I,K,1)-UC(I,K,1)*UC(I,K,1)
      XL12=COM(I,K,2)-UC(I,K,1)*UC(I,K,2)
      XL13=COM(I,K,3)-UC(I,K,1)*UC(I,K,3)
      XL22=COM(I,K,4)-UC(I,K,2)*UC(I,K,2)
      XL23=COM(I,K,5)-UC(I,K,2)*UC(I,K,3)
      XL33=COM(I,K,6)-UC(I,K,3)*UC(I,K,3)
      
      A1=A1+2.*XL12*SR1(I,K,2)+XL11*SR1(I,K,1)
     >     +2.*XL13*SR1(I,K,3)+XL22*SR1(I,K,4)
     >     +2.*XL23*SR1(I,K,5)+XL33*SR1(I,K,6)

      B1=B1+2.*XL12*COM1(I,K,2)+XL11*COM1(I,K,1)
     >     +2.*XL13*COM1(I,K,3)+XL22*COM1(I,K,4)
     >     +2.*XL23*COM1(I,K,5)+XL33*COM1(I,K,6)

   40 CONTINUE
      
      A1=-8.0*A1*NXZ
      B1=-2.0*B1*NXZ

      C1=0.0
      D1=0.0
      E1=0.0
!$omp parallel do reduction(+:C1,D1,E1) 
      DO 50 K=1,N3M
      DO 50 I=1,N1M

      C1=C1+2.*COM1(I,K,2)*COM1(I,K,2)+COM1(I,K,1)*COM1(I,K,1)
     >     +2.*COM1(I,K,3)*COM1(I,K,3)+COM1(I,K,4)*COM1(I,K,4)
     >     +2.*COM1(I,K,5)*COM1(I,K,5)+COM1(I,K,6)*COM1(I,K,6)
     
      D1=D1+2.*SR1(I,K,2)*SR1(I,K,2)+SR1(I,K,1)*SR1(I,K,1)
     >     +2.*SR1(I,K,3)*SR1(I,K,3)+SR1(I,K,4)*SR1(I,K,4)
     >     +2.*SR1(I,K,5)*SR1(I,K,5)+SR1(I,K,6)*SR1(I,K,6)
      
      E1=E1+2.*SR1(I,K,2)*COM1(I,K,2)+SR1(I,K,1)*COM1(I,K,1)
     >     +2.*SR1(I,K,3)*COM1(I,K,3)+SR1(I,K,4)*COM1(I,K,4)
     >     +2.*SR1(I,K,5)*COM1(I,K,5)+SR1(I,K,6)*COM1(I,K,6)
      
   50 CONTINUE
   
      C1=4.0*C1*NXZ
      D1=64.0*D1*NXZ
      E1=32.0*E1*NXZ
      
!---------A1,B1,C1,D1,E1---------------------------------------------      



!---------A2,B2,C2,D2,E2---------------------------------------------      

      CALL STRAIN(U,SR,J) 
c      CALL STRAIN1(U,SR,J)
C     SECOND TERM OF M_IJ AND |S_IJ|
!$omp parallel do private(DELTA2,SUMSR,ABSSR,DELABS)
      DO 60 K=1,N3M
      DO 60 I=1,N1M
      DELTA2=DXYZ**(2./3.)
      SUMSR=2.*SR(I,K,2)**2+SR(I,K,1)**2  !total 9
     >     +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >     +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSSR=SQRT(2.*SUMSR)
!      SGSVIS(I,J,K)=ABSSR
      DELABS=DELTA2*ABSSR       !*0.5 test scale-dependent cs
      COM1(I,K,1)=DELABS*SR(I,K,1)
      COM1(I,K,2)=DELABS*SR(I,K,2)        
      COM1(I,K,3)=DELABS*SR(I,K,3)
      COM1(I,K,4)=DELABS*SR(I,K,4)
      COM1(I,K,5)=DELABS*SR(I,K,5)
      COM1(I,K,6)=DELABS*SR(I,K,6)
   60 CONTINUE
      CALL FILTER4(COM1,J,6) !6 indicates number of variables

C     FILTERED S_IJ
      CALL FILTER4(SR,J,6)
                  
C     ABSOLUTE VALUES OF FILTERED S_IJ
!$omp parallel do private(SUMFSR,ABSFSR,FDXYZ,FDABS) 
      DO 70 K=1,N3M
      DO 70 I=1,N1M  
      SUMFSR=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >      +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >      +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSFSR=SQRT(2.*SUMFSR)
      FDXYZ=DXYZ**(2./3.)
      FDABS=FDXYZ*ABSFSR
      SR1(I,K,1)=FDABS*SR(I,K,1)
      SR1(I,K,2)=FDABS*SR(I,K,2)
      SR1(I,K,3)=FDABS*SR(I,K,3)
      SR1(I,K,4)=FDABS*SR(I,K,4)
      SR1(I,K,5)=FDABS*SR(I,K,5)
      SR1(I,K,6)=FDABS*SR(I,K,6)  
   70 CONTINUE

C     U_I
C     U_I*U_J
      JP=J+1
!$omp parallel do private(KP,IP)        
      DO 80 K=1,N3M
      KP=KPA(K)
      
      DO 80 I=1,N1M
      IP=IPA(I)
      
      UC(I,K,1)=(U(I,J,K,1)+U(IP,J ,K ,1))*0.5
      UC(I,K,2)=(U(I,J,K,2)+U(I ,JP,K ,2))*0.5
      UC(I,K,3)=(U(I,J,K,3)+U(I ,J ,KP,3))*0.5
      COM(I,K,1)=UC(I,K,1)*UC(I,K,1)   
      COM(I,K,2)=UC(I,K,1)*UC(I,K,2)   
      COM(I,K,3)=UC(I,K,1)*UC(I,K,3)   
      COM(I,K,4)=UC(I,K,2)*UC(I,K,2)   
      COM(I,K,5)=UC(I,K,2)*UC(I,K,3)   
      COM(I,K,6)=UC(I,K,3)*UC(I,K,3)   
   80 CONTINUE

      CALL FILTER4(UC,J,3)
      CALL FILTER4(COM,J,6)
C     RESOLVED TURBULENT STRESS, L_IJ
C     OBTAIN DYNAMIC SMAGORINSKY CONSTANT, C AND SUBGRID-SCALE VISCOSITY.
C     C CAN BE A VERY BIG NEGATIVE VALUE, WHICH MEANS THE PRESENT VERSION
C     OF SGS MODEL WOULD HAVE A PROBLEM FOR 3D GEOMETRIES.
C     CLIPPING OPERATION IS ADDED IN 45-DO LOOP.

      A2=0.0
      B2=0.0
!$omp parallel do private(XQ11,XQ12,XQ13,XQ22,XQ23,XQ33) 
!$omp& reduction(+:A2,B2)     
      DO 90 K=1,N3M
      
      DO 90 I=1,N1M

      XQ11=COM(I,K,1)-UC(I,K,1)*UC(I,K,1)
      XQ12=COM(I,K,2)-UC(I,K,1)*UC(I,K,2)
      XQ13=COM(I,K,3)-UC(I,K,1)*UC(I,K,3)
      XQ22=COM(I,K,4)-UC(I,K,2)*UC(I,K,2)
      XQ23=COM(I,K,5)-UC(I,K,2)*UC(I,K,3)
      XQ33=COM(I,K,6)-UC(I,K,3)*UC(I,K,3)
      
      A2=A2+2.*XQ12*SR1(I,K,2)+XQ11*SR1(I,K,1)
     >     +2.*XQ13*SR1(I,K,3)+XQ22*SR1(I,K,4)
     >     +2.*XQ23*SR1(I,K,5)+XQ33*SR1(I,K,6)

      B2=B2+2.*XQ12*COM1(I,K,2)+XQ11*COM1(I,K,1)
     >     +2.*XQ13*COM1(I,K,3)+XQ22*COM1(I,K,4)
     >     +2.*XQ23*COM1(I,K,5)+XQ33*COM1(I,K,6)

   90 CONTINUE
      
      A2=-32.0*A2*NXZ
      B2=-2.0*B2*NXZ

      C2=0.0
      D2=0.0
      E2=0.0
   
!$omp parallel do reduction(+:C2,D2,E2)       
      DO 11 K=1,N3M
      DO 11 I=1,N1M

      C2=C2+2.*COM1(I,K,2)*COM1(I,K,2)+COM1(I,K,1)*COM1(I,K,1)
     >     +2.*COM1(I,K,3)*COM1(I,K,3)+COM1(I,K,4)*COM1(I,K,4)
     >     +2.*COM1(I,K,5)*COM1(I,K,5)+COM1(I,K,6)*COM1(I,K,6)
     
      D2=D2+2.*SR1(I,K,2)*SR1(I,K,2)+SR1(I,K,1)*SR1(I,K,1)
     >     +2.*SR1(I,K,3)*SR1(I,K,3)+SR1(I,K,4)*SR1(I,K,4)
     >     +2.*SR1(I,K,5)*SR1(I,K,5)+SR1(I,K,6)*SR1(I,K,6)
      
      E2=E2+2.*SR1(I,K,2)*COM1(I,K,2)+SR1(I,K,1)*COM1(I,K,1)
     >     +2.*SR1(I,K,3)*COM1(I,K,3)+SR1(I,K,4)*COM1(I,K,4)
     >     +2.*SR1(I,K,5)*COM1(I,K,5)+SR1(I,K,6)*COM1(I,K,6)
      
   11 CONTINUE
   
      C2=4.0*C2*NXZ
      D2=1024.0*D2*NXZ
      E2=128.0*E2*NXZ
      
!---------A2,B2,C2,D2,E2--------------------------------------------- 

      CALL COEFFICIENT(A1,B1,C1,D1,E1,A2,B2,C2,D2,E2,PA)     
      
      CALL DQRRT(PA,5,BETA,1.0E-5) !DSRRT(PA,BETA,6,5)
      
      BETAZ(J)=AMAX1(BETA,0.125)
      

      
!---------CALCULATING THE SCALE-DEPENDENT SGSVIS --------------------
     
      CALL STRAIN(U,SR,J) 
c      CALL STRAIN1(U,SR,J)
C     SECOND TERM OF M_IJ AND |S_IJ|
!$omp parallel do private(DELTA2,SUMSR,ABSSR,DELABS) 
      DO 21 K=1,N3M
      DO 21 I=1,N1M
      DELTA2=DXYZ**(2./3.)
      SUMSR=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >     +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >     +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSSR=SQRT(2.*SUMSR)
      SGSVIS(I,J,K)=ABSSR
      DELABS=DELTA2*ABSSR
      COM(I,K,1)=DELABS*SR(I,K,1)
      COM(I,K,2)=DELABS*SR(I,K,2)        
      COM(I,K,3)=DELABS*SR(I,K,3)
      COM(I,K,4)=DELABS*SR(I,K,4)
      COM(I,K,5)=DELABS*SR(I,K,5)
      COM(I,K,6)=DELABS*SR(I,K,6)
   21 CONTINUE
      
      CALL FILTER(COM,J,6)

C     FILTERED S_IJ
      CALL FILTER(SR,J,6)
                  
C     ABSOLUTE VALUES OF FILTERED S_IJ
C     M_IJ
!$omp parallel do private(SUMFSR,ABSFSR,FDXYZ,FDABS) 
      DO 31 K=1,N3M
      DO 31 I=1,N1M  
      SUMFSR=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >      +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >      +2.*SR(I,K,5)**2+SR(I,K,6)**2
      ABSFSR=SQRT(2.*SUMFSR)
      FDXYZ=4.0*DXYZ**(2./3.)*BETA !!SCALE-DEPENDENT
      FDABS=FDXYZ*ABSFSR
      SR(I,K,1)=FDABS*SR(I,K,1)-COM(I,K,1)
      SR(I,K,2)=FDABS*SR(I,K,2)-COM(I,K,2)
      SR(I,K,3)=FDABS*SR(I,K,3)-COM(I,K,3)
      SR(I,K,4)=FDABS*SR(I,K,4)-COM(I,K,4)
      SR(I,K,5)=FDABS*SR(I,K,5)-COM(I,K,5)
      SR(I,K,6)=FDABS*SR(I,K,6)-COM(I,K,6)     
   31 CONTINUE

C     U_I
C     U_I*U_J
      JP=J+1
!$omp parallel do private(KP,IP)      
      DO 41 K=1,N3M
      KP=KPA(K)
      
      DO 41 I=1,N1M
      IP=IPA(I)
      
      UC(I,K,1)=(U(I,J,K,1)+U(IP,J ,K ,1))*0.5
      UC(I,K,2)=(U(I,J,K,2)+U(I ,JP,K ,2))*0.5
      UC(I,K,3)=(U(I,J,K,3)+U(I ,J ,KP,3))*0.5
      COM(I,K,1)=UC(I,K,1)*UC(I,K,1)   
      COM(I,K,2)=UC(I,K,1)*UC(I,K,2)   
      COM(I,K,3)=UC(I,K,1)*UC(I,K,3)   
      COM(I,K,4)=UC(I,K,2)*UC(I,K,2)   
      COM(I,K,5)=UC(I,K,2)*UC(I,K,3)   
      COM(I,K,6)=UC(I,K,3)*UC(I,K,3)   
   41 CONTINUE
      
      CALL FILTER(UC,J,3)
      CALL FILTER(COM,J,6)
                              
C     RESOLVED TURBULENT STRESS, L_IJ
C     OBTAIN DYNAMIC SMAGORINSKY CONSTANT, C AND SUBGRID-SCALE VISCOSITY.
C     C CAN BE A VERY BIG NEGATIVE VALUE, WHICH MEANS THE PRESENT VERSION
C     OF SGS MODEL WOULD HAVE A PROBLEM FOR 3D GEOMETRIES.
C     CLIPPING OPERATION IS ADDED IN 45-DO LOOP.
      SUMNUM=0.
      SUMDEN=0.
!$omp parallel do private(XL11,XL12,XL13,XL22,XL23,XL33,SQMM,SQLM)
!$omp& reduction(+:SUMDEN,SUMNUM)      
      DO 51 K=1,N3M
      DO 51 I=1,N1M
      XL11=COM(I,K,1)-UC(I,K,1)*UC(I,K,1)
      XL12=COM(I,K,2)-UC(I,K,1)*UC(I,K,2)
      XL13=COM(I,K,3)-UC(I,K,1)*UC(I,K,3)
      XL22=COM(I,K,4)-UC(I,K,2)*UC(I,K,2)
      XL23=COM(I,K,5)-UC(I,K,2)*UC(I,K,3)
      XL33=COM(I,K,6)-UC(I,K,3)*UC(I,K,3)
      SQMM=2.*SR(I,K,2)**2+SR(I,K,1)**2
     >    +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >    +2.*SR(I,K,5)**2+SR(I,K,6)**2
      SQLM=2.*XL12*SR(I,K,2)+XL11*SR(I,K,1)
     >    +2.*XL13*SR(I,K,3)+XL22*SR(I,K,4)
     >    +2.*XL23*SR(I,K,5)+XL33*SR(I,K,6)
      SUMDEN=SUMDEN+SQMM
      SUMNUM=SUMNUM+SQLM
   51 CONTINUE
      
      SUMNUM=SUMNUM*NXZ
      SUMDEN=SUMDEN*NXZ
      SMAGC(J)=-0.5*SUMNUM/SUMDEN ! C=-0.5*(L_IJ*M_IJ)/(M_IJ*M_IJ)
!      SMAGC(J)=0.1**2          ! Constant Smagorinsky constant

!!      CALL DAMPINGFUNC_PORTE(0.1,DELTA2,CS,J) 
!      SMAGC(J)=CS**2          ! Constant Smagorinsky constant
!      print*,j,SMAGC(J)
!      pause

!      print*,beta,SMAGC(J)
!      pause
            

!      CALL DAMPINGFUNC(CS_DAMPING,J) ! Constant Smagorinsky constant
      CS_DAMPING=1.0                ! Constant Smagorinsky constant
!$omp parallel do private(TURBV) 
      DO 45 K=1,N3M
      DO 45 I=1,N1M
      TURBV=SMAGC(J)*DXYZ**(2./3.)*SGSVIS(I,J,K)        
      SGSVIS(I,J,K)=AMAX1(TURBV,-RENINV)*CS_DAMPING  ! Clipping
C      SGSVIS(I,J,K)=TURBV                 ! No clipping
   45 CONTINUE


!      IF (J.LT.N2M) GO TO 1000
1000  CONTINUE

      IF(MOD(NTIME,100).EQ.0)THEN
      open(1000,file='BETA.plt')
      DO J=1,N2M
      WRITE(1000,*)Y(J),BETAZ(J)
      ENDDO
      CLOSE(1000)
      ENDIF
      
      RETURN
      END
                    
      
      
      SUBROUTINE WALL_MODEL(U,P,PRESG,PRESG3,SGSVIS,NTIME)    
      INCLUDE 'COMMON1.FI'
     
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)

      COMMON/SGSMEAN/STRM(M2,6),SGSNUT(M2),SGSDS(M2),SMAGC(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)

      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/ROUGHNESS/Y0

      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL P(M1,M2,M3)
      REAL COM(M1,M3,6),SR(M1,M3,6)
      REAL COM1(M1,M3,6),SR1(M1,M3,6)
      REAL NXZ,CS_DAMPING
      REAL SGSH(M2)
      REAL PRESG,PRESG3
      INTEGER NTIME
      
      IF(WMODEL.EQ.0.0) THEN
          IFWALL=0.0
          CALL WALLSS(U,P,TIME)
      ELSEIF(WMODEL.EQ.1.0)THEN
c        CALL WALLMODEL_SCHUMANN(U,TAUW,TAUB,SGSVIS) 
          IF(Y0.EQ.0.0) THEN
              CALL WALLMODEL_INSTANT(U,TAUW,TAUB,SGSVIS)
          ELSE
              CALL WALLMODEL_ROUGH(U,TAUW,TAUB,SGSVIS)
          ENDIF
          IFWALL=1.0
      ELSEIF(WMODEL.EQ.2.0)THEN
c        CALL WALLMODEL_GROTZBACH(U,TAUW,TAUB,SGSVIS) 
          IF(Y0.EQ.0.0) THEN
              CALL WALLMODEL_IWMLES_SMOOTH(U,P,PRESG,PRESG3,TAUW,TAUB,
     >                                    SGSVIS,NTIME)
          ELSE
              CALL WALLMODEL_IWMLES_ROUGH(U,P,PRESG,PRESG3,TAUW,TAUB,
     >                                    SGSVIS,NTIME)
      ENDIF
          IFWALL=1.0      
      ELSEIF(WMODEL.EQ.3.0)THEN      
          CALL WALLMODEL_MARUSIC1(U,TAUW,TAUB,SGSVIS) 
          IFWALL=1.0    
      ELSEIF(WMODEL.EQ.4.0)THEN      
          CALL  WALLMODEL_WERNER(U,TAUW,TAUB,SGSVIS,NTIME) 
          IFWALL=1.0      
      ELSEIF(WMODEL.EQ.5.0)THEN
          CALL WALLMODEL_PIOMELLI(U,TAUW,TAUB,SGSVIS) 
          IFWALL=1.0          
      ELSE      
          PRINT*,'WRONG! NO CORRESPONDING WALL MODEL'      
          IFWALL=0.0
      ENDIF
      
      IF(MOD(NTIME,100).EQ.0.0)THEN
          OPEN(1111,file='SGSMAG.plt')
          DO J=1,N2M
          WRITE(1111,*)J,Y(J),SMAGC(J)
          ENDDO
          CLOSE(1111)
      ENDIF

 
      RETURN
      END
      
      
       SUBROUTINE SOURCE_PAR(U,RUH1,RUH2,RUH3,NTIME)  
       INCLUDE 'COMMON1.FI'   
  
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
   
      REAL RUH1(0:M1,0:M2,0:M3)   
      REAL RUH2(0:M1,0:M2,0:M3)       
      REAL RUH3(0:M1,0:M2,0:M3)
      REAL U(0:M1,0:M2,0:M3,3)

	INTEGER IPP,JPP,KPP
	INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UP,VP,WP
    	REAL    UF,VF,WF,G
	REAL    QU,QV,QW

	QQ=1000.0
	G=9.8/(RE*1.5/100000/HHH)**2
       DO II=1,N_PAR
	  IIP=II 
          XP=PAR(1,IIP)
	  YP=PAR(2,IIP)
	  ZP=PAR(3,IIP)
	  UP=PAR(4,IIP)
	  VP=PAR(5,IIP)
	  WP=PAR(6,IIP) 
       	  CALL INDEX_IJK (XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
 	  CALL FORCE_PAR (UP,VP,WP,UF,VF,WF,QU,QV,QW)
        IPP=IP
        JPP=JP
        KPP=KP
       IF(XP.GT.XC(IP))IPP=IPA(IP)
       IF(YP.GT.YC(JP))JPP=JPA(JP)
       IF(ZP.GT.ZC(KP))KPP=KPA(KP)       
        RUH1(IPP,JP,KP) =RUH1(IPP,JP,KP)-Mass_P*QQ*QU*DX1*DX3/H(JP)
        RUH2(IP,JPP,KP) =RUH2(IP,JPP,KP)-Mass_P*QQ*(QV+G)*DX1*DX3/H(JPP)
        RUH3(IP,JP,KPP) =RUH3(IP,JP,KPP)-Mass_P*QQ*QW*DX1*DX3/H(JP)

       ENDDO
       
       RETURN
       END
     
C--------------------COEFFICIENT OF POLYNOMIAL----------------------------------
C     THIS SUBROUTINE IS TO CALCULATE THE SCALE-DEPENDENT BETA.
C     PA
      SUBROUTINE COEFFICIENT(A1,B1,C1,D1,E1,A2,B2,C2,D2,E2,PA)
      REAL PA(5)
      
!      PA(6)=(B2*C1-B1*C2)/(A1*D2)
!      PA(5)=(A1*C2-B2*E1)/(A1*D2)
!      PA(4)=(B2*D1+B1*E2-A2*C1)/(A1*D2)
!      PA(3)=(A2*E1-A1*E2)/(A1*D2)
!      PA(2)=-(A2*D1-B1*D2)/(A1*D2)
!      PA(1)=1.0

      PA(5)=(B2*C1-B1*C2)/(A1*D2)
      PA(4)=(A1*C2-B2*E1)/(A1*D2)
      PA(3)=(B2*D1+B1*E2-A2*C1)/(A1*D2)
      PA(2)=(A2*E1-A1*E2)/(A1*D2)
      PA(1)=-(A2*D1+B1*D2)/(A1*D2)
!      print*,pa
!      pause

      END   
C********************** DAMPINGFUNC ***************************
C     THIS SUBROUTINE IS TO NEAR WALL DAMPING FUNCTION.
C     SMAGORINSKY CS MODIFICATION WITH CS*[1-EXP(-YPLUS/APLUS)]
C     ONLY FOR SMAGROINGSKY MODEL CS=0.1, SEE SUBROUTINE SGS

      SUBROUTINE DAMPINGFUNC(CS_DAMPING,J)  
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/MESH2/Y(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS

      COMMON/WSMEM/WSM(3),WSMO(3)  
      COMMON/WALLTAO/TAOW

      REAL CS_DAMPING 
      REAL APLUS

      APLUS=17 !26
      
      JP=J+1
      IF(J.LT.(M2-1)/2)THEN
      YH=0.5*(Y(J)+Y(JP))
      ELSE
      YH=1.0-0.5*(Y(J)+Y(JP))
      ENDIF
      
      IF(WMODEL.EQ.0)THEN
      UTAO=SQRT(SQRT(WSM(1)**2+WSM(3)**2))
      ELSE
      UTAO=SQRT(TAOW)
      ENDIF
            
      YPLUS=YH*UTAO*RE
      CS_DAMPING =1-EXP(-YPLUS/APLUS)      ! Cui's book
!      CS_DAMPING=1-EXP(-(YPLUS/APLUS)**3) !Piomelli and Balaras 2002
!      CS_DAMPING =(1-EXP(-YPLUS/APLUS))**2  ! Larsson and Kawai, 2010
      RETURN
      END
      
C********************** WALLMODEL_SCHUMANN ***************************
C     THIS SUBROUTINE IS TO NEAR WALL MODEL.  
C     EQUILIBRIUM LAWS OF SCHUMANN (1975) Txy,w=<Tw>*u/<u>
C                                         Tyz,w=viscos*w/<w>


      SUBROUTINE WALLMODEL_SCHUMANN(U,TAUW,TAUB,SGSVIS)  ! Modified by Ruifeng Hu, Limin Wang and Yangyue Zhang
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
c      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      
      REAL U1UP,U1DOWN,W1UP,W1DOWN
      REAL DELTAP,TAUW1,TAUW2,LAMBDA
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      INTEGER I,K 
      REAL KAPA,B,TEMP
      REAL RETUW,U1,U2,REINV
      real TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      
      U1DOWN=0.0
      U1UP=0.0
      W1DOWN=0.0
      W1UP=0.0
      REINV=1.0/DBLE(RE)
      
      DO K=1,N3M
      DO I=1,N1M 
      U1DOWN=U1DOWN+U(I,3,K,1)
      W1DOWN=W1DOWN+U(I,3,K,3)
      U1UP=U1UP+U(I,N2M-2,K,1)  
      W1UP=W1UP+U(I,N2M-2,K,3)
c      U1DOWN=U1DOWN+U(I,1,K,1)
c      W1DOWN=W1DOWN+U(I,1,K,3)
c      U1UP=U1UP+U(I,N2M,K,1)  
c      W1UP=W1UP+U(I,N2M,K,3)
      ENDDO
      ENDDO
      
      U1DOWN=U1DOWN/DBLE(N1M*N3M)
      U1UP=U1UP/DBLE(N1M*N3M)
      W1DOWN=W1DOWN/DBLE(N1M*N3M)
      W1UP=W1UP/DBLE(N1M*N3M)
      
C      U1DOWN=VMO(1,1)
C      U1UP=VMO(N2M,1)
c      U1DOWN=VMO(2,1)
c      U1UP=VMO(N2M-1,1)
 
!SOLVED BY NEWTON ITERATIONS
!!!!!!!!!!!!!!!!!	>>>>>>>>>>>>>>>>>>>
      KAPA=0.41

c     bottom wall
c      DELTAP=DY(1)/2.0 
      DELTAP=DY(1)+DY(2)+DY(3)/2.0  
      
c      U_ave=SQRT(U1DOWN**2.+W1DOWN**2.)   
      U_ave=U1DOWN  
      B=5.2   
      
      !U_tau2=0.2
      !U_tau1=0.000001
      !DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      !TEMP=0.5*U_tau2+0.5*U_tau1
      !FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      !FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B      
      !TEMP=FUtau1*FUtau2
      !IF(TEMP.LE.0.0)THEN
      !U_tau2=0.5*U_tau2+0.5*U_tau1
      !ELSE    
      !U_tau1=0.5*U_tau2+0.5*U_tau1
      !ENDIF      
      !ENDDO
       
      
      U_tau1=1.0
      U_tau2=0.000001
      FUtau3=1.0
      FUtau1=0.1
      LAMBDA=1.0
      
      DO WHILE(abs(U_tau2-U_tau1).GT.1.0E-6)
      U_tau1=U_tau2
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=-U_ave/U_tau1**2-1./KAPA/U_tau1
      U_tau2=U_tau1-LAMBDA*FUtau1/FUtau2
      FUtau3=U_ave/U_tau2-ALOG(U_tau2*DELTAP*RE)/KAPA-B      
      IF (ABS(FUtau3).GT.ABS(FUtau1)) THEN
      LAMBDA=0.5*LAMBDA
      ENDIF      
      ENDDO 
      
     
      TAUW1=U_tau1**2

!      write(*,*)'1',U_tau1,VRMS(1,4),RANSVIS,DUDY_TEMP
      DO J=1,N3M
      DO I=1,N1M 
      TAUW(I,J,1)=TAUW1*U(I,3,J,1)/U_ave
c      TAUB(I,J,1)=1.0/RE*U(I,3,J,3)/DELTAP
      TAUB(I,J,1)=TAUW1*U(I,3,J,3)/U_ave
      ENDDO
      ENDDO

      TAOW=(U_tau1)**2
      WSM(1)=TAOW
  
c      WRITE(*,*) WSM(1)
      
c     top wall
c      DELTAP=DY(N2M)/2.0   
      DELTAP=DY(N2M)+DY(N2M-1)+DY(N2M-2)/2.0 
      
c      U_ave= SQRT(U1UP**2+W1UP**2)   
      U_ave= U1UP
      B=5.2   
      
c      write(*,*) U_ave, DELTAP
      
      !U_tau2=0.2
      !U_tau1=0.000001
      !DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      !TEMP=0.5*U_tau2+0.5*U_tau1
      !FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      !FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B      
      !TEMP=FUtau1*FUtau2
      !IF(TEMP.LE.0.0)THEN
      !U_tau2=0.5*U_tau2+0.5*U_tau1
      !ELSE    
      !U_tau1=0.5*U_tau2+0.5*U_tau1
      !ENDIF   
      !ENDDO
      
      U_tau1=1.0
      U_tau2=0.000001
      FUtau3=1.0
      FUtau1=0.1
      LAMBDA=1.0
      
      DO WHILE(abs(U_tau2-U_tau1).GT.1.0E-6)
      U_tau1=U_tau2
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=-U_ave/U_tau1**2-1./KAPA/U_tau1
      U_tau2=U_tau1-FUtau1/FUtau2   
      FUtau3=U_ave/U_tau2-ALOG(U_tau2*DELTAP*RE)/KAPA-B      
      IF (ABS(FUtau3).GT.ABS(FUtau1)) THEN
      LAMBDA=0.5*LAMBDA
      ENDIF 
      ENDDO 
      
c      write(*,*) U_tau2

      TAUW1=U_tau1**2
      
      DO J=1,N3M
      DO I=1,N1M 
      TAUW(I,J,2)=TAUW1*U(I,N2M-2,J,1)/U_ave
c      TAUB(I,J,2)=1.0/RE*U(I,N2M-2,J,3)/DELTAP
      TAUB(I,J,2)=TAUW1*U(I,N2M-2,J,3)/U_ave
      ENDDO
      ENDDO
      
      RETURN
      END
      
cccc  Instantaneous wall model      
      SUBROUTINE WALLMODEL_INSTANT(U,TAUW,TAUB,SGSVIS)  ! Modified by Ruifeng Hu, Limin Wang and Yangyue Zhang
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/ROUGHNESS/Y0
c      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      
      REAL U1DOWN(M1,M3),W1DOWN(M1,M3)
      REAL DELTAP,TAUW1,TAUW2,LAMBDA
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      INTEGER I,K 
      REAL KAPA,B,TEMP(M1,M3),TEMP1
      REAL RETUW,U1,U2,REINV
      real TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      
      REINV=1.0/DBLE(RE)
      
c      DO K=1,N3M
c      DO I=1,N1M 
c      U1DOWN=U1DOWN+U(I,3,K,1)
c      W1DOWN=W1DOWN+U(I,3,K,3)
c      ENDDO
c      ENDDO
      
c      U1DOWN=U1DOWN/DBLE(N1M*N3M)
c      U1UP=U1UP/DBLE(N1M*N3M)
c      W1DOWN=W1DOWN/DBLE(N1M*N3M)
c      W1UP=W1UP/DBLE(N1M*N3M)
      
c      U1DOWN=VMO(1,1)
c      U1UP=VMO(N2M,1)
c      U1DOWN=VMO(3,1)
c      U1UP=VMO(N2M-2,1)
 
!SOLVED BY NEWTON ITERATIONS
!!!!!!!!!!!!!!!!!	>>>>>>>>>>>>>>>>>>>
      KAPA=0.41
      B=5.2 

c     bottom wall
c      DELTAP=DY(1)/2.0 
      DELTAP=DY(1)+DY(2)+DY(3)/2.0 
      WSM(1)=0.0

CCC   FILTERING VELOCITY AT THE FIRST GRID
c!$OMP PARALLEL DO PRIVATE(IP,IM,FILVM,FILVC,FILVP)
c      DO 10 K=1,N3M
c      DO 10 I=1,N1M
c      IP=IPA(I)
c      IM=IMA(I)
c      FILVM=U(IM,3,K,1)
c      FILVC=U(I ,3,K,1)
c      FILVP=U(IP,3,K,1)      
cC     Simpson's 1/3 rule
c      TEMP(I,K)=(FILVM+4.0*FILVC+FILVP)/6.0
c   10 CONTINUE
c
cC     FILTERING IN X3
c!$OMP PARALLEL DO PRIVATE(KP,KM,FILVM,FILVC,FILVP)
c      DO 20 K=1,N3M
c      KP=KPA(K)
c      KM=KMA(K)
c      DO 20 I=1,N1M
c      FILVM=TEMP(I,KM)
c      FILVC=TEMP(I,K )
c      FILVP=TEMP(I,KP)
cc     Simpson rule
c      U1DOWN(I,K)=(FILVM+4.0*FILVC+FILVP)/6.0
c   20 CONTINUE
c
c!$OMP PARALLEL DO PRIVATE(IP,IM,FILVM,FILVC,FILVP)
c      DO 30 K=1,N3M
c      DO 30 I=1,N1M
c      IP=IPA(I)
c      IM=IMA(I)
c      FILVM=U(IM,3,K,3)
c      FILVC=U(I ,3,K,3)
c      FILVP=U(IP,3,K,3)      
cC     Simpson's 1/3 rule
c      TEMP(I,K)=(FILVM+4.0*FILVC+FILVP)/6.0
c   30 CONTINUE
c
cC     FILTERING IN X3
c!$OMP PARALLEL DO PRIVATE(KP,KM,FILVM,FILVC,FILVP)
c      DO 40 K=1,N3M
c      KP=KPA(K)
c      KM=KMA(K)
c      DO 40 I=1,N1M
c      FILVM=TEMP(I,KM)
c      FILVC=TEMP(I,K )
c      FILVP=TEMP(I,KP)
cc     Simpson rule
c      W1DOWN(I,K)=(FILVM+4.0*FILVC+FILVP)/6.0
c   40 CONTINUE
      

c!$OMP PARALLEL DO PRIVATE(U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
c!$OMP& ,LAMBDA,TAUW1) REDUCTION(+:WSM)      
      DO K=1,N3M
      DO I=1,N1M
          
c      U1DOWN=U(I,3,K,1)
c      W1DOWN=U(I,3,K,3)
      
c      U_ave=SQRT(U1DOWN(I,K)**2.+W1DOWN(I,K)**2.)   
C      U_ave=U1DOWN(I,K)
      U_ave=SQRT(U(I,3,K,1)**2.+U(I,3,K,3)**2.)   
      
     
      U_tau1=1.0
      U_tau2=0.000001
      FUtau3=1.0
      FUtau1=0.1
      LAMBDA=1.0
     
      DO WHILE(abs(U_tau2-U_tau1).GT.1.0E-6)
      U_tau1=U_tau2
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=-U_ave/U_tau1**2-1./KAPA/U_tau1
      U_tau2=U_tau1-LAMBDA*FUtau1/FUtau2
      FUtau3=U_ave/U_tau2-ALOG(U_tau2*DELTAP*RE)/KAPA-B      
      IF (ABS(FUtau3).GT.ABS(FUtau1)) THEN
      LAMBDA=0.5*LAMBDA
      ENDIF      
      ENDDO       

C      U_tau2=1.0
C      U_tau1=0.000001
C      DO WHILE(abs(U_tau2-U_tau1).GT.0.000001)
C      TEMP1=0.5*U_tau2+0.5*U_tau1
C      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
C      FUtau2=U_ave/TEMP1-ALOG(TEMP1*DELTAP*RE)/KAPA-B
C      TEMP1=FUtau1*FUtau2
C      IF(TEMP1.LE.0.0)THEN
C      U_tau2=0.5*U_tau2+0.5*U_tau1
C      ELSE
C      U_tau1=0.5*U_tau2+0.5*U_tau1
C      ENDIF
C      ENDDO
     
      TAUW1=U_tau1**2

C      TAUW1=(KAPA/ALOG(DELTAP/Y0))*(U1DOWN(I,J)**2+W1DOWN(I,J)**2)
      
      WSM(1)=WSM(1)+TAUW1

!      write(*,*)'1',U_tau1,VRMS(1,4),RANSVIS,DUDY_TEMP
 
      TAUW(I,K,1)=TAUW1*U(I,3,K,1)/U_ave
      TAUB(I,K,1)=TAUW1*U(I,3,K,3)/U_ave
      
      ENDDO
      ENDDO

      WSM(1)=WSM(1)/REAL(N1M*N3M)
      TAOW=WSM(1)

      
      RETURN
      END
      
cccc  Rough wall model      
      SUBROUTINE WALLMODEL_ROUGH(U,TAUW,TAUB,SGSVIS)  ! Modified by Ruifeng Hu, Limin Wang and Yangyue Zhang
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/ROUGHNESS/Y0
c      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      
      REAL U1DOWN(M1,M3),W1DOWN(M1,M3)
      REAL DELTAP,TAUW1,TAUW2,LAMBDA
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      INTEGER I,K 
      REAL KAPA,B,TEMP(M1,M3),TEMP1
      REAL RETUW,U1,U2,REINV
      real TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      
      REINV=1.0/DBLE(RE)
 
      KAPA=0.41
      B=5.2 

c     bottom wall
      DELTAP=DY(1)+DY(2)+DY(3)/2.0 
      WSM(1)=0.0

 
      DO K=1,N3M
      DO I=1,N1M
          
      U_ave=SQRT(U(I,3,K,1)**2.+U(I,3,K,3)**2.)       

      TAUW1=(KAPA/ALOG(DELTAP/Y0)*U_ave)**2.0
      
      WSM(1)=WSM(1)+TAUW1

      TAUW(I,K,1)=TAUW1*U(I,3,K,1)/U_ave
      TAUB(I,K,1)=TAUW1*U(I,3,K,3)/U_ave
      
      ENDDO
      ENDDO

      WSM(1)=WSM(1)/REAL(N1M*N3M)
      
      TAOW=WSM(1)
      
c      write(*,*) ntime, taow
      
      RETURN
      END
      
cccc  Rough IWMLES model 
!*******************************************************************************
      SUBROUTINE WALLMODEL_IWMLES_ROUGH(U,P,PRESG,PRESG3,TAUW,TAUB,
     >                                    SGSVIS,NTIME)
      ! ORIGINATED FROM LESGO-JHU, MODIFIED BY RUI-FENG HU
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/ROUGHNESS/Y0
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL U1DOWN(M1,M3),W1DOWN(M1,M3)
      REAL DELTAP,TAUW1,TAUW2,LAMBDA
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL P(M1,M2,M3)
      REAL KAPA,B,TEMP(M1,M3),TEMP1
      REAL RETUW,U1,U2,REINV
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL PRESG,PRESG3
      REAL wm_dt
      INTEGER NTIME
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      
      REINV=1.0/DBLE(RE)
 
      KAPA = 0.4
      B    = 5.0 
      
CCCCC PARAMETERS FOR IWMLES
      iwm_dirx = 1  ! direction x
      iwm_diry = 2  ! direction z
      
      iwm_DN   = 2
      
      iwm_Lu   = 1   ! index for integral of u
      iwm_Luu  = 2   ! index for integral of uu
      iwm_Lv   = 3   ! etc.
      iwm_Lvv  = 4
      iwm_Luv  = 5
      
      iwm_LN   = 5  ! the total number of integrals that need to be calculated
      
      iwm_ntime_skip = 5
      
      IF (NTIME.EQ.1) CALL iwm_init_rough(U)
      
      CALL iwm_wallstress_rough(U,P,PRESG,PRESG3,TAUW,TAUB,NTIME)   

c     bottom wall
      DELTAP=DY(1)+DY(2)+DY(3)/2.0 
      WSM(1)=0.0
      
      DO K=1,N3M
      DO I=1,N1M
          
C      U_ave=SQRT(U(I,3,K,1)**2.+U(I,3,K,3)**2.)
C      TAUW1=(KAPA/ALOG(DELTAP/Y0)*U_ave)**2.0 
      TAUW1=SQRT(TAUW(I,K,1)**2.+TAUW(I,K,3)**2.)
      WSM(1)=WSM(1)+TAUW1
C      TAUW(I,K,1)=TAUW1*U(I,3,K,1)/U_ave
C      TAUB(I,K,1)=TAUW1*U(I,3,K,3)/U_ave
      
      ENDDO
      ENDDO
      
      WSM(1)=WSM(1)/REAL(N1M*N3M)      
      TAOW=WSM(1)
      
c      write(*,*) ntime, taow
      
      RETURN
      END
      
!******************************************************************************* 
      subroutine iwm_wallstress_rough(U,P,PRESG,PRESG3,TAUW,TAUB,NTIME)      
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL PRESG,PRESG3   
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      
C      use sim_param , only : dudz, dvdz, txz, tyz
C      implicit none

      ! Calculate the time step used in the integral wall model
      !! DO NOT USE iwm_ntime_skip=1 !! !! this number is hard coded to prevent any
      !! mis-use...
      if (mod(NTIME,iwm_ntime_skip)==1) then
          wm_dt = DT
      else
          wm_dt = wm_dt+DT
      end if

      ! Compute the wall stress
      if(mod(NTIME,iwm_ntime_skip)==0) then
      ! gather flow status, update the integrated unsteady term, convective term,
      ! turbulent diffusion term etc.
      call iwm_calc_lhs_rough(U,P,PRESG,PRESG3)
      ! the subroutine to calculate wall stress
      call iwm_calc_wallstress_rough
      ! this is to monitor any quantity from the iwm, useful debugging tool
      call iwm_monitor_rough(NTIME)
      end if

      ! Imposing txz, tyz, dudz, dvdz every time step even iwm_* are not computed
      ! every time step.
      do I=1,N1M
      do K=1,N3M
          
      ! wall stress, use the value calculated in iwm, note the negative sign
      TAUW(I,K,1) = wm_tauwx(I,K)
      TAUW(I,K,3) = wm_tauwy(I,K)

      ! Note: Using the value calculated in iwm, note the positive sign
C      dudz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
C      dvdz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_diry)
      end do
      end do

      end subroutine iwm_wallstress_rough
      
!*******************************************************************************
      subroutine iwm_init_rough(U)
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/ROUGHNESS/Y0
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL vonk
      REAL usinit, uinit, vinit, Dzp
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      REAL U(0:M1,0:M2,0:M3,3)
      
      vonk = 0.4
      
      ! at the height of the first grid point.
      Dzp=DY(1)+DY(2)+DY(3)/2.0 
C      Dzp=dz/2._rprec
      
      ! initial value for the x-velocity at first grid point
c      uinit = usinit/vonk*log(Dzp/Y0)
      uinit = 0.0
      DO K=1,N3M
      DO I=1,N1M
      uinit = uinit + U(I,3,K,1)
      ENDDO
      ENDDO
      uinit = uinit/real(n1m*n3m)

      ! initial value for us (the friction velocity)
      usinit= uinit*vonk/log(Dzp/Y0)
      
      ! initial value for the y-velocity at first grid point
      vinit = 0.0 
      
      DO K=1,N3M
      DO I=1,N1M
          
      ! us in x, y directions
      wm_utx(I,K) = usinit
      wm_uty(I,K) = 0.0
      
      ! wall stress in x, y directions
      wm_tauwx(I,K) = usinit**2.0
      wm_tauwy(I,K) = 0.0
      
      ! filitered velocity at the first grid point in x, y directions
      wm_flt_tagvel  (I,K,iwm_dirx) = uinit
      wm_flt_tagvel  (I,K,iwm_diry) = vinit
      wm_flt_tagvel_m(I,K,iwm_dirx) = uinit
      wm_flt_tagvel_m(I,K,iwm_diry) = vinit
      
      ! pressure at first grid point
      wm_flt_p(I,K) = 0.0
      
      ! integrals of Lu, Lv, etc.
      wm_inte(I,K,iwm_Lu)    = uinit*Dzp
      wm_inte(I,K,iwm_Lv)    = 0.0
      wm_inte(I,K,iwm_Luu)   = uinit**2.0*Dzp
      wm_inte(I,K,iwm_Lvv)   = 0.0
      wm_inte(I,K,iwm_Luv)   = 0.0
      wm_inte_m(I,K,iwm_Lu)  = uinit*Dzp
      wm_inte_m(I,K,iwm_Lv)  = 0.0
      wm_inte_m(I,K,iwm_Luu) = uinit**2.0*Dzp
      wm_inte_m(I,K,iwm_Lvv) = 0.0
      wm_inte_m(I,K,iwm_Luv) = 0.0
      
      ! each term in the integral equation and top/bottom derivatives
      wm_unsdy(I,K,iwm_dirx)   = 0.0
      wm_unsdy(I,K,iwm_diry)   = 0.0
      wm_conv(I,K,iwm_dirx)    = 0.0
      wm_conv(I,K,iwm_diry)    = 0.0
      wm_PrsGrad(I,K,iwm_dirx) = 0.0
      wm_PrsGrad(I,K,iwm_diry) = 0.0
      wm_diff(I,K,iwm_dirx)    = 0.0
      wm_diff(I,K,iwm_diry)    = 0.0
      wm_LHS(I,K,iwm_dirx)     = -uinit*Dzp
      wm_LHS(I,K,iwm_diry)     = -uinit*Dzp
      wm_dudzT(I,K,iwm_dirx) = usinit/vonk/Dzp
      wm_dudzT(I,K,iwm_diry) = 0.0
      wm_dudzB(I,K,iwm_dirx) = usinit/vonk/Y0
      wm_dudzB(I,K,iwm_diry) = 0.0
      
      ! filtered friction velocity and the filtering time scale, tR<1
      wm_flt_us(I,K) = usinit
      wm_tR(I,K)     = (CFLMAX*ALX/N1M/uinit)/(Dzp/vonk/usinit)
      
      ! cell height and imposed roughness length
      wm_Dz(I,K) = Dzp
      wm_z0(I,K) = Y0   
      
      ! linear correction to the log profile
      wm_Ax(I,K) = 0.0
      wm_Ay(I,K) = 0.0
      
      ENDDO
      ENDDO       
      
      ! time step seen by the iwm
c      wm_dt=iwm_ntime_skip*CFLMAX*ALX/N1M/uinit  
      wm_dt=iwm_ntime_skip*DT
      
      end subroutine iwm_init_rough

!*******************************************************************************
      subroutine iwm_calc_lhs_rough(U,P,PRESG,PRESG3)
!*******************************************************************************
! Ths subroutine calculates the left hand side of the iwm system.
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/ROUGHNESS/Y0
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL PRESG,PRESG3   
      REAL u_inst(0:M1,0:M3),v_inst(0:M1,0:M3),w_inst(0:M1,0:M3)
      REAL p_inst(0:M1,0:M3)
      REAL uux, uvx, uvy, vvy, ux, vy
      REAL phip, phim   
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      real wm_dt,p_bar

      
      do I=1,N1M
      do K=1,N3M
          ! update the u, v for previous time step
          wm_flt_tagvel_m(I,K,iwm_dirx) =                     
     >         wm_flt_tagvel(I,K,iwm_dirx)
          wm_flt_tagvel_m(I,K,iwm_diry) =                     
     >         wm_flt_tagvel(I,K,iwm_diry)
          
          ! get the instantaneous field
          u_inst(I,K) = U(I,3,K,1)
          v_inst(I,K) = U(I,3,K,3)
          w_inst(I,K) = U(I,2,K,2)*0.25
          p_inst(I,K) = P(I,1,K)
c          write(*,*) i,k,p_inst(I,K)
      end do
      end do
      
      ! obtain the pressure fluctuations
      p_bar = 0.0
      do I = 1, N1M
      do K = 1, N3M
          p_bar = p_bar+p_inst(I,K)
c	  if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c      	  write(*,*) p_bar
c          end if
      end do
      end do
      p_bar=p_bar/REAL(N1M)/REAL(N3M)

      do I = 1, N1M
      do K = 1, N3M
          p_inst(I,K) = p_inst(I,K)-p_bar
c          if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c      	  write(*,*) p_bar,p_inst(i,k)
c          end if
      end do
      end do

      ! all the data enters must be filtered (see Anderson & Meneveau 2011 JFM)
      !call test_filter ( u_inst )
      !call test_filter ( v_inst )
      !call test_filter ( w_inst )
      !call test_filter ( p_inst )

      !temporal filtering
      do I=1,N1M
      do K=1,N3M
          wm_flt_tagvel(I,K,iwm_dirx) =                                
     >     wm_flt_tagvel(I,K,iwm_dirx)*(1.0-wm_tR(I,K))    
     >         + u_inst(I,K)*wm_tR(I,K)
          wm_flt_tagvel(I,K,iwm_diry) =                                 
     >         wm_flt_tagvel(I,K,iwm_diry)*(1.0-wm_tR(I,K))    
     >         + v_inst(I,K)*wm_tR(I,K)
          wm_flt_p(I,K) = wm_flt_p(I,K)                            
     >         * (1.0-wm_tR(I,K))                                       
     >         + p_inst(I,K)*wm_tR(I,K)
c      if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c      write(*,*) p_inst(i,k),                                 
c     >     wm_flt_p(i,k)
c      end if
      end do
      end do

      ! calculate LHS, calculation of the integrals is done from the last time step
      ! in the subroutine iwm_calc_wallstress, so is iwm_diff
      do K = 1, N3M
          KP=KPA(K)
          KM=KMA(K)
      do I = 1, N1M
          IP=IPA(I)
          IM=IMA(I)
          ! the unsteady term
          wm_unsdy(I,K,iwm_dirx) =                                      
     >         (wm_inte(I,K,iwm_Lu)-wm_inte_m(I,K,iwm_Lu))/wm_dt
          wm_unsdy(I,K,iwm_diry) =                                      
     >         (wm_inte(I,K,iwm_Lv)-wm_inte_m(I,K,iwm_Lv))/wm_dt

          ! the convective term
          phip = wm_inte(IP,K,iwm_Luu)
          phim = wm_inte(IM,K,iwm_Luu)
          uux = (phip-phim)*DX1/2.0
          phip = wm_inte(I,KP,iwm_Luv)
          phim = wm_inte(I,KM,iwm_Luv)
          uvy = (phip-phim)*DX3/2.0
          phip = wm_inte(IP,K,iwm_Luv)
          phim = wm_inte(IM,K,iwm_Luv)
          uvx = (phip-phim)*DX1/2.0
          phip = wm_inte(I,KP,iwm_Lvv)
          phim = wm_inte(I,KM,iwm_Lvv)
          vvy = (phip-phim)*DX3/2.0
          phip = wm_inte(IP,K,iwm_Lu )
          phim = wm_inte(IM,K,iwm_Lu )
          ux = (phip-phim)*DX1/2.0
          phip = wm_inte(I,KP,iwm_Lv )
          phim = wm_inte(I,KM,iwm_Lv )
          vy = (phip-phim)*DX3/2.0
          wm_conv(I,K,iwm_dirx) = uux + uvy                             
     >          - wm_flt_tagvel_m(I,K,iwm_dirx)*(ux+vy)
          wm_conv(I,K,iwm_diry) = uvx + vvy                           
     >          - wm_flt_tagvel_m(I,K,iwm_diry)*(ux+vy)

          ! the pressure gradient term
          phip = wm_flt_p(IP,K)
          phim = wm_flt_p(IM,K)
          wm_PrsGrad(I,K,iwm_dirx) = (phip-phim)*DX1/2.0                
     >         * wm_Dz(I,K) + PRESG*wm_Dz(I,K)
          phip = wm_flt_p(I,KP)
          phim = wm_flt_p(I,KM)
          wm_PrsGrad(I,K,iwm_diry) = (phip-phim)*DX3/2.0   
     >         * wm_Dz(I,K) + PRESG3*wm_Dz(I,K)

          ! the left hand side
          ! this is the integrated momentum equation, except for the Lu term
          wm_lhs(I,K,iwm_dirx) = -wm_inte(I,K,iwm_Lu) 
     >         + wm_dt*( wm_conv(I,K,iwm_dirx)                 
     >         + wm_PrsGrad(I,K,iwm_dirx)                       
     >         - wm_diff(I,K,iwm_dirx) )
          ! this is the integrated momentum equation, except for the Lv term
          wm_lhs(I,K,iwm_diry) = -wm_inte(I,K,iwm_Lv) 
     >         + wm_dt*( wm_conv(I,K,iwm_diry)                 
     >         + wm_PrsGrad(I,K,iwm_diry)                       
     >         - wm_diff(I,K,iwm_diry) )
      end do
      end do
      

      end subroutine iwm_calc_lhs_rough     

      
!*******************************************************************************
      subroutine iwm_slv_rough(lhsx,lhsy,Ux,Uy,Dz,z0,utx,uty,fx,fy)
!*******************************************************************************
      INCLUDE 'PARAM.H'

      REAL vonk
      REAL lhsx, lhsy, Ux, Uy, Dz, z0, utx, uty
      REAL fx,fy
      REAL Ax, Ay, Vel, inteLu, inteLv
      
      vonk = 0.4

      Ax = (Ux-utx/vonk*log(Dz/z0))/((1.0-z0/Dz))
      Ay = (Uy-uty/vonk*log(Dz/z0))/((1.0-z0/Dz))
      Vel = sqrt(Ux**2.0+Uy**2.0)
      inteLu = 1.0/2.0*Dz*Ax*(1.0-z0/Dz)**2.0
     >            +1.0/vonk*utx*Dz*(z0/Dz-1.0+log(Dz/z0))
      inteLv = 1.0/2.0*Dz*Ay*(1.0-z0/Dz)**2.0
     >            +1.0/vonk*uty*Dz*(z0/Dz-1.0+log(Dz/z0))
      fx = inteLu+lhsx
      fy = inteLv+lhsy

      end subroutine iwm_slv_rough

!*******************************************************************************
      subroutine iwm_calc_wallstress_rough
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL vonk
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip      

      real fx, fy, fxp, fyp
      real wm_tol, wm_eps
      real a11, a12, a21, a22
      real wmutxP, wmutyP      
      real equilWMpara, equilutx, equiluty
      real wmpAx, wmpAy, wmputx, wmputy, wmpz0, wmpDz
      real utaup
      real dVelzT, dVelzB, Vel
      integer iter, MaxIter, equil_flag, div_flag
      
      vonk = 0.4

      MaxIter=1500

      wm_tol = 0.000001
      wm_eps = 0.000000001

      do i=1,n1m
      do k=1,n3m

      ! use Newton method to solve the system
c      wm_utx(i,k) = 1.0 * sign(1.0,wm_flt_tagvel(i,k,iwm_dirx))         
c      wm_uty(i,k) = 0.1 * sign(1.0,wm_flt_tagvel(i,k,iwm_diry))         

      call iwm_slv_rough(wm_lhs(i,k,iwm_dirx), wm_lhs(i,k,iwm_diry), 
     >     wm_flt_tagvel(i,k,iwm_dirx),                                 
     >     wm_flt_tagvel(i,k,iwm_diry),                                 
     >     wm_Dz(i,k), wm_z0(i,k),                              
     >     wm_utx(i,k), wm_uty(i,k),fx,fy )

c      if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c      write(*,*) wm_lhs(i,k,iwm_dirx), wm_lhs(i,k,iwm_diry), 
c     >     wm_flt_tagvel(i,k,iwm_dirx),                                 
c     >     wm_flt_tagvel(i,k,iwm_diry),                                 
c     >     wm_Dz(i,k), wm_z0(i,k),                              
c     >     wm_utx(i,k), wm_uty(i,k),fx,fy 
c	      end if

      iter = 0
      equil_flag = 0
      div_flag = 0
      do while (max(abs(fx),abs(fy))>wm_tol)

          wmutxP=wm_utx(i,k)+wm_eps
          wmutyP=wm_uty(i,k)
          call iwm_slv_rough(wm_lhs(i,k,iwm_dirx),                           
     >         wm_lhs(i,k,iwm_diry),wm_flt_tagvel(i,k,iwm_dirx),
     >         wm_flt_tagvel(i,k,iwm_diry),wm_Dz(i,k),          
     >         wm_z0(i,k), wmutxP, wmutyP, fxp, fyp )
          a11 = (fxp-fx)/wm_eps
	  a21 = 0.0
c          a21 = (fyp-fy)/wm_eps
c	  if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c          write(*,*) fxp,fyp,a11,a21
c          end if

          wmutxP = wm_utx(i,k)
          wmutyP = wm_uty(i,k)+wm_eps
          call iwm_slv_rough(wm_lhs(i,k,iwm_dirx),                            
     >         wm_lhs(i,k,iwm_diry),                                    
     >         wm_flt_tagvel(i,k,iwm_dirx),                             
     >         wm_flt_tagvel(i,k,iwm_diry),                             
     >         wm_Dz(i,k), wm_z0(i,k),                          
     >         wmutxP, wmutyP, fxp, fyp)
c          a12 = (fxp-fx)/wm_eps
	  a12 = 0.0
          a22 = (fyp-fy)/wm_eps

          wm_utx(i,k) = wm_utx(i,k)                            
     >         - ( a22*fx-a12*fy)/(a11*a22-a12*a21)
          wm_uty(i,k) = wm_uty(i,k)                            
     >         - (-a21*fx+a11*fy)/(a11*a22-a12*a21)
          call iwm_slv_rough(wm_lhs(i,k,iwm_dirx),                            
     >         wm_lhs(i,k,iwm_diry),                                    
     >         wm_flt_tagvel(i,k,iwm_dirx),                             
     >         wm_flt_tagvel(i,k,iwm_diry),                             
     >         wm_Dz(i,k), wm_z0(i,k),                          
     >         wm_utx(i,k), wm_uty(i,k), fx, fy)
          iter = iter+1
          
c	  if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c          write(*,*) wm_utx(i,k),wm_uty(i,k)
c          end if

          ! maximum iteration reached
          if (iter>MaxIter) then
              equil_flag = 1
              div_flag = 1;
              exit
          end if
      end do

      ! infinity check
      if (wm_utx(i,k)-1.0.eq.wm_utx(i,k) .or.                    
     >     wm_uty(i,k)-1.0.eq.wm_uty(i,k)) then
          equil_flag=1
          div_flag  =1
      end if



      ! calculate equilibrium us for equil_flag=1 use
      equilutx = vonk*wm_flt_tagvel(i,k,iwm_dirx)                       
     >     / log(wm_Dz(i,k)/wm_z0(i,k))
      equiluty = vonk*wm_flt_tagvel(i,k,iwm_diry)                       
     >     / log(wm_Dz(i,k)/wm_z0(i,k))
      if (equil_flag==1) then
          wm_utx(i,k) = equilutx
          wm_uty(i,k) = equiluty
      end if

      !calculate Ax, Ay
      if(equil_flag==1) then
          wmpAx = 0.0
          wmpAy = 0.0
      else
          ! eq. D2 in Yang et al. 2015
          wmpAx = ( wm_flt_tagvel(i,k,iwm_dirx)                        
     >         -wm_utx(i,k)/vonk*log(wm_Dz(i,k)                 
     >         /wm_z0(i,k)))                                             
     >         / ((1.0-wm_z0(i,k)/wm_Dz(i,k)))
          wmpAy = ( wm_flt_tagvel(i,k,iwm_diry)                        
     >         -wm_uty(i,k)/vonk*log(wm_Dz(i,k)                 
     >         /wm_z0(i,k)))                                            
     >         /((1.0-wm_z0(i,k)/wm_Dz(i,k)))
      end if



      ! check for excessive linear term correction
      ! after first 100 step this check is rarely invoked
c      if (abs(wmpAx)>1.0 .or. abs(wmpAy)>1.0) then
c          equil_flag = 1
c          wm_utx(i,k) = equilutx
c          wm_uty(i,k) = equiluty
c          wmpAx = 0.0
c          wmpAy = 0.0
c      end if

      ! store the linear correction
      wm_Ax(i,k) = wmpAx
      wm_Ay(i,k) = wmpAy

      ! update integral for last time step
      wm_inte_m(i,k,iwm_Lu ) = wm_inte(i,k,iwm_Lu )
      wm_inte_m(i,k,iwm_Lv ) = wm_inte(i,k,iwm_Lv )
      wm_inte_m(i,k,iwm_Luv) = wm_inte(i,k,iwm_Luv)
      wm_inte_m(i,k,iwm_Luu) = wm_inte(i,k,iwm_Luu)
      wm_inte_m(i,k,iwm_Lvv) = wm_inte(i,k,iwm_Lvv)

      !those temporary variables are used for convenient reference
      wmputx = wm_utx(i,k)
      wmputy = wm_uty(i,k)
      wmpDz  = wm_Dz (i,k)
      wmpz0  = wm_z0 (i,k)

      ! calculate the needed integrals

      ! Eq. D7 in Yang et al. 2015
      wm_inte(i,k,iwm_Lu) = 1.0/2.0*wmpDz*wmpAx                       
     >     *(1.0-wmpz0/wmpDz)**2.0                                       
     >     +1.0/vonk*wmputx*wmpDz*(wmpz0/wmpDz-1.0+log(wmpDz/wmpz0))
      wm_inte(i,k,iwm_Lv) = 1.0/2.0*wmpDz*wmpAy                       
     >     *(1.0-wmpz0/wmpDz)**2.0                                              
     >     +1.0/vonk*wmputy*wmpDz*(wmpz0/wmpDz-1.0+log(wmpDz/wmpz0))

      ! Eq. D8 in Yang et al 2015
      wm_inte(i,k,iwm_Luv) = 1.0/vonk**2.0*wmputx*wmputy*wmpDz       
     >     *(1.0-2*wmpz0/wmpDz+(1.0-log(wmpDz/wmpz0))**2.0)                  
     >     +1.0/3.0*wmpAx*wmpAy*wmpDz*(1.0 - wmpz0/wmpDz)**3.0               
     >     -1.0/4.0/vonk*(wmpAx*wmputy+wmpAy*wmputx)*wmpDz                   
     >     *(1.0-4.0*wmpz0/wmpDz+3.0*wmpz0**2.0/wmpDz**2.0                    
     >     -2.0*log(wmpDz/wmpz0)+4.0*wmpz0/wmpDz*log(wmpDz/wmpz0))
      wm_inte(i,k,iwm_Luu) = 1.0/vonk**2.0*wmputx**2.0*wmpDz          
     >     *((log(wmpDz/wmpz0)-1.0)**2.0-2.0*wmpz0/wmpDz+1.0)                 
     >     +1.0/3.0*wmpAx**2.0*wmpDz*(1.0-wmpz0/wmpDz)**3.0                   
     >     -1.0/2.0/vonk*wmputx*wmpAx*wmpDz                                    
     >     *(1.0-4.0*wmpz0/wmpDz+3.0*wmpz0**2.0/wmpDz**2.0                    
     >     -2.0*log(wmpDz/wmpz0)+4.0*wmpz0/wmpDz*log(wmpDz/wmpz0))

      ! Eq. D9 in Yang et al 2015
      wm_inte(i,k,iwm_Lvv) = 1.0/vonk**2.0*wmputy**2.0*wmpDz*         
     >     ((log(wmpDz/wmpz0)-1.0)**2.0-2.0*wmpz0/wmpDz+1.0)                  
     >     +1.0/3.0*wmpAy**2.0*wmpDz*(1.0-wmpz0/wmpDz)**3.0                   
     >     -1.0/2.0/vonk*wmputy*wmpAy*wmpDz                                    
     >     *(1.0-4.0*wmpz0/wmpDz-3.0*wmpz0**2.0/wmpDz**2.0                    
     >     -2.0*log(wmpDz/wmpz0)+4.0*wmpz0/wmpDz*log(wmpDz/wmpz0))

      ! calculate top and bottom derivatives
      ! Eq. D5 (a)
      wm_dudzT(i,k,iwm_dirx) = 1.0/wmpDz*(wmpAx+wmputx/vonk)
      wm_dudzT(i,k,iwm_diry) = 1.0/wmpDz*(wmpAy+wmputy/vonk)
      ! Eq. D5 (b)
      wm_dudzB(i,k,iwm_dirx) = 1.0/wmpDz*wmpAx+wmputx/vonk/wmpz0
      wm_dudzB(i,k,iwm_diry) = 1.0/wmpDz*wmpAy+wmputy/vonk/wmpz0

      ! calculte the turbulent diffusion term
      !total velocity
      Vel = sqrt(wm_flt_tagvel(i,k,iwm_dirx)**2.0                  
     >     +wm_flt_tagvel(i,k,iwm_diry)**2.0)
      ! Eq. D6
      dVelzT=abs(wm_flt_tagvel(i,k,iwm_dirx)/Vel                        
     >     *wm_dudzT(i,k,iwm_dirx)+wm_flt_tagvel(i,k,iwm_diry)  
     >     /Vel*wm_dudzT(i,k,iwm_diry))
      dvelzB = sqrt(wm_dudzB(i,k,iwm_dirx)**2.0                    
     >     +wm_dudzB(i,k,iwm_diry)**2.0)

      ! Eq. D4, the eddy viscosity is nu_T=(vonk*y)^2*dudy, hence the formula
      wm_diff(i,k,iwm_dirx) = (vonk*wmpDz)**2.0*dVelzT            
     >     *wm_dudzT(i,k,iwm_dirx)-(vonk*wmpz0)**2.0               
     >     *dVelzB*wm_dudzB(i,k,iwm_dirx)
      wm_diff(i,k,iwm_diry) = (vonk*wmpDz)**2.0*dVelzT            
     >     *wm_dudzT(i,k,iwm_diry)-(vonk*wmpz0)**2.0               
     >     *dVelzB*wm_dudzB(i,k,iwm_diry)

      ! calculate the wall stress
      if (equil_flag==1) then
          equilWMpara = sqrt(wm_flt_tagvel(i,k,iwm_dirx)**2.0      
     >         + wm_flt_tagvel(i,k,iwm_diry)**2.0)                  
     >         *vonk**2.0/(log(wm_Dz(i,k)                           
     >         /wm_z0(i,k)))**2.0
          wm_tauwx(i,k) = equilWMpara                                   
     >         *wm_flt_tagvel(i,k,iwm_dirx)
          wm_tauwy(i,k) = equilWMpara                                   
     >         *wm_flt_tagvel(i,k,iwm_diry)
      else
          ! Eq. D4
          wm_tauwx(i,k) = (vonk*wmpz0)**2.0*dVelzB                
     >         *wm_dudzB(i,k,iwm_dirx)
          wm_tauwy(i,k) = (vonk*iwmpz0)**2.0*dVelzB                
     >         *wm_dudzB(i,k,iwm_diry)
      end if

      ! calculate the friciton velocity
      ! definition of friction velocity
      utaup = (wm_tauwx(i,k)**2.0                                  
     >     +wm_tauwy(i,k)**2.0 )**0.25
      ! the filtered friction velocity used for filtering time scale
      wm_flt_us(i,k) = wm_flt_us(i,k)                          
     >     *(1.0-wm_tR(i,k))+utaup*wm_tR(i,k)

      ! update the filtering time scale
      ! Eq. 26
      wm_tR(i,k) = wm_dt/(wm_Dz(i,k)                          
     >     /wm_flt_us(i,k)/vonk)
      ! filtering time scale can only be larger than the time step,
      ! if not, then just use the instantaneous flow field to do the model
      if (wm_tR(i,k)>1.0) then
          wm_tR(i,k) = 1.0
      end if

      end do
      end do

      end subroutine iwm_calc_wallstress_rough

!*******************************************************************************
      subroutine iwm_monitor_rough(NTIME)
!*******************************************************************************
! This subroutine is to monitor the parameters at one point, do not call this
! subroutine if you are not interested in how the model works
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      INTEGER NTIME
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip

      integer iwm_i,iwm_j,dmpPrd,fid
      character*50 fname

      dmpPrd = iwm_ntime_skip
      iwm_i = int(n1m/2.0)
      iwm_j = int(n3m/2.0)
      write(fname,'(A,i5.5,A)') 'iwm_track.dat'
      open(1,file=fname, ACCESS='APPEND')
      if( mod(ntime,dmpPrd)==0) then
      write(1,*) ntime, wm_flt_tagvel(iwm_i,iwm_j,1), 
     >            wm_flt_tagvel(iwm_i,iwm_j,2),
     >            wm_utx(iwm_i,iwm_j),              
     >            wm_uty(iwm_i,iwm_j),  wm_Ax(iwm_i,iwm_j),                                
     >            wm_Ay(iwm_i,iwm_j), wm_tR(iwm_i,iwm_j)
      end if
      close(1)

      end subroutine iwm_monitor_rough
      

cccc  smooth-wall IWMLES model 
!*******************************************************************************
      subroutine wallmodel_iwmles_smooth (U,P,PRESG,PRESG3,TAUW,TAUB,
     >                                    SGSVIS,NTIME)
      ! BY RUI-FENG HU
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/ROUGHNESS/Y0
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL U1DOWN(M1,M3),W1DOWN(M1,M3)
      REAL DELTAP,TAUW1,TAUW2,LAMBDA
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL P(M1,M2,M3)
      REAL KAPA,B,TEMP(M1,M3),TEMP1
      REAL RETUW,U1,U2,REINV
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL PRESG,PRESG3
      REAL wm_dt
      INTEGER NTIME
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      
      REINV=1.0/DBLE(RE)
 
      KAPA = 0.4
      B    = 5.0 
      
CCCCC PARAMETERS FOR IWMLES
      iwm_dirx = 1  ! direction x
      iwm_diry = 2  ! direction z
      
      iwm_DN   = 2
      
      iwm_Lu   = 1   ! index for integral of u
      iwm_Luu  = 2   ! index for integral of uu
      iwm_Lv   = 3   ! etc.
      iwm_Lvv  = 4
      iwm_Luv  = 5
      
      iwm_LN   = 5  ! the total number of integrals that need to be calculated
      
      iwm_ntime_skip = 5
      
      IF (NTIME.EQ.1) CALL iwm_init_smooth(U)
      
      CALL iwm_wallstress_smooth(U,P,PRESG,PRESG3,TAUW,TAUB,NTIME)   

c     bottom wall
      DELTAP=DY(1)+DY(2)+DY(3)/2.0 
      WSM(1)=0.0
      
      DO K=1,N3M
      DO I=1,N1M
          
C      U_ave=SQRT(U(I,3,K,1)**2.+U(I,3,K,3)**2.)
C      TAUW1=(KAPA/ALOG(DELTAP/Y0)*U_ave)**2.0 
      TAUW1=SQRT(TAUW(I,K,1)**2.+TAUW(I,K,3)**2.)
      WSM(1)=WSM(1)+TAUW1
C      TAUW(I,K,1)=TAUW1*U(I,3,K,1)/U_ave
C      TAUB(I,K,1)=TAUW1*U(I,3,K,3)/U_ave
      
      ENDDO
      ENDDO
      
      WSM(1)=WSM(1)/REAL(N1M*N3M)      
      TAOW=WSM(1)
      
c      write(*,*) ntime, taow
      
      RETURN
      END

!******************************************************************************* 
      subroutine iwm_wallstress_smooth(U,P,PRESG,PRESG3,TAUW,TAUB,NTIME)      
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL PRESG,PRESG3   
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      
C      use sim_param , only : dudz, dvdz, txz, tyz
C      implicit none

      ! Calculate the time step used in the integral wall model
      !! DO NOT USE iwm_ntime_skip=1 !! !! this number is hard coded to prevent any
      !! mis-use...
      if (mod(NTIME,iwm_ntime_skip)==1) then
          wm_dt = DT
      else
          wm_dt = wm_dt+DT
      end if

      ! Compute the wall stress
      if(mod(NTIME,iwm_ntime_skip)==0) then
      ! gather flow status, update the integrated unsteady term, convective term,
      ! turbulent diffusion term etc.
      call iwm_calc_lhs_smooth(U,P,PRESG,PRESG3)
      ! the subroutine to calculate wall stress
      call iwm_calc_wallstress_smooth(U)
      ! this is to monitor any quantity from the iwm, useful debugging tool
      call iwm_monitor_smooth(NTIME)
      end if

      ! Imposing txz, tyz, dudz, dvdz every time step even iwm_* are not computed
      ! every time step.
      do I=1,N1M
      do K=1,N3M
          
      ! wall stress, use the value calculated in iwm, note the negative sign
      TAUW(I,K,1) = wm_tauwx(I,K)
      TAUW(I,K,3) = wm_tauwy(I,K)

      ! Note: Using the value calculated in iwm, note the positive sign
C      dudz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_dirx)
C      dvdz(iwm_i,iwm_j,1) = iwm_dudzT(iwm_i,iwm_j,iwm_diry)
      end do
      end do

      end subroutine iwm_wallstress_smooth


!*******************************************************************************
      subroutine iwm_init_smooth(U)
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/ROUGHNESS/Y0
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/PARA/RE
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL vonk,KAPA,B
      REAL usinit, uinit, vinit, Dzp
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      REAL U(0:M1,0:M2,0:M3,3)
      REAL REINV
      REAL u_tau1,u_tau2,futau3,futau1,futau2,lambda
      
      vonk = 0.4
      KAPA = 0.4
      B = 5.0
      REINV = 1.0/DBLE(RE)
      
      ! at the height of the first grid point.
      Dzp=DY(1)+DY(2)+DY(3)/2.0 
C      Dzp=dz/2._rprec
      
      ! initial value for the x-velocity at first grid point
c      uinit = usinit/vonk*log(Dzp/Y0)
      uinit = 0.0
      DO K=1,N3M
      DO I=1,N1M
      uinit = uinit + U(I,3,K,1)
      ENDDO
      ENDDO
      uinit = uinit/real(n1m*n3m)

      ! initial value for us (the friction velocity)
      U_tau1=1.0
      U_tau2=0.000001
      FUtau3=1.0
      FUtau1=0.1
      LAMBDA=1.0
     
      DO WHILE(abs(U_tau2-U_tau1).GT.1.0E-6)
      U_tau1=U_tau2
      FUtau1=uinit/U_tau1-ALOG(U_tau1*Dzp*RE)/KAPA-B
      FUtau2=-uinit/U_tau1**2.0-1./KAPA/U_tau1
      U_tau2=U_tau1-LAMBDA*FUtau1/FUtau2
      FUtau3=uinit/U_tau2-ALOG(U_tau2*Dzp*RE)/KAPA-B      
      IF (ABS(FUtau3).GT.ABS(FUtau1)) THEN
      LAMBDA=0.5*LAMBDA
      ENDIF      
      ENDDO 

      usinit= U_tau1

c      write(*,*) usinit
      
      ! initial value for the y-velocity at first grid point
      vinit = 0.0 
      
      
      DO K=1,N3M
      DO I=1,N1M
          
      ! us in x, y directions
      wm_utx(I,K) = usinit
      wm_uty(I,K) = 0.0
      
      ! wall stress in x, y directions
      wm_tauwx(I,K) = usinit**2.0
      wm_tauwy(I,K) = 0.0
      
      ! filitered velocity at the first grid point in x, y directions
      wm_flt_tagvel  (I,K,iwm_dirx) = uinit
      wm_flt_tagvel  (I,K,iwm_diry) = vinit
      wm_flt_tagvel_m(I,K,iwm_dirx) = uinit
      wm_flt_tagvel_m(I,K,iwm_diry) = vinit
      
      ! pressure at first grid point
      wm_flt_p(I,K) = 0.0
      
      ! integrals of Lu, Lv, etc.
      wm_inte(I,K,iwm_Lu)    = uinit*Dzp
      wm_inte(I,K,iwm_Lv)    = 0.0
      wm_inte(I,K,iwm_Luu)   = uinit**2.0*Dzp
      wm_inte(I,K,iwm_Lvv)   = 0.0
      wm_inte(I,K,iwm_Luv)   = 0.0
      wm_inte_m(I,K,iwm_Lu)  = uinit*Dzp
      wm_inte_m(I,K,iwm_Lv)  = 0.0
      wm_inte_m(I,K,iwm_Luu) = uinit**2.0*Dzp
      wm_inte_m(I,K,iwm_Lvv) = 0.0
      wm_inte_m(I,K,iwm_Luv) = 0.0
      
      ! each term in the integral equation and top/bottom derivatives
      wm_unsdy(I,K,iwm_dirx)   = 0.0
      wm_unsdy(I,K,iwm_diry)   = 0.0
      wm_conv(I,K,iwm_dirx)    = 0.0
      wm_conv(I,K,iwm_diry)    = 0.0
      wm_PrsGrad(I,K,iwm_dirx) = 0.0
      wm_PrsGrad(I,K,iwm_diry) = 0.0
      wm_diff(I,K,iwm_dirx)    = 0.0
      wm_diff(I,K,iwm_diry)    = 0.0
      wm_LHS(I,K,iwm_dirx)     = -uinit*Dzp
      wm_LHS(I,K,iwm_diry)     = -uinit*Dzp
      wm_dudzT(I,K,iwm_dirx) = usinit/vonk/Dzp
      wm_dudzT(I,K,iwm_diry) = 0.0
      wm_dudzB(I,K,iwm_dirx) = 0.0
      wm_dudzB(I,K,iwm_diry) = 0.0
      
      ! filtered friction velocity and the filtering time scale, tR<1
      wm_flt_us(I,K) = usinit
      wm_tR(I,K)     = (CFLMAX*ALX/N1M/uinit)/(Dzp/vonk/usinit)
      
      ! cell height and length scales
      wm_Dz(I,K)   = Dzp
      wm_deli(I,K) = min(11.0*REINV/usinit,Dzp)
      wm_delv(I,K) = REINV/usinit
c      if (i.eq.(n1m/2).and.k.eq.(n3m/2)) then
c      write(*,*) wm_Dz(I,K),wm_deli(I,K),wm_delv(I,K)
c      end if
      
      ! linear correction to the log profile
      wm_Ax(I,K) = 0.0
      wm_Ay(I,K) = 0.0
      wm_Cx(I,K) = 0.0
      wm_Cy(I,K) = 0.0
      
      ENDDO
      ENDDO       
      
      ! time step seen by the iwm
c      wm_dt=iwm_ntime_skip*CFLMAX*ALX/N1M/uinit  
      wm_dt=iwm_ntime_skip*DT
      
      end subroutine iwm_init_smooth

      
!*******************************************************************************
      subroutine iwm_calc_lhs_smooth(U,P,PRESG,PRESG3)
!*******************************************************************************
! Ths subroutine calculates the left hand side of the iwm system.
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/ROUGHNESS/Y0
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/PARA/RE
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL PRESG,PRESG3   
      REAL u_inst(0:M1,0:M3),v_inst(0:M1,0:M3),w_inst(0:M1,0:M3)
      REAL p_inst(0:M1,0:M3)
      REAL uux, uvx, uvy, vvy, ux, vy
      REAL phip, phim   
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip
      real wm_dt,p_bar
      real reinv

      REINV = 1.0/DBLE(RE)
      
      do I=1,N1M
      do K=1,N3M
          ! update the u, v for previous time step
          wm_flt_tagvel_m(I,K,iwm_dirx) =                     
     >         wm_flt_tagvel(I,K,iwm_dirx)
          wm_flt_tagvel_m(I,K,iwm_diry) =                     
     >         wm_flt_tagvel(I,K,iwm_diry)
          
          ! get the instantaneous field
          u_inst(I,K) = U(I,3,K,1)
          v_inst(I,K) = U(I,3,K,3)
          w_inst(I,K) = U(I,2,K,2)*0.25
          p_inst(I,K) = P(I,1,K)
      end do
      end do
      
      ! obtain the pressure fluctuations
      p_bar = 0.0
      do I = 1, N1M
      do K = 1, N3M
          p_bar = p_bar+p_inst(I,K)
      end do
      end do
      p_bar=p_bar/REAL(N1M)/REAL(N3M)

      do I = 1, N1M
      do K = 1, N3M
          p_inst(I,K) = p_inst(I,K)-p_bar
      end do
      end do

      ! all the data enters must be filtered (see Anderson & Meneveau 2011 JFM)
      !call test_filter ( u_inst )
      !call test_filter ( v_inst )
      !call test_filter ( w_inst )
      !call test_filter ( p_inst )

      !temporal filtering
      do I=1,N1M
      do K=1,N3M
          wm_flt_tagvel(I,K,iwm_dirx) =                                
     >     wm_flt_tagvel(I,K,iwm_dirx)*(1.0-wm_tR(I,K))    
     >         + u_inst(I,K)*wm_tR(I,K)
          wm_flt_tagvel(I,K,iwm_diry) =                                 
     >         wm_flt_tagvel(I,K,iwm_diry)*(1.0-wm_tR(I,K))    
     >         + v_inst(I,K)*wm_tR(I,K)
          wm_flt_p(I,K) = wm_flt_p(I,K)                            
     >         * (1.0-wm_tR(I,K))                                       
     >         + p_inst(I,K)*wm_tR(I,K)
      end do
      end do

      ! calculate LHS, calculation of the integrals is done from the last time step
      ! in the subroutine iwm_calc_wallstress, so is iwm_diff
      do K = 1, N3M
          KP=KPA(K)
          KM=KMA(K)
      do I = 1, N1M
          IP=IPA(I)
          IM=IMA(I)
          ! the unsteady term
          wm_unsdy(I,K,iwm_dirx) =                                      
     >         (wm_inte(I,K,iwm_Lu)-wm_inte_m(I,K,iwm_Lu))/wm_dt
          wm_unsdy(I,K,iwm_diry) =                                      
     >         (wm_inte(I,K,iwm_Lv)-wm_inte_m(I,K,iwm_Lv))/wm_dt

          ! the convective term
          phip = wm_inte(IP,K,iwm_Luu)
          phim = wm_inte(IM,K,iwm_Luu)
          uux = (phip-phim)*DX1/2.0
          phip = wm_inte(I,KP,iwm_Luv)
          phim = wm_inte(I,KM,iwm_Luv)
          uvy = (phip-phim)*DX3/2.0
          phip = wm_inte(IP,K,iwm_Luv)
          phim = wm_inte(IM,K,iwm_Luv)
          uvx = (phip-phim)*DX1/2.0
          phip = wm_inte(I,KP,iwm_Lvv)
          phim = wm_inte(I,KM,iwm_Lvv)
          vvy = (phip-phim)*DX3/2.0
          phip = wm_inte(IP,K,iwm_Lu )
          phim = wm_inte(IM,K,iwm_Lu )
          ux = (phip-phim)*DX1/2.0
          phip = wm_inte(I,KP,iwm_Lv )
          phim = wm_inte(I,KM,iwm_Lv )
          vy = (phip-phim)*DX3/2.0
          wm_conv(I,K,iwm_dirx) = uux + uvy                             
     >          - wm_flt_tagvel_m(I,K,iwm_dirx)*(ux+vy)
          wm_conv(I,K,iwm_diry) = uvx + vvy                           
     >          - wm_flt_tagvel_m(I,K,iwm_diry)*(ux+vy)

          ! the pressure gradient term
          phip = wm_flt_p(IP,K)
          phim = wm_flt_p(IM,K)
          wm_PrsGrad(I,K,iwm_dirx) = (phip-phim)*DX1/2.0                
     >         * wm_Dz(I,K) + PRESG*wm_Dz(I,K)
          phip = wm_flt_p(I,KP)
          phim = wm_flt_p(I,KM)
          wm_PrsGrad(I,K,iwm_diry) = (phip-phim)*DX3/2.0   
     >         * wm_Dz(I,K) + PRESG3*wm_Dz(I,K)

          ! the left hand side
          ! this is the integrated momentum equation, except for the Lu term
          wm_lhs(I,K,iwm_dirx) = -wm_inte(I,K,iwm_Lu) 
     >         + wm_dt*( wm_conv(I,K,iwm_dirx)                 
     >         + wm_PrsGrad(I,K,iwm_dirx)                       
     >         - wm_diff(I,K,iwm_dirx) )
          ! this is the integrated momentum equation, except for the Lv term
          wm_lhs(I,K,iwm_diry) = -wm_inte(I,K,iwm_Lv) 
     >         + wm_dt*( wm_conv(I,K,iwm_diry)                 
     >         + wm_PrsGrad(I,K,iwm_diry)                       
     >         - wm_diff(I,K,iwm_diry) )
      end do
      end do
      

      end subroutine iwm_calc_lhs_smooth 
      
      
!*******************************************************************************
      subroutine iwm_slv_smooth(lhsx,lhsy,Ux,Uy,Dz,utx,uty,fx,fy,i,k)
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M

      REAL vonk
      REAL lhsx, lhsy, Ux, Uy, Dz, utx, uty, us
      REAL fx,fy
      REAL Ax, Ay, Cx, Cy, Vel, inteLu, inteLv, delv, deli
      real reinv
      integer i,k
      
      vonk  = 0.4
      reinv = 1.0/DBLE(RE) 
      
      us   = (utx**4.0+uty**4.0)**0.25
c      delv = sqrt(reinv**2.0*(utx**2.0+uty**2.0)
c     >        /(utx**4.0+uty**4.0))
      delv = reinv/us
      deli = min(11.0*reinv/us,Dz)
      
      Vel = sqrt(Ux**2.0+Uy**2.0)
      Ax  = (Ux/us+Ux/Vel/vonk*log(deli/Dz)-deli/delv*utx/us)
     >        /(1.0-deli/Dz)
      Ay  = (Uy/us+Uy/Vel/vonk*log(deli/Dz)-deli/delv*uty/us)
     >        /(1.0-deli/Dz)
      Cx  = Ux/us - Ax
      Cy  = Uy/us - Ay      
      
      inteLu = 0.5*utx*deli**2.0/delv+us*Dz*(0.5*Ax
     >            *(1.0-deli**2.0/Dz**2.0)+Cx*(1.0-deli/Dz)
     >            -1.0/vonk*Ux/Vel*(1.0-deli/Dz+deli/Dz*log(deli/Dz)))
      inteLv = 0.5*uty*deli**2.0/delv+us*Dz*(0.5*Ay
     >            *(1.0-deli**2.0/Dz**2.0)+Cy*(1.0-deli/Dz)
     >            -1.0/vonk*Uy/Vel*(1.0-deli/Dz+deli/Dz*log(deli/Dz)))
      
      fx = inteLu+lhsx
      fy = inteLv+lhsy

c      if (i.eq.(n1m/2).and.k.eq.(n3m/2)) then
c      write(*,*) us,delv,deli,Ax,Ay,Cx,Cy
c      end if

      end subroutine iwm_slv_smooth

      
!*******************************************************************************
      subroutine iwm_calc_wallstress_smooth(U)
!*******************************************************************************
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      REAL vonk,reinv,kapa,B,deltap
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip      

      real fx, fy, fxp, fyp
      real wm_tol, wm_eps
      real*8 a11, a12, a21, a22
      real wmutxP, wmutyP      
      real equilWMpara, equilutx, equiluty
      real wmpAx, wmpAy, wmpCx, wmpCy, wmputx, wmputy, wmpDz
      real wmpdelv, wmpdeli
      real utaup
      real dVelzT, dVelzB, Vel
      integer iter, MaxIter, equil_flag, div_flag
      REAL u_tau1,u_tau2,futau3,futau1,futau2,lambda,u_ave
      
      REAL U(0:M1,0:M2,0:M3,3)
      
      vonk = 0.4
      kapa = 0.4
      B = 5.0
      reinv = 1.0/dble(RE)
      deltap = DY(1)+DY(2)+DY(3)/2.0

      MaxIter=1500

      wm_tol = 0.000001
      wm_eps = 0.000000001

      do i=1,n1m
      do k=1,n3m

      ! use Newton method to solve the system
c      wm_utx(i,k) = wm_utx(i,k) * sign(1.0,wm_flt_tagvel(i,k,iwm_dirx))         
c      wm_uty(i,k) = wm_uty(i,k) * sign(1.0,wm_flt_tagvel(i,k,iwm_diry))         

      call iwm_slv_smooth(wm_lhs(i,k,iwm_dirx), wm_lhs(i,k,iwm_diry), 
     >     wm_flt_tagvel(i,k,iwm_dirx),                                 
     >     wm_flt_tagvel(i,k,iwm_diry),                                 
     >     wm_Dz(i,k),                               
     >     wm_utx(i,k), wm_uty(i,k),fx,fy,i,k )

c      if (i.eq.(n1m/2).and.k.eq.(n3m/2)) then
c      write(*,*) wm_utx(i,k), wm_uty(i,k),fx,fy
c      end if

      iter = 0
      equil_flag = 0
      div_flag = 0
      do while (max(abs(fx),abs(fy))>wm_tol)

          wmutxP=wm_utx(i,k)+wm_eps
          wmutyP=wm_uty(i,k)
          call iwm_slv_smooth(wm_lhs(i,k,iwm_dirx),                           
     >         wm_lhs(i,k,iwm_diry),wm_flt_tagvel(i,k,iwm_dirx),
     >         wm_flt_tagvel(i,k,iwm_diry),wm_Dz(i,k),          
     >         wmutxP, wmutyP, fxp, fyp ,i,k)
          a11 = (fxp-fx)/wm_eps
          a21 = (fyp-fy)/wm_eps
c	  if (i.eq.int(n1m/2).and.k.eq.int(n3m/2)) then
c          write(*,*) fx,fy,fxp,fyp,a11,a21
c          end if

          wmutxP = wm_utx(i,k)
          wmutyP = wm_uty(i,k)+wm_eps
          call iwm_slv_smooth(wm_lhs(i,k,iwm_dirx),                            
     >         wm_lhs(i,k,iwm_diry),                                    
     >         wm_flt_tagvel(i,k,iwm_dirx),                             
     >         wm_flt_tagvel(i,k,iwm_diry),                             
     >         wm_Dz(i,k),                           
     >         wmutxP, wmutyP, fxp, fyp,i,k)
          a12 = (fxp-fx)/wm_eps
          a22 = (fyp-fy)/wm_eps

          wm_utx(i,k) = wm_utx(i,k)                            
     >         - ( a22*fx-a12*fy)/(a11*a22-a12*a21)
          wm_uty(i,k) = wm_uty(i,k)                            
     >         - (-a21*fx+a11*fy)/(a11*a22-a12*a21)

c	  wm_utx(i,k)=abs(wm_utx(i,k))
c	  wm_uty(i,k)=abs(wm_uty(i,k))

          call iwm_slv_smooth(wm_lhs(i,k,iwm_dirx),                            
     >         wm_lhs(i,k,iwm_diry),                                    
     >         wm_flt_tagvel(i,k,iwm_dirx),                             
     >         wm_flt_tagvel(i,k,iwm_diry),                             
     >         wm_Dz(i,k),                           
     >         wm_utx(i,k), wm_uty(i,k), fx, fy,i,k)
          iter = iter+1
          
c	  if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c          write(*,*) wm_utx(i,k), wm_uty(i,k), fx, fy
c          end if

          ! maximum iteration reached
          if (iter>MaxIter) then
              equil_flag = 1
              div_flag = 1;
              exit
          end if
      end do

      ! infinity check
      if (wm_utx(i,k)-1.0.eq.wm_utx(i,k) .or.                    
     >     wm_uty(i,k)-1.0.eq.wm_uty(i,k)) then
          equil_flag=1
          div_flag  =1
      end if


      ! calculate equilibrium us for equil_flag=1 use
c      equilutx = vonk*wm_flt_tagvel(i,k,iwm_dirx)                       
c     >     / log(wm_Dz(i,k)/wm_z0(i,k))
c      equiluty = vonk*wm_flt_tagvel(i,k,iwm_diry)                       
c     >     / log(wm_Dz(i,k)/wm_z0(i,k))
      
      if (equil_flag.eq.1) then
          
          U_ave=SQRT(wm_flt_tagvel(i,k,iwm_dirx)**2.
     >		+wm_flt_tagvel(i,k,iwm_diry)**2.)     
     
          U_tau1=1.0
          U_tau2=0.000001
          FUtau3=1.0
          FUtau1=0.1
          LAMBDA=1.0
     
          DO WHILE(abs(U_tau2-U_tau1).GT.1.0E-6)
          U_tau1=U_tau2
          FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
          FUtau2=-U_ave/U_tau1**2-1./KAPA/U_tau1
          U_tau2=U_tau1-LAMBDA*FUtau1/FUtau2
          FUtau3=U_ave/U_tau2-ALOG(U_tau2*DELTAP*RE)/KAPA-B      
          IF (ABS(FUtau3).GT.ABS(FUtau1)) THEN
          LAMBDA=0.5*LAMBDA
          ENDIF      
          ENDDO
     
          TAUW1=U_tau1**2

c	  if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c          write(*,*) u_tau1
c          end if
 
          equilutx=TAUW1*wm_flt_tagvel(i,k,iwm_dirx)/U_ave
          equiluty=TAUW1*wm_flt_tagvel(i,k,iwm_diry)/U_ave
          
          wm_utx(i,k) = sqrt(abs(equilutx))
          wm_uty(i,k) = sqrt(abs(equiluty))

c	  if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c          write(*,*) equilutx,equiluty,wm_utx(i,k),wm_uty(i,k)
c          end if
          
      end if

      !those temporary variables are used for convenient reference
      wmputx  = wm_utx(i,k)
      wmputy  = wm_uty(i,k)
      wmpDz   = wm_Dz (i,k)      
      utaup   = (wmputx**4.0+wmputy**4.0)**0.25
c      wmpdelv = sqrt(reinv**2.0*(wmputx**2.0+wmputy**2.0)
c     >                /(wmputx**4.0+wmputy**4.0))
      wmpdelv = reinv/utaup
      wmpdeli = min(11.0*reinv/utaup,wmpDz)
      wmpUx   = wm_flt_tagvel(i,k,iwm_dirx)
      wmpUy   = wm_flt_tagvel(i,k,iwm_diry)
      wmpVel  = sqrt(wmpUx**2.0+wmpUy**2.0)

      
      !calculate Ax, Ay, Cx, Cy
      if(equil_flag==1)then
          wmpAx = 0.0
          wmpAy = 0.0
          wmpCx = 0.0
          wmpCy = 0.0
      else
          wmpAx = (wmpUx/utaup+wmpUx/wmpVel/vonk*log(wmpdeli/wmpDz)
     >            -wmpdeli/wmpdelv*wmputx/utaup)
     >            /(1.0-wmpdeli/wmpDz)
          wmpAy = (wmpUy/utaup+wmpUy/wmpVel/vonk*log(wmpdeli/wmpDz)
     >            -wmpdeli/wmpdelv*wmputy/utaup)
     >            /(1.0-wmpdeli/wmpDz)
          wmpCx = wmpUx/utaup - wmpAx
          wmpCy = wmpUy/utaup - wmpAy
      end if

c      if (i.eq.int(n1m/2).and.k.eq.(n3m/2)) then
c      write(*,*) wmpAx,wmpAy,wmpCx,wmpCy
c      end if


      ! check for excessive linear term correction
      ! after first 100 step this check is rarely invoked
      if (abs(wmpAx)>1.0 .or. abs(wmpAy)>1.0) then
          equil_flag = 1
          wm_utx(i,k) = sqrt(abs(equilutx))
          wm_uty(i,k) = sqrt(abs(equiluty))
          wmpAx = 0.0
          wmpAy = 0.0
          wmpCx = 0.0
          wmpCy = 0.0
      end if

      ! store the linear correction
      wm_Ax(i,k) = wmpAx
      wm_Ay(i,k) = wmpAy
      wm_Cx(i,k) = wmpCx
      wm_Cy(i,k) = wmpCy

      ! update integral for last time step
      wm_inte_m(i,k,iwm_Lu ) = wm_inte(i,k,iwm_Lu )
      wm_inte_m(i,k,iwm_Lv ) = wm_inte(i,k,iwm_Lv )
      wm_inte_m(i,k,iwm_Luv) = wm_inte(i,k,iwm_Luv)
      wm_inte_m(i,k,iwm_Luu) = wm_inte(i,k,iwm_Luu)
      wm_inte_m(i,k,iwm_Lvv) = wm_inte(i,k,iwm_Lvv)

      ! calculate the needed integrals

      ! Eq. C19 in Yang et al. 2015
      wm_inte(i,k,iwm_Lu) = 0.5*wmputx*wmpdeli**2.0/wmpdelv
     >     +utaup*wmpDz*(0.5*wmpAx*(1.0-wmpdeli**2.0/wmpDz**2.0)
     >     +wmpCx*(1.0-wmpdeli/wmpDz)-1.0/vonk*wmpUx/wmpVel
     >     *(1.0-wmpdeli/wmpDz+wmpdeli/wmpDz*log(wmpdeli/wmpDz)))
      wm_inte(i,k,iwm_Lv) = 0.5*wmputy*wmpdeli**2.0/wmpdelv
     >     +utaup*wmpDz*(0.5*wmpAy*(1.0-wmpdeli**2.0/wmpDz**2.0)
     >     +wmpCy*(1.0-wmpdeli/wmpDz)-1.0/vonk*wmpUy/wmpVel
     >     *(1.0-wmpdeli/wmpDz+wmpdeli/wmpDz*log(wmpdeli/wmpDz)))

      ! Eq. C21 in Yang et al 2015
      wm_inte(i,k,iwm_Luv) = 1.0/3.0*wmputx*wmputy*wmpdeli**3.0
     >    /wmpdelv**2.0+utaup**2.0*wmpDz*(-1.0/vonk*(wmpAx*wmpUy/wmpVel
     >    +wmpAy*wmpUx/wmpVel)*(0.25-0.25*wmpdeli**2.0/wmpDz**2.0
     >    +0.5*wmpdeli**2.0/wmpDz*log(wmpdeli/wmpDz))
     >    -1.0/vonk*(wmpCx*wmpUy/wmpVel+wmpCy*wmpUx/wmpVel)
     >    *(1.0-wmpdeli/wmpDz+wmpdeli/wmpDz*log(wmpdeli/wmpDz))
     >    -wmpUx*wmpUy/wmpVel**2.0/vonk**2.0*(wmpdeli/wmpDz-2.0
     >    +wmpdeli/wmpDz*(log(wmpdeli/wmpDz)-1.0)**2.0)
     >    +1.0/3.0*wmpAx*wmpAy*(1.0-wmpdeli**3.0/wmpDz**3.0)
     >    +0.5*(wmpAx*wmpCy+wmpAy*wmpCx)*(1.0-wmpdeli**2.0/wmpDz**2.0)
     >    +wmpCx*wmpCy*(1.0-wmpdeli/wmpDz))
      
      ! Eq. C20 in Yang et al 2015
      wm_inte(i,k,iwm_Luu) = 1.0/3.0*wmputx**2.0*wmpdeli**3.0
     >    /wmpdelv**2.0+utaup**2.0*wmpDz*(-wmpAx*wmpUx/wmpVel/vonk
     >    *wmpdeli**2.0/wmpDz**2.0*log(wmpdeli/wmpDz)
     >    +wmpAx*(wmpCx-1.0/2.0/vonk*wmpUx/wmpVel)
     >    *(1.0-wmpdeli**2.0/wmpDz**2.0)
     >    +wmpAx**2.0/3.0*(1.0-wmpdeli**3.0/wmpDz**3.0)
     >    +(wmpCx-1.0/vonk*wmpUx/wmpVel)**2.0-wmpdeli/wmpDz
     >    *(wmpCx-wmpUx/wmpVel/vonk+wmpUx/wmpVel/vonk
     >    *log(wmpdeli/wmpDz))**2.0+1.0/vonk**2.0*wmpUx**2.0
     >    /wmpVel**2.0*(1.0-wmpdeli/wmpDz))      
      wm_inte(i,k,iwm_Lvv) = 1.0/3.0*wmputy**2.0*wmpdeli**3.0
     >    /wmpdelv**2.0+utaup**2.0*wmpDz*(-wmpAy*wmpUy/wmpVel/vonk
     >    *wmpdeli**2.0/wmpDz**2.0*log(wmpdeli/wmpDz)
     >    +wmpAy*(wmpCy-1.0/2.0/vonk*wmpUy/wmpVel)
     >    *(1.0-wmpdeli**2.0/wmpDz**2.0)
     >    +wmpAy**2.0/3.0*(1.0-wmpdeli**3.0/wmpDz**3.0)
     >    +(wmpCy-1.0/vonk*wmpUy/wmpVel)**2.0-wmpdeli/wmpDz
     >    *(wmpCy-wmpUy/wmpVel/vonk+wmpUy/wmpVel/vonk
     >    *log(wmpdeli/wmpDz))**2.0+1.0/vonk**2.0*wmpUy**2.0
     >    /wmpVel**2.0*(1.0-wmpdeli/wmpDz)) 
      
      ! calculate top and bottom derivatives
      if (wmpdeli.lt.wmpDz) then
          wm_dudzT(i,k,iwm_dirx) = utaup/wmpDz*(wmpAx
     >            +1.0/vonk*wmpUx/wmpVel)
          wm_dudzT(i,k,iwm_diry) = utaup/wmpDz*(wmpAy
     >            +1.0/vonk*wmpUy/wmpVel)
      else
          wm_dudzT(i,k,iwm_dirx) = wmputx/wmpdelv
          wm_dudzT(i,k,iwm_diry) = wmputy/wmpdelv
      end if
      
      wm_dudzB(i,k,iwm_dirx) = wmputx/wmpdelv
      wm_dudzB(i,k,iwm_diry) = wmputy/wmpdelv

      ! calculte the turbulent diffusion term
      dVelzT=abs(wmpUx/wmpVel*wm_dudzT(i,k,iwm_dirx)                        
     >     +wmpUy/wmpVel*wm_dudzT(i,k,iwm_diry))  
      dvelzB = sqrt(wm_dudzB(i,k,iwm_dirx)**2.0                    
     >     +wm_dudzB(i,k,iwm_diry)**2.0)

      ! the eddy viscosity is nu_T=(vonk*y)^2*dudy, hence the formula
      if (wmpdeli.lt.wmpDz) then
          wm_diff(i,k,iwm_dirx) = (reinv+(vonk*wmpDz)**2.0*dVelzT)            
     >        *wm_dudzT(i,k,iwm_dirx)-reinv*wm_dudzB(i,k,iwm_dirx)
          wm_diff(i,k,iwm_diry) = (reinv+(vonk*wmpDz)**2.0*dVelzT)            
     >        *wm_dudzT(i,k,iwm_diry)-reinv*wm_dudzB(i,k,iwm_diry)
      else
          wm_diff(i,k,iwm_dirx) = 0.0
          wm_diff(i,k,iwm_diry) = 0.0
      end if

      ! calculate the wall stress
      if (equil_flag==1) then          
          wm_tauwx(i,k) = equilutx
          wm_tauwy(i,k) = equiluty
      else
          if (wmpdeli.lt.wmpDz) then
              wm_tauwx(i,k) = wmputx**2.0
              wm_tauwy(i,k) = wmputy**2.0
          else
              wm_tauwx(i,k) = reinv*wmpUx/wmpDz
              wm_tauwy(i,k) = reinv*wmpUy/wmpDz
          end if
      end if
      
      ! the filtered friction velocity used for filtering time scale
      wm_flt_us(i,k) = wm_flt_us(i,k)                          
     >     *(1.0-wm_tR(i,k))+utaup*wm_tR(i,k)

      ! update the filtering time scale
      ! Eq. 26
      wm_tR(i,k) = wm_dt/(wm_Dz(i,k)                          
     >     /wm_flt_us(i,k)/vonk)
      ! filtering time scale can only be larger than the time step,
      ! if not, then just use the instantaneous flow field to do the model
      if (wm_tR(i,k)>1.0) then
          wm_tR(i,k) = 1.0
      end if

      end do
      end do

      end subroutine iwm_calc_wallstress_smooth



!*******************************************************************************
      subroutine iwm_monitor_smooth(NTIME)
!*******************************************************************************
! This subroutine is to monitor the parameters at one point, do not call this
! subroutine if you are not interested in how the model works
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      
      COMMON/IWMLES/wm_utx(M1,M3), wm_uty(M1,M3),wm_tauwx(M1,M3), 
     >            wm_tauwy(M1,M3),wm_flt_tagvel(M1,M3,2), 
     >            wm_flt_tagvel_m(M1,M3,2),wm_flt_p(M1,M3),
     >            wm_inte(M1,M3,5), wm_inte_m(M1,M3,5),
     >            wm_unsdy(M1,M3,2), wm_conv(M1,M3,2), 
     >            wm_PrsGrad(M1,M3,2),wm_diff(M1,M3,2), 
     >            wm_LHS(M1,M3,2), wm_dudzT(M1,M3,2), wm_dudzB(M1,M3,2),
     >            wm_flt_us(M1,M3), wm_tR(M1,M3), wm_Dz(M1,M3), 
     >            wm_z0(M1,M3), wm_Ax(M1,M3), wm_Ay(M1,M3),
     >            wm_delv(M1,M3),wm_deli(M1,M3),
     >            wm_Cx(M1,M3),wm_Cy(M1,M3) 
      
      COMMON/IWMLES_PARA/iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv,
     >                    iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip,wm_dt
      
      INTEGER NTIME
      INTEGER iwm_dirx,iwm_diry,iwm_DN,iwm_Lu,iwm_Luu,iwm_Lv
      INTEGER iwm_Lvv,iwm_Luv,iwm_LN,iwm_ntime_skip

      integer iwm_i,iwm_j,dmpPrd,fid
      character*50 fname

      dmpPrd = iwm_ntime_skip
      iwm_i = int(n1m/2.0)
      iwm_j = int(n3m/2.0)
      write(fname,'(A,i5.5,A)') 'iwm_track.dat'
      open(1,file=fname, ACCESS='APPEND')
      if( mod(ntime,dmpPrd)==0) then
      write(1,*) ntime, wm_flt_tagvel(iwm_i,iwm_j,1), 
     >            wm_flt_tagvel(iwm_i,iwm_j,2),
     >            wm_utx(iwm_i,iwm_j),              
     >            wm_uty(iwm_i,iwm_j),  wm_Ax(iwm_i,iwm_j),                                
     >            wm_Ay(iwm_i,iwm_j), wm_Cx(iwm_i,iwm_j),                                
     >            wm_Cy(iwm_i,iwm_j), wm_tR(iwm_i,iwm_j)
      end if
      close(1)

      end subroutine iwm_monitor_smooth

C********************** WALLMODEL_MARUSIC1 ***************************
C     THIS SUBROUTINE IS TO NEAR WALL MODEL.  !by wang 0910
C     EQUILIBRIUM LAWS OF MARUSIC (2014) Txy,w=<Tw>+0.1*(u(x,y(OL),z)-u(x,y,z))
C                                         Tyz,w=<Tw>*w(x+delts)/<u>

      SUBROUTINE WALLMODEL_MARUSIC1(U,TAUW,TAUB,SGSVIS)  !BY WANG 0901
      INCLUDE 'PARAM.H'
      
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
!      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)     
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M   
      INTEGER IP,JP,KP
      INTEGER I,J,K
      INTEGER ZS
      REAL  U1UP,U1DOWN
      REAL  DELTAP,TAUW1,TAUW2,LAMBDA
      REAL  U_ave,U_tau1,U_tau2,FUtau1,FUtau2,FUtau3
      REAL  RETUW,U1,U2,REINV
      REAL  UTAU1,UTAU2
      REAL  KAPA,B,TEMP,DTEMP,UMTEMP
	REAL  XP,YP,ZP
	REAL  UF,VF,WF
      REAL  ZO,B1,B2,B3,C
     
      
      REAL  U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL  TAUWW(0:M1,0:M3,4)
      REAL  TAUW(0:M1,0:M3,2)
      REAL  TAUB(0:M1,0:M3,2)
     
      !SOLVED BY NEWTON ITERATIONS      
      KAPA=0.41
      
      IFWALL=1.0
      PI=ACOS(-1.0)
      ZO=0.06  ! THE MIDDLE OF THE LOG-LAYER
      
      
      DO J=1,N2M
      IF((ABS(Y(J)-ZO)).LE.(ABS(Y(J+1)-ZO)).AND.
     >     (ABS(Y(J)-ZO)).LE.(ABS(Y(J-1)-ZO))) THEN
      ZS = J
      ENDIF
      ENDDO     !THE LOCATION OF OUTER SIGNAL   
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
      U1DOWN=0.0
      U1UP=0.0
      REINV=1.0/RE
      DO J=1,N3M
      DO I=1,N1M 
      U1DOWN=U1DOWN+U(I,ZS,J,1)  
      U1UP=U1UP+U(I,N2M-ZS+1,J,1) 
      ENDDO
      ENDDO         

      U1DOWN=U1DOWN/DBLE(N1M*N3M)      !THE OUTER MEAN STREAMWISE VELOCITY
      U1UP=U1UP/DBLE(N1M*N3M)
      
      DELTAP= YC(ZS)    
      U_tau1=1.0
      U_tau2=0.000001
      U_ave= U1DOWN
      B=5.2 
      
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
      U_tau2=0.5*U_tau2+0.5*U_tau1
      ELSE
      U_tau1=0.5*U_tau2+0.5*U_tau1
      ENDIF
      ENDDO
      
      UTAU2=U_tau1
      TAOW=(UTAU2)**2
      WSM(1)= TAOW
      TAUW1=(UTAU2)**2
      
      YPLUSML=DELTAP*UTAU2*RE
      YPLUS1=DY(1)/2.0*UTAU2*RE
      C=1.0 + 0.05 * (ALOG(YPLUS1/YPLUSML))**2
      DELTS=DELTAP/TAN(16.0*PI/180)
      
      YP=DELTAP
      DO I=1,N1M
      DO K=1,N3M
      XP=XG(I)+DELTS
      IF(XP.GT.ALX) XP= XP-ALX
      ZP=ZC(K)
      CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
      CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
      
      TAUW(I,K,1)=TAUW1+C*0.1*UTAU2*(UF-U1DOWN)
      
      XP=XC(I)+DELTS
      IF(XP.GT.ALX) XP= XP-ALX
      ZP=ZG(K)
      CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
      CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
      TAUB(I,K,1)=TAUW1*WF/U1DOWN
      ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      DELTAP=YC(ZS)     
      
      U_tau1=1.0
      U_tau2=0.000001
      U_ave= U1UP
      
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
      U_tau2=0.5*U_tau2+0.5*U_tau1
      ELSE
      U_tau1=0.5*U_tau2+0.5*U_tau1
      ENDIF
      ENDDO
      
      UTAU2=U_tau1
      TAUW1=U_tau1**2     
     
      YPLUSML=DELTAP*UTAU2*RE
      YPLUS1=DY(1)/2.0*UTAU2*RE
      C=1.0 + 0.05 * (ALOG(YPLUS1/YPLUSML))**2
      DELTS=DELTAP/TAN(16.0*PI/180)
 
      YP=YC(N2M-ZS+1)
      DO I=1,N1M
      DO K=1,N3M
      XP=XG(I)+DELTS
      IF(XP.GT.ALX) XP= XP-ALX
      ZP=ZC(K)
      CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
      CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
      TAUW(I,K,2)=TAUW1+C*0.1*UTAU2*(UF-U1UP)
      XP=XC(I)+DELTS
      IF(XP.GT.ALX) XP= XP-ALX
      ZP=ZG(K)
      CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
      CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
      TAUB(I,K,2)=TAUW1*WF/U1UP
      ENDDO
      ENDDO    
      RETURN
      END
      
      
C********************** WALLMODEL_PIOMELLI ***************************
C     THIS SUBROUTINE IS TO NEAR WALL MODEL.  !by wang 0910
C     EQUILIBRIUM LAWS OF PIOMELLI (1989) Txy,w=<Tw>*u(x+delts)/<u>
C                                         Tyz,w=<Tw>*w(x+delts)/<u>


      SUBROUTINE WALLMODEL_GROTZBACH(U,TAUW,TAUB,SGSVIS)   !BY WANG 0901
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
!      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
      COMMON/WALLTAO/TAOW
      COMMON/WSMEM/WSM(3),WSMO(3)   
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M   
      COMMON/RANSVISO/RANSVIS,RANSVIS2
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      REAL  U1UP,U1DOWN
      REAL DELTAP,TAUW1,TAUW2
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2
      REAL  RETUW,U1,U2,REINV
      REAL UTAU1,UTAU2
      INTEGER I,K 
      REAL KAPA,B,TEMP,DTEMP,UMTEMP
      INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UF,VF,WF
   
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL   TAUWW(0:M1,0:M3,4)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      KAPA=0.41
      B=5.2
      IFWALL=1.0
      N2M=M2-1  
      PI=ACOS(-1.0)
      U1DOWN=0.0
      U1UP=0.0
      REINV=1.0/RE
      DO J=1,N3M
      DO I=1,N1M 
      U1DOWN=U1DOWN+U(I,1,J,1)
      U1UP=U1UP+U(I,N2M,J,1)      
      ENDDO
      ENDDO
      U1DOWN=U1DOWN/DBLE(N1M*N3M)
      U1UP=U1UP/DBLE(N1M*N3M)
	      
      DELTAP=DY(1)/2.0             
      U_ave= ABS(U1DOWN)
      U_tau2=0.6
      U_tau1=0.000001
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
        TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
       U_tau2=0.5*U_tau2+0.5*U_tau1
       ELSE
       U_tau1=0.5*U_tau2+0.5*U_tau1
       ENDIF
      ENDDO
      UTAU2=U_tau1
      
      TAOW=(UTAU2)**2
      WSM(1)= TAOW
      
      TAUW1=(UTAU2)**2/ABS(U1DOWN)
      TAUB1=2./RE/DY(1) 
        YPLUS=DELTP*UTAU2*RE
        IF(YPLUS.GT.30.0.AND.YPLUS.LT.50.0)THEN
         DELTS=DELTP/TAN(8.0*PI/180)
        ELSE IF(YPLUS.GE.50)THEN
         DELTS=DELTP/TAN(13.0*PI/180)
        ELSE
         DELTS=0.0
        ENDIF
 
	  YP=DY(1)
       DO I=1,N1M
       DO K=1,N3M
           XP=XG(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZC(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
         TAUW(I,K,1)=TAUW1*UF
           XP=XC(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZG(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
       TAUB(I,K,1)=TAUB1*WF
      ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      DELTAP=DY(N2M)/2.0             
      U_ave= ABS(U1UP)
      U_tau2=0.2
      U_tau1=0.000001
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
         U_tau2=0.5*U_tau2+0.5*U_tau1
       ELSE
       U_tau1=0.5*U_tau2+0.5*U_tau1
       ENDIF
      ENDDO

      UTAU2=U_tau1
        
      TAUW1=(UTAU2)**2/U_ave 
         TAUB1=2./RE/DY(N2M)    
        YPLUS=DELTP*UTAU2*RE
        IF(YPLUS.GT.30.0.AND.YPLUS.LT.50.0)THEN
         DELTS=DELTP/TAN(8.0*PI/180)
        ELSE IF(YPLUS.GE.50)THEN
         DELTS=DELTP/TAN(13.0*PI/180)
        ELSE
         DELTS=0.0
        ENDIF
 
	  YP=YC(N2M)
       DO I=1,N1M
       DO K=1,N3M
           XP=XG(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZC(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
         TAUW(I,K,1)=TAUW1*UF
           XP=XC(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZG(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
       TAUB(I,K,1)=TAUB1*WF
      ENDDO
      ENDDO    


       !     WRITE(*,*)UTAU2,WSM(1)     
      RETURN
      END

C********************** WALLMODEL_MARUSIC ***************************
C     THIS SUBROUTINE IS TO NEAR WALL MODEL.  !by wang 0910
C     EQUILIBRIUM LAWS OF MARUSIC (2001) Txy,w=<Tw>+0.1*(u(x+delts)-u(x))
C                                         Tyz,w=<Tw>*w(x+delts)/<u>


      SUBROUTINE WALLMODEL_MARUSIC(U,TAUW,TAUB,SGSVIS)  !BY WANG 0901
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
!      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
       COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)     
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M   
      REAL  U1UP,U1DOWN
      REAL DELTAP,TAUW1,TAUW2
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2
      REAL  RETUW,U1,U2,REINV
      REAL UTAU1,UTAU2
      REAL KAPA,B,TEMP,DTEMP,UMTEMP
      INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UF,VF,WF
   
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL   TAUWW(0:M1,0:M3,4)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      INTEGER I,K 
      !SOLVED BY NEWTON ITERATIONS
      
      KAPA=0.41
      B=5.2
      IFWALL=1.0
      PI=ACOS(-1.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
      U1DOWN=0.0
      U1UP=0.0
      REINV=1.0/RE
      DO J=1,N3M
      DO I=1,N1M 
      U1DOWN=U1DOWN+U(I,1,J,1)
      U1UP=U1UP+U(I,N2M,J,1)      
      ENDDO
      ENDDO

      U1DOWN=U1DOWN/DBLE(N1M*N3M)
      U1UP=U1UP/DBLE(N1M*N3M)


      DELTAP=DY(1)/2.0             
      U_ave= U1DOWN
      U_tau2=0.2
      U_tau1=0.000001
       B=5.2   
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
       TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
         U_tau2=0.5*U_tau2+0.5*U_tau1
       ELSE
       U_tau1=0.5*U_tau2+0.5*U_tau1
       ENDIF
      ENDDO
      UTAU2=U_tau1
      TAOW=(UTAU2)**2
      WSM(1)= TAOW
      
      TAUW1=(UTAU2)**2
        YPLUS=DELTP*UTAU2*RE
        IF(YPLUS.GT.30.0.AND.YPLUS.LT.50.0)THEN
         DELTS=DELTP/TAN(8.0*PI/180)
        ELSE IF(YPLUS.GE.50)THEN
         DELTS=DELTP/TAN(13.0*PI/180)
        ELSE
         DELTS=0.0
        ENDIF
 
	  YP=H(1)
       DO I=1,N1M
       DO K=1,N3M
           XP=XG(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZC(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
         TAUW(I,K,1)=TAUW1+0.1*UTAU2*(UF-U1DOWN)
           XP=XC(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZG(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
       TAUB(I,K,1)=TAUW1*WF/U1DOWN
      ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      DELTAP=DY(N2M)/2.0             
      U_ave= U1UP
      U_tau2=0.2
      U_tau1=0.000001
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
         U_tau2=0.5*U_tau2+0.5*U_tau1
       ELSE
       U_tau1=0.5*U_tau2+0.5*U_tau1
       ENDIF
   
      ENDDO

      UTAU2=U_tau1
      
      TAUW1=(UTAU2)**2
     
        YPLUS=DELTP*UTAU1*RE
        IF(YPLUS.GT.30.0.AND.YPLUS.LT.50.0)THEN
         DELTS=DELTP/TAN(8.0*PI/180)
        ELSE IF(YPLUS.GE.50)THEN
         DELTS=DELTP/TAN(13.0*PI/180)
        ELSE
         DELTS=0.0
        ENDIF
 
	  YP=YC(N2M)
       DO I=1,N1M
       DO K=1,N3M
           XP=XG(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZC(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
        TAUW(I,K,2)=TAUW1+0.1*UTAU2*(UF-U1UP)
        XP=XC(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZG(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
       TAUB(I,K,2)=TAUW1*WF/U1UP
      ENDDO
      ENDDO    
      RETURN
      END
C********************** WALLMODEL_WERNER ***************************


      SUBROUTINE WALLMODEL_WERNER(U,TAUW,TAUB,SGSVIS,NTIME) 
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2) 
       COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)

      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLTAO/TAOW

      COMMON/RANSVISO/RANSVIS,RANSVIS2
       REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3) 
      REAL  RETUW,U1,U2,TEMP
      REAL  U1UP,U1DOWN
      real TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL TAUWE(M1,M3,4)
      REAL A,B,AB1,AB2,AB3,AB4,AB5
	REAL DELTAP,TAUW1,TAUW2
	REAL KAMEN
	  KAPA=0.41

      IFWALL=1.0
	   A=8.3
	   B=0.1428571429
	 AB1=7.202125274  !0.5*(1-B)*A**((1+B)/(1-B))
	 AB2=1.142857143  !1+B
	 AB3=69.74057973  !0.5*A**(2/(1-B))
	 AB4=0.1376936317 !(1+B)/A
	 AB5=1.75         !2/(1+B)                          

      DO J=1,N3M
      KP=KPA(J)
      DO I=1,N1M 
	
        IP=IPA(J)
	 DELTAP=DY(1) !/2
        U1=0.5*(U(I,1,J,1)+U(IP,1,J,1))
        U2=0.5*(U(I,1,J,3)+U(I,1,KP,3))
	  VELOCP=SQRT( U1**2+U2**2)
          TEMP=SQRT(ABS(WSM(1)))/11.81
!         TEMP=1./RE/DELTAP
         VCASE=AB3*TEMP

	  IF    ( VELOCP .LE. VCASE )THEN

	    TAUW1=2.*VELOCP/DELTAP/RE

	  ELSEIF( VELOCP .GT. VCASE )THEN

	    TAUW1=AB1*(1.0/DELTAP/RE)**(1.+B)
	    TAUW1=TAUW1+AB4*(1.0/DELTAP/RE)**B*VELOCP
	    TAUW1=TAUW1**AB5
	  ENDIF

C---------------------------------CENTER ---
	TAUWE(I,J,1)=TAUW1*( U1/VELOCP)
	TAUWE(I,J,2)=TAUW1*( U2/VELOCP )
C---------------------------------

      ENDDO
      ENDDO
      
      
      DO J=1,N3M
      KP=KPA(J)
      DO I=1,N1M 
	
        IP=IPA(J)
	 DELTAP=DY(N2M) !/2
        U1=0.5*(U(I,N2M,J,1)+U(IP,N2M,J,1))
        U2=0.5*(U(I,N2M,J,3)+U(I,N2M,KP,3))
	  VELOCP=SQRT( U1**2+U2**2)
          TEMP=SQRT(ABS(WSM(1)))/11.81
!          TEMP=1./RE/DELTAP
         VCASE=AB3*TEMP

	  IF    ( VELOCP .LE. VCASE )THEN

	    TAUW1=2.*VELOCP/DELTAP/RE

	  ELSEIF( VELOCP .GT. VCASE )THEN

	    TAUW1=AB1*(1.0/DELTAP/RE)**(1.+B)
	    TAUW1=TAUW1+AB4*(1.0/DELTAP/RE)**B*VELOCP
	    TAUW1=TAUW1**AB5
	  ENDIF

C---------------------------------CENTER ---
	TAUWE(I,J,3)=TAUW1*( U1/VELOCP)
	TAUWE(I,J,4)=TAUW1*( U2/VELOCP )
C---------------------------------

      ENDDO
      ENDDO
          
      
      
      
      
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      DO J=1,N3M
      KP=KMA(J)
      DO I=1,N1M 
	
        IP=IMA(J)
        TAUW(I,J,1)=0.5*(TAUWE(I,J,1)+TAUWE(IP,J,1))      
        TAUB(I,J,1)=0.5*(TAUWE(I,J,2)+TAUWE(I,KP,2))
        
        TAUW(I,J,2)=0.5*(TAUWE(I,J,3)+TAUWE(IP,J,3))      
        TAUB(I,J,2)=0.5*(TAUWE(I,J,4)+TAUWE(I,KP,4))     
        
        
       ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!WRITE-BY-FENG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!CCCCCCCCCCCCCCCCCCC!!!!!!!!!!!!!!!!!!!!
        WSM(1)=0.
        WSM(2)=0.
        WSM(3)=0.
         
      DO J=1,N3M
      DO I=1,N1M 
		
  	WSM(1)=WSM(1)+TAUWE(I,J,1)
      ENDDO
      ENDDO
   	WSM(1)=WSM(1)/DBLE(N1M*N3M)

      IS=0
      
      TAOW=WSM(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CCCCCCCCCCCCCCCCCCCCCCC!!!!!!!!!!!!!!!!
      RETURN
      END
      
       SUBROUTINE WALLMODEL_PIOMELLI(U,TAUW,TAUB,SGSVIS)  !BY WANG 0901
      INCLUDE 'PARAM.H'
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
!      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLTAO/TAOW
      COMMON/RANSVISO/RANSVIS,RANSVIS2
       COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH02/ XG(0:M1),YG(0:M2),ZG(0:M3),
     >	           XC(0:M1),YC(0:M2),ZC(0:M3)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)     
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M   
      REAL  U1UP,U1DOWN
      REAL DELTAP,TAUW1,TAUW2
      REAL U_ave,U_tau1,U_tau2,FUtau1,FUtau2
      REAL  RETUW,U1,U2,REINV
      REAL UTAU1,UTAU2
      REAL KAPA,B,TEMP,DTEMP,UMTEMP
      INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UF,VF,WF
   
      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL   TAUWW(0:M1,0:M3,4)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      INTEGER I,K 
      !SOLVED BY NEWTON ITERATIONS
      
      KAPA=0.41
      B=5.2
      IFWALL=1.0
      PI=ACOS(-1.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
      U1DOWN=0.0
      U1UP=0.0
      REINV=1.0/RE
      DO J=1,N3M
      DO I=1,N1M 
      U1DOWN=U1DOWN+U(I,1,J,1)
      U1UP=U1UP+U(I,N2M,J,1)      
      ENDDO
      ENDDO

      U1DOWN=U1DOWN/DBLE(N1M*N3M)
      U1UP=U1UP/DBLE(N1M*N3M)


      DELTAP=DY(1)/2.0             
      U_ave= U1DOWN
      U_tau2=0.2
      U_tau1=0.000001
       B=5.2   
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
       TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
         U_tau2=0.5*U_tau2+0.5*U_tau1
       ELSE
       U_tau1=0.5*U_tau2+0.5*U_tau1
       ENDIF
      ENDDO
      UTAU2=U_tau1
      TAOW=(UTAU2)**2
      WSM(1)= TAOW
      
      TAUW1=(UTAU2)**2
        YPLUS=DELTP*UTAU2*RE
        IF(YPLUS.GT.30.0.AND.YPLUS.LT.50.0)THEN
         DELTS=DELTP/TAN(8.0*PI/180)
        ELSE IF(YPLUS.GE.50)THEN
         DELTS=DELTP/TAN(13.0*PI/180)
        ELSE
         DELTS=0.0
        ENDIF
 
	  YP=H(1)
       DO I=1,N1M
       DO K=1,N3M
           XP=XG(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZC(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
         TAUW(I,K,1)=TAUW1+0.1*UTAU2*VF
           XP=XC(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZG(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
       TAUB(I,K,1)=TAUW1*WF/U1DOWN
      ENDDO
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      DELTAP=DY(N2M)/2.0             
      U_ave= U1UP
      U_tau2=0.2
      U_tau1=0.000001
      DO WHILE(abs(U_tau2-U_tau1).GT.0.0000001)
      TEMP=0.5*U_tau2+0.5*U_tau1
      FUtau1=U_ave/U_tau1-ALOG(U_tau1*DELTAP*RE)/KAPA-B
      FUtau2=U_ave/TEMP-ALOG(TEMP*DELTAP*RE)/KAPA-B
      TEMP=FUtau1*FUtau2
      IF(TEMP.LE.0.0)THEN
         U_tau2=0.5*U_tau2+0.5*U_tau1
       ELSE
       U_tau1=0.5*U_tau2+0.5*U_tau1
       ENDIF
   
      ENDDO

      UTAU2=U_tau1
      
      TAUW1=(UTAU2)**2
     
        YPLUS=DELTP*UTAU1*RE
        IF(YPLUS.GT.30.0.AND.YPLUS.LT.50.0)THEN
         DELTS=DELTP/TAN(8.0*PI/180)
        ELSE IF(YPLUS.GE.50)THEN
         DELTS=DELTP/TAN(13.0*PI/180)
        ELSE
         DELTS=0.0
        ENDIF
 
	  YP=YC(N2M)
       DO I=1,N1M
       DO K=1,N3M
           XP=XG(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZC(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)     
        TAUW(I,K,2)=TAUW1+0.1*UTAU2*VF
           XP=XC(I)+DELTS
          IF(XP.GT.ALX) XP= XP-ALX
      	  ZP=ZG(K)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U) 
       TAUB(I,K,2)=TAUW1*WF/U1UP
      ENDDO
      ENDDO    
      RETURN
      END     
  


      
c      SUBROUTINE NEARWALLTREAD(U,SGSVIS,NTIME)  !BY WANG 0901
c      INCLUDE 'PARAM.H'
c      COMMON/PARA/RE
c      COMMON/MESH2/Y(0:M2)
c      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
c      COMMON/WSMEM/WSM(3),WSMO(3)
c      COMMON/NEARWAL/NEARWALL
c      COMMON/RANSVISO/RANSVIS,RANSVIS2
c      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
c      COMMON/SIZE/ALX,ALY,ALZ,VOL
c      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
c      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
c      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
c      REAL DELTAP,DUDY_TEMP
c      REAL U_tau1,KAPA
c      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
c      REAL XX(4),YY(4)
c      INTEGER I,K ,J,JJ
c      REINV=1.0/RE
c	      KAPA=0.41
c      IF(NEARWALL.EQ.1)THEN 
c       U_tau1=sqrt(WSM(1))
c       DELTAP=DY(1)/2.0  
c        DUDY_TEMP=U_tau1/KAPA/DELTAP 
c        RANSVIS=U_tau1*KAPA*DELTAP     
c      DO K=1,N3M
c      DO I=1,N1M
c        SGSVIS(I,1,K)=AMAX1(RANSVIS 
c     >               +VRMS(1,4)/DUDY_TEMP,-REINV)
c        SGSVIS(I,N2M,K)=SGSVIS(I,N2M,K)*(2.-ALY)
c     >           +(ALY-1.)*AMAX1(RANSVIS-VRMS(N2M,4)/DUDY_TEMP,-REINV)  
c      ENDDO
c      ENDDO
c 
c !         write(*,*)U_tau1,VRMS(1,4),RANSVIS,DUDY_TEMP,KAPA
c       ELSEIF(NEARWALL.EQ.2)THEN
c 
c      DO I=1,N1M
c      DO K=1,N3M
c      DO J=1,4
c      XX(J)=0.5*(Y(J)+Y(J+1))
c      YY(J)=SGSVIS(I,J,K) 
c      ENDDO
c      CALL fit(XX,YY,4,XA,XB)  
c   
c      SGSVIS(I,1,K)=XA+XB*0.5*(Y(1)+Y(2))
c       
c      ENDDO
c      ENDDO
c       DO I=1,N1M
c      DO K=1,N3M
c      DO J=1,4
c          JJ=J+N2M-4
c      XX(J)=0.5*(Y(JJ)+Y(JJ+1))
c      YY(J)=SGSVIS(I,JJ,K) 
c      ENDDO
c      CALL fit(XX,YY,4,XA,XB)  
c   
c      SGSVIS(I,N2M,K)=SGSVIS(I,N2M,K)*(2.-ALY)
c     >               +(ALY-1.)*(XA+XB*0.5*(Y(N2M)+Y(N2)))
c       
c      ENDDO
c      ENDDO
c   
c        
c       ELSE
c       ENDIF
c
c      RETURN
c      END 
C********************** FILTER ***************************
C     THIS SUBROUTINE IS TO FILTER FLOW VARIABLES.
C     TEST FILTER USED IN THIS ROUTINE IS BOX FILTER.
C     ^DX1=2*DX1; ^DX3=2*DX3

      SUBROUTINE FILTER(FILV,J,N)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      REAL FILV(M1,M3,N),TEMP(M1,M3)
   
C     FILTERING IN X1
      DO 5 L=1,N
      DO 10 K=1,N3M
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      FILVM=FILV(IM,K,L)
      FILVC=FILV(I ,K,L)
      FILVP=FILV(IP,K,L)
C     Trapezoidal rule
c      TEMP(I,K)=(FILVM+2.0*FILVC+FILVP)/4.0
C     Simpson's 1/3 rule
      TEMP(I,K)=(FILVM+4.0*FILVC+FILVP)/6.0
   10 CONTINUE

C     FILTERING IN X3
      DO 20 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 20 I=1,N1M
      FILVM=TEMP(I,KM)
      FILVC=TEMP(I,K )
      FILVP=TEMP(I,KP)
C     Trapezoidal rule
c      FILV(I,K,L)=(FILVM+2.0*FILVC+FILVP)*0.25
c      FILV(I,K,L)=(FILVM+2.0*FILVC+FILVP)/4.0
c     Simpson rule
      FILV(I,K,L)=(FILVM+4.0*FILVC+FILVP)/6.0
   20 CONTINUE
    5 CONTINUE

      RETURN
      END
      
CCCC  SPECTRAL CUT-OFF FILTER      
      SUBROUTINE FILTERS(FILV,J,N)      
      INCLUDE 'PARAM.H'
      
      PARAMETER (M1M=M1-1,M3M=M3-1,M3MH=M3M/2+1,M3MD=M3M+2)
      PARAMETER (MM=M1M*M3M)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q      
      
      REAL FILV(M1,M3,N),TEMP(M1,M3)
      
      INTEGER*8 fwd , bwd 
      REAL FRDP(M1,M3MD), FDP(M1,M3MD)
      COMPLEX*16 COEF(M1M,M3M)

      INTEGER fftw_forward
      PARAMETER (fftw_forward=-1)
      INTEGER fftw_backward
      PARAMETER (fftw_backward=+1)

      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)

      INTEGER fftw_estimate
      PARAMETER (fftw_estimate=64)
      
      integer iret
      
      nthds = omp_get_max_threads()
      
      call dfftw_init_threads(iret)      
      call dfftw_plan_with_nthreads(nthds)
      
      call dfftw_plan_dft_2d(fwd, n1m, n3m, COEF, COEF, 
     >                       FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(bwd, n1m, n3m, COEF, COEF, 
     >                      FFTW_BACKWARD, FFTW_ESTIMATE)


      N1MH=N1M/2+1
      N3MH=N3M/2+1      
      
      PI=ACOS(-1.)
      
      DELTA1=2./DX1
      DELTA3=2./DX3
      N1MC=N1M/4+1
      N3MC=N3M/4+1   

      DO 5 L=1,N
       
!$OMP PARALLEL DO      
      DO 1 K=1,N3M
      DO 1 I=1,N1M
      COEF(I,K)=CMPLX(FILV(I,K,L),0.0)         ! MAKE COMPLEX VARIABLE
    1 CONTINUE  
      
      call dfftw_execute_dft(fwd,coef,coef)
      

!     RATIO=DBLE(N1M*N3M) ! DEC ONLY USE
!     REDUCE FOURIER COMPONENT N3M->N3MH
!$OMP PARALLEL DO PRIVATE (KR,KI)
      DO 11 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 11 I=1,N1M
      FRDP(I,KR)=REAL(COEF(I,K))            
      FRDP(I,KI)=AIMAG(COEF(I,K))
   11 CONTINUE

!$OMP PARALLEL DO PRIVATE (KR,KI)      
      DO 110 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 110 I=1,N1M
      IF (I.GE.N1MC.AND.I.LE.(N1M-N1MC+2).AND.K.GE.N3MC) THEN
      FRDP(I,KR)=0.0
      FRDP(I,KI)=0.0
      ENDIF
  110 CONTINUE  
      
!     DO THE INVERSE FFT.
!$OMP PARALLEL DO PRIVATE (KR,KI)      
      DO 21 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 21 I=1,N1M
      COEF(I,K)=CMPLX(FRDP(I,KR),FRDP(I,KI))
   21 CONTINUE

!$OMP PARALLEL DO PRIVATE (KK,II)      
      DO 22 K=N3MH+1,N3M
      KK=N3M-K+2
      DO 22 I=1,N1M
      II=N1M-I+2
      IF(I.EQ.1) II=1
      COEF(I,K)=CONJG(COEF(II,KK))
   22 CONTINUE

      call dfftw_execute_dft(bwd, COEF, COEF)

!$OMP PARALLEL DO       
      DO 23 K=1,N3M
      DO 23 I=1,N1M
      FILV(I,K,L)=REAL(COEF(I,K))/FLOAT(N1M*N3M)
   23 CONTINUE
            
    5 CONTINUE
      
      call dfftw_destroy_plan(FWD)
      call dfftw_destroy_plan(BWD)

      RETURN
      END      
      
CCCC  2ND SPECTRAL CUT-OFF FILTER      
      SUBROUTINE FILTER4S(FILV,J,N)      
      INCLUDE 'PARAM.H'
      
      PARAMETER (M1M=M1-1,M3M=M3-1,M3MH=M3M/2+1,M3MD=M3M+2)
      PARAMETER (MM=M1M*M3M)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q      
      
      REAL FILV(M1,M3,N),TEMP(M1,M3)
      
      INTEGER*8 fwd , bwd 
      REAL FRDP(M1,M3MD), FDP(M1,M3MD)
      COMPLEX*16 COEF(M1M,M3M)

      INTEGER fftw_forward
      PARAMETER (fftw_forward=-1)
      INTEGER fftw_backward
      PARAMETER (fftw_backward=+1)

      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)

      INTEGER fftw_estimate
      PARAMETER (fftw_estimate=64)
      
      integer iret
      
      nthds = omp_get_max_threads()
      
      call dfftw_init_threads(iret)      
      call dfftw_plan_with_nthreads(nthds)
      
      call dfftw_plan_dft_2d(fwd, n1m, n3m, COEF, COEF, 
     >                       FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(bwd, n1m, n3m, COEF, COEF, 
     >                      FFTW_BACKWARD, FFTW_ESTIMATE)


      N1MH=N1M/2+1
      N3MH=N3M/2+1      
      
      PI=ACOS(-1.)
      
      DELTA1=4./DX1
      DELTA3=4./DX3
      N1MC=N1M/8+1
      N3MC=N3M/8+1   

      DO 5 L=1,N
       
!$OMP PARALLEL DO      
      DO 1 K=1,N3M
      DO 1 I=1,N1M
      COEF(I,K)=CMPLX(FILV(I,K,L),0.0)         ! MAKE COMPLEX VARIABLE
    1 CONTINUE  
      
      call dfftw_execute_dft(fwd,coef,coef)
      

!     RATIO=DBLE(N1M*N3M) ! DEC ONLY USE
!     REDUCE FOURIER COMPONENT N3M->N3MH
!$OMP PARALLEL DO PRIVATE (KR,KI)
      DO 11 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 11 I=1,N1M
      FRDP(I,KR)=REAL(COEF(I,K))            
      FRDP(I,KI)=AIMAG(COEF(I,K))
   11 CONTINUE

!$OMP PARALLEL DO PRIVATE (KR,KI)      
      DO 110 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 110 I=1,N1M
      IF (I.GE.N1MC.AND.I.LE.(N1M-N1MC+2).AND.K.GE.N3MC) THEN
      FRDP(I,KR)=0.0
      FRDP(I,KI)=0.0
      ENDIF
  110 CONTINUE  
      
!     DO THE INVERSE FFT.
!$OMP PARALLEL DO PRIVATE (KR,KI)      
      DO 21 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 21 I=1,N1M
      COEF(I,K)=CMPLX(FRDP(I,KR),FRDP(I,KI))
   21 CONTINUE

!$OMP PARALLEL DO PRIVATE (KK,II)      
      DO 22 K=N3MH+1,N3M
      KK=N3M-K+2
      DO 22 I=1,N1M
      II=N1M-I+2
      IF(I.EQ.1) II=1
      COEF(I,K)=CONJG(COEF(II,KK))
   22 CONTINUE

      call dfftw_execute_dft(bwd, COEF, COEF)

!$OMP PARALLEL DO       
      DO 23 K=1,N3M
      DO 23 I=1,N1M
      FILV(I,K,L)=REAL(COEF(I,K))/FLOAT(N1M*N3M)
   23 CONTINUE
            
    5 CONTINUE
      
      call dfftw_destroy_plan(FWD)
      call dfftw_destroy_plan(BWD)

      RETURN
      END         
      
C********************** FILTER4 ***************************
C     THIS SUBROUTINE IS TO FILTER FLOW VARIABLES.
C     TEST FILTER USED IN THIS ROUTINE IS BOX FILTER.
C     ^DX1=4*DX1; ^DX3=4*DX3

      SUBROUTINE FILTER4(FILV,J,N)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      REAL FILV(M1,M3,N),TEMP(M1,M3)
      DIMENSION IPA2(N1M),IMA2(N1M),KPA2(N3M),KMA2(N3M)
      
      DO IC=1,N1M
      IPA2(IC)=IC+2
      IMA2(IC)=IC-2
      ENDDO
      IPA2(N1M)=2
      IPA2(N1M-1)=1
      IMA2(1)=N1M-1
      IMA2(2)=N1M
      
      DO KC=1,N3M
      KPA2(KC)=KC+2
      KMA2(KC)=KC-2
      ENDDO
      KPA2(N3M)=2
      KPA2(N3M-1)=1
      KMA2(1)=N3M-1
      KMA2(2)=N3M
      
C     FILTERING IN X1
      DO 5 L=1,N
      DO 10 K=1,N3M
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      IP2=IPA2(I)
      IM2=IMA2(I)
      FILVM2=FILV(IM2,K,L)
      FILVM=FILV(IM,K,L)
      FILVC=FILV(I ,K,L)
      FILVP=FILV(IP,K,L)
      FILVP2=FILV(IP2,K,L)
C     Trapezoidal rule
c      TEMP(I,K)=(FILVM2+2.0*FILVM+2.0*FILVC+2.0*FILVP+FILVP2)/8.0
c     Simpson rule
      TEMP(I,K)=(FILVM2+4.0*FILVM+2.0*FILVC+4.0*FILVP+FILVP2)/12.0
   10 CONTINUE

C     FILTERING IN X3
      DO 20 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      KP2=KPA2(K)
      KM2=KMA2(K)
      DO 20 I=1,N1M
      FILVM2=TEMP(I,KM2)
      FILVM=TEMP(I,KM)
      FILVC=TEMP(I,K )
      FILVP=TEMP(I,KP)
      FILVP2=TEMP(I,KP2)
C     Simpson's 1/3 rule
c      FILV(I,K,L)=(FILVM2+2.0*FILVM+2.0*FILVC+2.0*FILVP+FILVP2)/8.0
c     Simpson rule
      FILV(I,K,L)=(FILVM2+4.0*FILVM+2.0*FILVC+4.0*FILVP+FILVP2)/12.0
   20 CONTINUE
      
    5 CONTINUE

      RETURN
      END      

      
C************************ STRAIN ****************************
C     THIS SUBROUTINE IS TO CALCULATE THE STRAIN-RATE TENSOR.
C     SR(I,J,K,1) = S_{11}
C     SR(I,J,K,2) = S_{12}
C     SR(I,J,K,3) = S_{13}
C     SR(I,J,K,4) = S_{22}
C     SR(I,J,K,5) = S_{23}
C     SR(I,J,K,6) = S_{33}

      SUBROUTINE STRAIN(U,SR,J)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS

      REAL U(0:M1,0:M2,0:M3,3),SR(M1,M3,6)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL WMODEL,IFWALL
      INTEGER JP,JM,KP,KM,IP,IM,FJUM,FJUP
      REAL U2,U1,V2,V1,W2,W1
      REAL UJP,UJC,UJM,VIC,VIP,VIM
      REAL UKP,UKC,UKM,WIP,WIC,WIM
      REAL WJP,WJC,WJM,VKP,VKC,VKM
      
C     THE STRAIN-RATES ARE EVALUATED AT THE CELL CENTER.
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      
      U2=U(IP,J,K,1)      
      U1=U(I ,J,K,1)      
      SR(I,K,1)=DX1*(U2-U1)
      
      UJP=0.5*(U(I,JP,K,1)+U(IP,JP,K,1))
      UJC=0.5*(U(I,J ,K,1)+U(IP,J ,K,1)) 
      UJM=0.5*(U(I,JM,K,1)+U(IP,JM,K,1)) 
      IF (WMODEL.EQ.0.) THEN
      U2=FJUP*1.0/H(JP)*(DY(J)/2.0*UJP+DY(JP)/2.0*UJC)
     >   +(1.-FJUP)*0.5*(U(I,N2,K,1)+U(IP,N2,K,1))
      U1=FJUM*1.0/H(J )*(DY(J)/2.0*UJM+DY(JM)/2.0*UJC)
     >   +(1.-FJUM)*0.5*(U(I,0,K,1)+U(IP,0,K,1))
      ELSE
      U2=FJUP*1.0/H(JP)*(DY(J)/2.0*UJP+DY(JP)/2.0*UJC)
     >   +(1.-FJUP)*UJC
      U1=FJUM*1.0/H(J )*(DY(J)/2.0*UJM+DY(JM)/2.0*UJC)
     >   +(1.-FJUM)*UJC
      ENDIF          
      VIP=0.5*(U(IP,J,K,2)+U(IP,JP,K,2))
      VIC=0.5*(U(I ,J,K,2)+U(I ,JP,K,2)) 
      VIM=0.5*(U(IM,J,K,2)+U(IM,JP,K,2)) 
      V2=0.5*(VIP+VIC)
      V1=0.5*(VIM+VIC)
      IF (WMODEL.EQ.0.) THEN
      SR(I,K,2)=0.5*(DX1*(V2-V1)+1.0/DY(J)*(U2-U1))
      ENDIF
      IF (WMODEL.NE.0.) THEN
      IF (J.EQ.1) THEN
      SR(I,K,2)=0.5*(DX1*(V2-V1)+2.0/DY(J)*(U2-U1))
      ELSEIF (J.EQ.N2M) THEN
      SR(I,K,2)=0.5*(DX1*(V2-V1)+2.0/DY(J)*(U2-U1))
      ELSE
      SR(I,K,2)=0.5*(DX1*(V2-V1)+1.0/DY(J)*(U2-U1))
      ENDIF
      ENDIF
!      SR(I,K,2)=0.5*DX1*(V2-V1)+0.5*(1.0/DY(J)*(U2-U1))*FJUM*FJUP
!     >       +0.5*(1.0/DY(J)*(U2-U1)
!     >          +0.09/DY(J)*VM(1,1))*(1-FJUM)
!     >       +0.5*(1.0/DY(J)*(U2-U1)
!     >          -0.09/DY(J)*VM(N2M,1))*(1-FJUP) !Porte-Agel et al (2000)
      
      UKP=0.5*(U(I,J,KP,1)+U(IP,J,KP,1))
      UKC=0.5*(U(I,J,K ,1)+U(IP,J,K ,1)) 
      UKM=0.5*(U(I,J,KM,1)+U(IP,J,KM,1)) 
      U2=0.5*(UKP+UKC)
      U1=0.5*(UKM+UKC)
      WIP=0.5*(U(IP,J,K,3)+U(IP,J,KP,3))
      WIC=0.5*(U(I ,J,K,3)+U(I ,J,KP,3)) 
      WIM=0.5*(U(IM,J,K,3)+U(IM,J,KP,3)) 
      W2=0.5*(WIP+WIC)
      W1=0.5*(WIM+WIC)
      SR(I,K,3)=0.5*(DX1*(W2-W1)+DX3*(U2-U1))
      
      V2=U(I,JP,K,2)      
      V1=U(I,J ,K,2)      
      SR(I,K,4)=1.0/DY(J)*(V2-V1)
      
      WJP=0.5*(U(I,JP,K,3)+U(I,JP,KP,3))
      WJC=0.5*(U(I,J ,K,3)+U(I,J ,KP,3)) 
      WJM=0.5*(U(I,JM,K,3)+U(I,JM,KP,3)) 
      IF (WMODEL.EQ.0.) THEN
      W2=FJUP*1.0/H(JP)*(DY(J)/2.0*WJP+DY(JP)/2.0*WJC)
     >   +(1.-FJUP)*0.5*(U(I,N2,K,3)+U(I,N2,KP,3))
      W1=FJUM*1.0/H(J )*(DY(J)/2.0*WJM+DY(JM)/2.0*WJC)
     >   +(1.-FJUM)*0.5*(U(I,0,K,3)+U(I,0,KP,3))
      ELSE
      W2=FJUP*1.0/H(JP)*(DY(J)/2.0*WJP+DY(JP)/2.0*WJC)
     >   +(1.-FJUP)*WJC
      W1=FJUM*1.0/H(J )*(DY(J)/2.0*WJM+DY(JM)/2.0*WJC)
     >   +(1.-FJUM)*WJC
      ENDIF
      VKP=0.5*(U(I,J,KP,2)+U(I,JP,KP,2))
      VKC=0.5*(U(I,J,K ,2)+U(I,JP,K ,2)) 
      VKM=0.5*(U(I,J,KM,2)+U(I,JP,KM,2)) 
      V2=0.5*(VKP+VKC)
      V1=0.5*(VKM+VKC)
      IF (WMODEL.EQ.0.) THEN
      SR(I,K,5)=0.5*(1.0/DY(J)*(W2-W1)+DX3*(V2-V1))
      ENDIF
      IF (WMODEL.NE.0.) THEN
      IF (J.EQ.1) THEN
      SR(I,K,5)=0.5*(2.0/DY(J)*(W2-W1)+DX3*(V2-V1))
      ELSEIF (J.EQ.N2M) THEN
      SR(I,K,5)=0.5*(2.0/DY(J)*(W2-W1)+DX3*(V2-V1))
      ELSE
      SR(I,K,5)=0.5*(1.0/DY(J)*(W2-W1)+DX3*(V2-V1))
      ENDIF
      ENDIF
      
      W2=U(I,J,KP,3)      
      W1=U(I,J,K ,3)      
      SR(I,K,6)=DX3*(W2-W1)
      
   10 CONTINUE
   
      RETURN
      END
C************************ STRAIN1 ****************************
C     THIS SUBROUTINE IS TO CALCULATE THE STRAIN-RATE TENSOR.
C     SR(I,J,K,1) = S_{11}
C     SR(I,J,K,2) = S_{12}
C     SR(I,J,K,3) = S_{13}
C     SR(I,J,K,4) = S_{22}
C     SR(I,J,K,5) = S_{23}
C     SR(I,J,K,6) = S_{33}

      SUBROUTINE STRAIN1(U,SR,J)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS

      REAL U(0:M1,0:M2,0:M3,3),SR(M1,M3,6)
      REAL TAUW(0:M1,0:M3,2)
      REAL TAUB(0:M1,0:M3,2)
      REAL WMODEL,IFWALL
      INTEGER JP,JM,KP,KM,IP,IM,FJUM,FJUP
      REAL U2,U1,V2,V1,W2,W1
      REAL UJP,UJC,UJM,VIC,VIP,VIM
      REAL UKP,UKC,UKM,WIP,WIC,WIM
      REAL WJP,WJC,WJM,VKP,VKC,VKM
      REAL A11(M1-1),B11(M1-1),C11(M1-1),D11(M1-1),X11(M1-1),X13(M1-1)
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1),X1(M1-1),X12(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1),X3(M3-1)
      REAL A33(M3-1),B33(M3-1),C33(M3-1),D33(M3-1),X33(M3-1)
      REAL UK(M1-1,M3-1),WI(M1-1,M3-1),X22(M3-1),X23(M3-1)

C     THE STRAIN-RATES ARE EVALUATED AT THE CELL CENTER.
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)

!$omp parallel do
      DO I=1,N1M
          A11(I)=1.0
          B11(I)=22.0
          C11(I)=1.0
      ENDDO 

!$omp parallel do
      DO I=1,N1M
          A1(I)=1.0
          B1(I)=6.0
          C1(I)=1.0
      ENDDO

!$omp parallel do
      DO K=1,N3M
          A3(K)=1.0
          B3(K)=6.0
          C3(K)=1.0
      ENDDO

!$omp parallel do
      DO K=1,N3M
          A33(K)=1.0
          B33(K)=22.0
          C33(K)=1.0
      ENDDO

CCC      
!$omp parallel do private(IP,IM,D11,X11)
      DO  K=1,N3M

      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      D11(I)=24.*DX1*(U(IP,J,K,1)-U(I,J,K,1))
      ENDDO

      CALL CTDMA11(A11,B11,C11,D11,X11)

      DO I=1,N1M
          SR(I,K,1)=X11(I)
      ENDDO

      ENDDO

CCC   
!$omp parallel do private(IP,IM,D1,D11,X1,X11,X12,X13,
!$omp&  UJP,UJC,UJM,U2,U1,V2,V1)
      DO K=1,N3M 

      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      D1(I)=4.*(U(I,JP,K,1)+U(IP,JP,K,1))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      D1(I)=4.*(U(I,J,K,1)+U(IP,J,K,1))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X12)

      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      D1(I)=4.*(U(I,JM,K,1)+U(IP,JM,K,1))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X13)

      DO I=1,N1M    
      UJP=X1(I)
      UJC=X12(I)
      UJM=X13(I)
      IF (WMODEL.EQ.0.) THEN
      U2=FJUP*1.0/H(JP)*(DY(J)/2.0*UJP+DY(JP)/2.0*UJC)
     >   +(1.-FJUP)*0.5*(U(I,N2,K,1)+U(IP,N2,K,1))
      U1=FJUM*1.0/H(J )*(DY(J)/2.0*UJM+DY(JM)/2.0*UJC)
     >   +(1.-FJUM)*0.5*(U(I,0,K,1)+U(IP,0,K,1))
      ELSE
      U2=FJUP*1.0/H(JP)*(DY(J)/2.0*UJP+DY(JP)/2.0*UJC)
     >   +(1.-FJUP)*UJC
      U1=FJUM*1.0/H(J )*(DY(J)/2.0*UJM+DY(JM)/2.0*UJC)
     >   +(1.-FJUM)*UJC
      ENDIF 
      IF (WMODEL.EQ.0.) THEN
      SR(I,K,2)=0.5*1.0/DY(J)*(U2-U1)
      ENDIF
      IF (WMODEL.NE.0.) THEN
      IF (J.EQ.1.) THEN
      SR(I,K,2)=0.5*2.0/DY(J)*(U2-U1)
      ELSE IF (J.EQ.N2M) THEN
      SR(I,K,2)=0.5*2.0/DY(J)*(U2-U1) 
      ELSE
      SR(I,K,2)=0.5*1.0/DY(J)*(U2-U1)
      ENDIF
      ENDIF
      ENDDO
      
      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)    
      D1(I)=4.*(U(IM,J,K,2)+U(I,J,K,2))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X1)
      
      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)    
      D1(I)=4.*(U(IM,JP,K,2)+U(I,JP,K,2))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X12)
      
      DO I=1,N1M
      IP=IPA(I)
      V2=X12(IP)+X1(IP)
      V1=X12(I)+X1(I)
      D11(I)=24.*DX1*(V2-V1)
      ENDDO

      CALL CTDMA11(A11,B11,C11,D11,X11)

      DO I=1,N1M
      SR(I,K,2)=SR(I,K,2)+X11(I)/2.
      ENDDO

      ENDDO

CCC

!$omp parallel do private(KP,KM,IP,IM,D1,X1)   
      DO  K=1,N3M
      KP=KPA(K)
      KM=KMA(K) 
    
      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      D1(I)=4.*(U(IP,J,K,1)+U(I,J,K,1))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      UK(I,K)=X1(I)
      ENDDO

      ENDDO

!$omp parallel do private(IP,IM,KP,KM,D3,X3,D33,X33)
      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)   

      DO  K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      D3(K)=4.*(UK(I,KM)+UK(I,K))
      ENDDO

      CALL CTDMA33(A3,B3,C3,D3,X3)
         
      DO  K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      D33(K)=24.*DX3*(X3(KP)-X3(K))
      ENDDO

      CALL CTDMA33(A33,B33,C33,D33,X33)

      DO K=1,N3M
      SR(I,K,3)=X33(K)/2.
      ENDDO
      
      DO K=1,N3M
      KP=KPA(K)
      D3(K)=4.*(U(I,J,KP,3)+U(I,J,K,3))
      ENDDO

      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      WI(I,K)=X3(K)
      ENDDO

      ENDDO
      
!$omp parallel do private(IM,IP,D1,X1,D11,X11)
      DO K=1,N3M

      DO I=1,N1M
      IM=IMA(I)
      IP=IPA(I)
      D1(I)=4.*(WI(IM,K)+WI(I,K))
      ENDDO

      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      IM=IMA(I)
      IP=IPA(I)
      D11(I)=24.*(X1(IP)-X1(I))
      ENDDO

      CALL CTDMA11(A11,B11,C11,D11,X11)

      DO I=1,N1M
      SR(I,K,3)=SR(I,K,3)+X11(K)/2.
      ENDDO

      ENDDO

ccc      
!$omp parallel do private(KP,KM,IP,IM,V2,V1)
      DO  K=1,N3M
      KP=KPA(K)
      KM=KMA(K)

      DO  I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      V2=U(I,JP,K,2)      
      V1=U(I,J ,K,2)      
      SR(I,K,4)=1.0/DY(J)*(V2-V1)
      ENDDO

      ENDDO
CCC
!$omp parallel do private(KP,KM,D3,D33,X3,X22,X23,X33,
!$omp& WJP,WJC,WJM,W2,W1,VKC,VKM)
      DO I=1,N1M   
  
      DO K=1,N3M
      KP=KPA(K)
      D3(K)=4.*(U(I,JP,K,3)+U(I,JP,KP,3))
      ENDDO 

      CALL CTDMA33(A3,B3,C3,D3,X3)
      
      DO K=1,N3M
      KP=KPA(K)
      D3(K)=4.*(U(I,J,K,3)+U(I,J,KP,3))
      ENDDO
 
      CALL CTDMA33(A3,B3,C3,D3,X22)

      DO K=1,N3M
      KP=KPA(K)
      D3(K)=4.*(U(I,JM,K,3)+U(I,JM,KP,3))
      ENDDO 

      CALL CTDMA33(A3,B3,C3,D3,X23)

      DO K=1,N3M
      WJP=X3(K)
      WJC=X22(K)
      WJM=X23(K)
      IF (WMODEL.EQ.0.) THEN
      W2=FJUP*1.0/H(JP)*(DY(J)/2.0*WJP+DY(JP)/2.0*WJC)
     >   +(1.-FJUP)*0.5*(U(I,N2,K,3)+U(I,N2,KP,3))
      W1=FJUM*1.0/H(J )*(DY(J)/2.0*WJM+DY(JM)/2.0*WJC)
     >   +(1.-FJUM)*0.5*(U(I,0,K,3)+U(I,0,KP,3))
      ELSE
      W2=FJUP*1.0/H(JP)*(DY(J)/2.0*WJP+DY(JP)/2.0*WJC)
     >   +(1.-FJUP)*WJC
      W1=FJUM*1.0/H(J )*(DY(J)/2.0*WJM+DY(JM)/2.0*WJC)
     >   +(1.-FJUM)*WJC
      ENDIF
      IF (WMODEL.EQ.0.) THEN
      SR(I,K,5)=0.5*1.0/DY(J)*(W2-W1)
      ENDIF
      IF (WMODEL.NE.0.) THEN
      IF (J.EQ.1) THEN
      SR(I,K,5)=0.5*2.0/DY(J)*(W2-W1)
      ELSE IF (J.EQ.N2M) THEN
      SR(I,K,5)=0.5*2.0/DY(J)*(W2-W1)
      ELSE
      SR(I,K,5)=0.5*1.0/DY(J)*(W2-W1)   
      ENDIF
      ENDIF
      ENDDO

      DO K=1,N3M
      KM=KMA(K)
      VKC=0.5*(U(I,J,K ,2)+U(I,JP,K ,2)) 
      VKM=0.5*(U(I,J,KM,2)+U(I,JP,KM,2))
      D3(K)=4.*(VKC+VKM)
      ENDDO

      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      KP=KPA(K)
      D33(K)=24.*DX3*(X3(KP)-X3(K))
      ENDDO

      CALL CTDMA33(A33,B33,C33,D33,X33)

      DO K=1,N3M
      SR(I,K,5)=SR(I,K,5)+X33(K)/2.
      ENDDO

      ENDDO

CCC      
!$omp parallel do private(KP,D33,X33)
      DO I=1,N1M

      DO K=1,N3M
      KP=KPA(K)
      D33(K)=24.*DX3*(U(I,J,KP,3)-U(I,J,K,3))
      ENDDO

      CALL CTDMA33(A33,B33,C33,D33,X33)

      DO K=1,N3M     
      SR(I,K,6)=X33(K)
      ENDDO

      ENDDO
   
      RETURN
      END

c***************** RHS1 ***********************  
      SUBROUTINE RHS1(U,P,HH,HH2,SS,SS2,RUH1,PRESG,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC(0:M1,0:M3,3,2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/FINDX3/FUP(M2),FDOWN(M2)

      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL HH2(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL SS2(0:M1,0:M2,0:M3,3)
      REAL RUH1(0:M1,0:M2,0:M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL NUTI,NUTJ,NUTK
      
      REAL X1(M1-1),X11(M1-1),X3(M3-1),X33(M3-1)
      REAL UU1(M1-1),WW1(M1-1,M3-1),UU3(M1-1,M3-1),VV1(M1-1,M2-1)
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1)
      REAL A11(M1-1),B11(M1-1),C11(M1-1),D11(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1)
      REAL A33(M3-1),B33(M3-1),C33(M3-1),D33(M3-1)
      REAL SGSVIS2(M1-1,M3-1)
      REAL SGSDUY,SGSVIS1(M1-1,M2-1)
      REAL SGSDUXY,SGSVIS12(M1-1,M2-1),UU12(M1-1,M2-1)
      REAL SGSVIS13(M1-1,M3-1),WW13(M1-1,M3-1)      

!$omp  parallel do private(IP,IM,JP,JM,KP,KM,FJUM,FJUP,
!$omp& VISCOS,PRESSG1,BC_DOWN,BC_UP,BC)
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

C     VISCOUS TERM: (0.5/RE)*D2U/DY2 
      VISCOS=
     >      (1./RE)*(HP(J)*U(I,JP,K,1)
     >                -HC(J)*U(I,J ,K,1)
     >     +HM(J)*U(I,JM,K,1))
      RUH1(I,J,K)=VISCOS

      
      
C      RUH1(I,J,K)= RUH1(I,J,K)+BC
      
 10   CONTINUE
      
 

CCC   VISCOS TERM: (1/RE)*D2U/DX2
!$OMP PARALLEL DO
      DO 15 I=1,N1M
      A1(I)=1.0
      B1(I)=10.0
      C1(I)=1.0
 15   CONTINUE      

!$OMP PARALLEL DO PRIVATE (D1,X1)
      DO 20 K=1,N3M
      DO 20 J=1,N2M    

      DO I=2,N1M-1
      D1(I)=12.*DX1Q*(U(I-1,J,K,1)-2.*U(I,J,K,1)+U(I+1,J,K,1))
      ENDDO
      D1(1)=12.*DX1Q*(U(N1M,J,K,1)-2.*U(1,J,K,1)+U(2,J,K,1))
      D1(N1M)=12.*DX1Q*(U(N1M-1,J,K,1)-2.*U(N1M,J,K,1)+U(1,J,K,1))

      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      RUH1(I,J,K)=RUH1(I,J,K)+(1./RE)*X1(I)
      ENDDO      
      
  20  CONTINUE
      
      
CCC   VISCOS TERM: (1/RE+NUTI)*D2U/DZ2      
!$OMP PARALLEL DO
      DO 25 K=1,N3M
      A3(K)=1.0
      B3(K)=10.0
      C3(K)=1.0
 25   CONTINUE     

      
!$OMP PARALLEL DO PRIVATE (D3,X3)
      DO 30 I=1,N1M
      DO 30 J=1,N2M     
      
      DO K=2,N3M-1
      D3(K)=12.*DX3Q*(U(I,J,K-1,1)-2.*U(I,J,K,1)+U(I,J,K+1,1))
      ENDDO
      D3(1)=12.*DX3Q*(U(I,J,N3M,1)-2.*U(I,J,1,1)+U(I,J,2,1))
      D3(N3M)=12.*DX3Q*(U(I,J,N3M-1,1)-2.*U(I,J,N3M,1)+U(I,J,1,1))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      RUH1(I,J,K)=RUH1(I,J,K)+(1./RE)*X3(K)   
      ENDDO      
      
  30  CONTINUE
      
CCC   PRESSURE GRADIENT TERM
!$OMP PARALLEL DO
      DO 35 I=1,N1M
      A1(I)=1.0
      B1(I)=22.0
      C1(I)=1.0
 35   CONTINUE      
      
!$OMP PARALLEL DO PRIVATE (D1,X1)
      DO 40 K=1,N3M
      DO 40 J=1,N2M     
    
      DO I=2,N1M
      D1(I)=24.*DX1*(-P(I-1,J,K)+P(I,J,K))
      ENDDO      
      D1(1)=24.*DX1*(-P(N1M,J,K)+P(1,J,K))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      RUH1(I,J,K)=RUH1(I,J,K)-X1(I)-PRESG 
      ENDDO      
      
  40  CONTINUE   
      
      
      
CCC   CONVECTION TERM HH 
      
C     D(UV)/DY
!$OMP PARALLEL DO
      DO I=1,N1M
      A1(I)=1.0
      B1(I)=6.0
      C1(I)=1.0
      ENDDO 

      DO K=1,N3M
      KP=KPA(K)
      KM=KMA(K)    
     
!$omp parallel do private(D1,X1)       
      DO J=1,N2M
      
      DO I=2,N1M
      D1(I) = 4.0*(U(I-1,J,K,2)+U(I,J,K,2))
      ENDDO        
      D1(1) = 4.0*(U(N1M,J,K,2)+U(1,J,K,2))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)
      
      DO I=1,N1M
      VV1(I,J)=X1(I)
      ENDDO       
      
      ENDDO

!$omp parallel do private(JP,JM,FJUM,FJUP,IP,IM,U2,U1,V2,V1)      
      DO J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

C     D(UV)/DY
      U2 = FJUP*1.0/H(JP)*(DY(J )/2.*U(I,JP,K,1)+DY(JP)/2.*U(I,J ,K,1))
     >  +(1.-FJUP)*U(I,N2,K,1)
      U1 = FJUM*1.0/H(J )*(DY(JM)/2.*U(I,J ,K,1)+DY(J )/2.*U(I,JM,K,1))
     >  +(1.-FJUM)*U(I,0,K,1)
      
      V2 = VV1(I,JP)
      V1 = VV1(I,J)
      
      IF (J.EQ.N2M) V2=0.0
      IF (J.EQ.1)   V1=0.0
      
      HH(I,J,K,1)=HH(I,J,K,1)+1./DY(J)*(V2*U2-V1*U1)  
      
C     (1/2)V(DU/DY)
C      V3 = 0.5*(V2+V1)
C      HH(I,J,K,1)=HH(I,J,K,1)+V3*1./DY(J)*(U2-U1)/2.
      
      ENDDO
      ENDDO
      
      ENDDO
      
      
C     D(UU)/DX 
!$OMP PARALLEL DO
      DO I=1,N1M
      A1(I)=1.0
      B1(I)=6.0
      C1(I)=1.0
      ENDDO 
      
!$OMP PARALLEL DO
      DO I=1,N1M
      A11(I)=1.0
      B11(I)=22.0
      C11(I)=1.0
      ENDDO      
      
!$OMP PARALLEL DO PRIVATE(D1,X1,UU1,D11,X11)
      DO 50 K=1,N3M
      DO 50 J=1,N2M     
      
      ! INTERPOLATION OF U IN X DIRECTION     
      DO I=1,N1M-1
      D1(I)  =4.0*(U(I,J,K,1)+U(I+1,J,K,1))
      ENDDO        
      D1(N1M)=4.0*(U(N1M,J,K,1)+U(1,J,K,1))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      UU1(I)=X1(I)  
      ENDDO        
      
C     D(UU)/DX 
      DO I=2,N1M
      D11(I)=24.*DX1*(UU1(I)**2-UU1(I-1)**2)
      ENDDO      
      D11(1)=24.*DX1*(UU1(1)**2-UU1(N1M)**2)
      
      CALL CTDMA11(A11,B11,C11,D11,X11)

      DO I=1,N1M
      HH(I,J,K,1)=HH(I,J,K,1)+X11(I)
      ENDDO  
      
C     (1/2)U(DU/DX)
C      DO I=2,N1M
C      D11(I)=24.*DX1*(UU1(I)-UU1(I-1))
C      ENDDO      
C      D11(1)=24.*DX1*(UU1(1)-UU1(N1M))
      
C      CALL CTDMA11(A11,B11,C11,D11,X11)
      
C      DO I=1,N1M
C      HH(I,J,K,1)=HH(I,J,K,1)+U(I,J,K,1)*X11(I)/2.
C      ENDDO 
      
  50  CONTINUE   
      

C     D(UW)/DZ
!$OMP PARALLEL DO
      DO K=1,N3M
      A3(K)=1.0
      B3(K)=6.0
      C3(K)=1.0
      ENDDO  
      
!$OMP PARALLEL DO
      DO I=1,N1M
      A1(I)=1.0
      B1(I)=6.0
      C1(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
      DO K=1,N3M
      A33(K)=1.0
      B33(K)=22.0
      C33(K)=1.0
      ENDDO
      

      DO 60 J=1,N2M      
      
!$OMP PARALLEL DO PRIVATE (D3,X3)      
      DO 55 I=1,N1M      
      
      ! INTERPOLATION OF U IN Z DIRECTION     
      DO K=2,N3M
      D3(K)=4.*(U(I,J,K-1,1)+U(I,J,K,1))
      ENDDO      
      D3(1)=4.*(U(I,J,N3M,1)+U(I,J,1,1))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      UU3(I,K)=X3(K)      
      ENDDO 
      
  55  CONTINUE
      
      ! INTERPOLATION OF W IN X DIRECTION
!$OMP PARALLEL DO PRIVATE (D1,X1)      
      DO 57 K=1,N3M      

      DO I=2,N1M
      D1(I)=4.*(U(I-1,J,K,3)+U(I,J,K,3))
      ENDDO      
      D1(1)=4.*(U(N1M,J,K,3)+U(1,J,K,3))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      WW1(I,K)=X1(I)      
      ENDDO 
      
  57  CONTINUE
      
C     (DUW/DZ)
!$OMP PARALLEL DO PRIVATE (D33,X33)      
      DO 59 I=1,N1M

      DO K=1,N3M-1
      D33(K)=24.*DX3*(-UU3(I,K)*WW1(I,K)+UU3(I,K+1)*WW1(I,K+1))
      ENDDO      
      D33(N3M)=24.*DX3*(-UU3(I,N3M)*WW1(I,N3M)+UU3(I,1)*WW1(I,1))
      
      CALL CTDMA33(A33,B33,C33,D33,X33)      

      DO K=1,N3M
      HH(I,J,K,1)=HH(I,J,K,1)+X33(K)  
c      WRITE(*,*) HH(I,J,K,1)
      ENDDO
      
C     (1/2)W(DU/DZ)
C      DO K=1,N3M-1
C      D3(K)=4.*(WW1(I,K)+WW1(I,K+1))
C      ENDDO      
C      D3(N3M)=4.*(WW1(I,N3M)+WW1(I,1))
      
C      CALL CTDMA33(A3,B3,C3,D3,X3)
      
C      DO K=1,N3M-1
C      D33(K)=24.*DX3*(-UU3(I,K)+UU3(I,K+1))
C      ENDDO      
C      D33(N3M)=24.*DX3*(-UU3(I,N3M)+UU3(I,1))
      
C      CALL CTDMA33(A33,B33,C33,D33,X33) 
      
C      DO K=1,N3M
C      HH(I,J,K,1)=HH(I,J,K,1)+X3(K)*X33(K)/2.
C      ENDDO
     
  59  CONTINUE
      
  60  CONTINUE 
      
      
C-----------------D(SGS*(DU/DX))/DX---------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
	ENDDO
      
!$OMP PARALLEL DO PRIVATE (D11,X11)
      DO 61 K=1,N3M
	DO 61 J=1,N2M
	
      DO I=1,N1M-1
      D11(I)=24.*DX1*(U(I+1,J,K,1)-U(I,J,K,1))
	ENDDO
	D11(N1M)=24.*DX1*(U(1,J,K,1)-U(N1M,J,K,1))

	CALL CTDMA11(A11,B11,C11,D11,X11)

	DO I=2,N1M
      D11(I)=24.*DX1*(SGSVIS(I,J,K)*X11(I)-SGSVIS(I-1,J,K)*X11(I-1))
      ENDDO
	D11(1)=24.*DX1*(SGSVIS(1,J,K)*X11(1)-SGSVIS(N1M,J,K)*X11(N1M))

	CALL CTDMA11(A11,B11,C11,D11,X11)	

      DO I=1,N1M
      HH(I,J,K,1)=HH(I,J,K,1)-2.*X11(I)
      ENDDO
         	
  61  CONTINUE
      
C------------------D(SGS*(DU/DZ))/DZ------------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
      ENDDO
      

      DO 65 J=1,N2M
      
!$OMP PARALLEL DO PRIVATE (D1,X1)      
      DO 66 K=1,N3M
          
	DO  I=2,N1M
	D1(I)=4.*(SGSVIS(I,J,K)+SGSVIS(I-1,J,K))
	ENDDO
	D1(1)=4.*(SGSVIS(1,J,K)+SGSVIS(N1M,J,K))
      
	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M
	SGSVIS2(I,K)=X1(I)
	ENDDO

  66  CONTINUE
      
!$OMP PARALLEL DO PRIVATE (D3,X3,D33,X33)      
	DO 67 I=1,N1M
	
      DO K=2,N3M
	D3(K)=4.*(SGSVIS2(I,K)+SGSVIS2(I,K-1))
	ENDDO
	D3(1)=4.*(SGSVIS2(I,1)+SGSVIS2(I,N3M))
      
	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=2,N3M
	D33(K)=24.*DX3*(U(I,J,K,1)-U(I,J,K-1,1))
	ENDDO
	D33(1)=24.*DX3*(U(I,J,1,1)-U(I,J,N3M,1))
      
	CALL CTDMA33(A33,B33,C33,D33,X33)

C     D(SGS*(DU/DZ))/DZ
	DO K=1,N3M-1
	D33(K)=24.*DX3*(X33(K+1)*X3(K+1)-X33(K)*X3(K))
	ENDDO
	D33(N3M)=24.*DX3*(X33(1)*X3(1)-X33(N3M)*X3(N3M))
      
      CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
      HH(I,J,K,1)=HH(I,J,K,1)-X33(K)
      ENDDO
    
 67   CONTINUE
      
 65   CONTINUE
      
C---------------------D(SGS*(DU/DY))/DY-----------------------------
!$OMP PARALLEL DO
	DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
      C1(I)=1.0
      ENDDO
      

      DO 62 K=1,N3M
      
!$OMP  PARALLEL DO  PRIVATE (D1,X1)
	DO 63 J=1,N2M
      
	DO I=2,N1M
	D1(I)=4.*(SGSVIS(I,J,K)+SGSVIS(I-1,J,K))
	ENDDO
	D1(1)=4.*(SGSVIS(1,J,K)+SGSVIS(N1M,J,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M
	SGSVIS1(I,J)=X1(I)
      ENDDO
      
  63  CONTINUE

!$OMP  PARALLEL DO  PRIVATE(IP,IM,JP,JM,FJUM,FJUP,SGSJ2,SGSJ1
!$omp& ,UU22,UU11,SGSDUY)      
	DO 64 I=1,N1M
      IP=IPA(I)
	IM=IMA(I)
	
      DO J=1,N2M
	JP=J+1
	JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)

      SGSJ2=FJUP*1./H(JP)*(DY(J )/2.*SGSVIS1(I,JP)
     >  +DY(JP)/2.*SGSVIS1(I,J))
     >  +(1.-FJUP)*0.5*(SGSVIS(I,N2,K)+SGSVIS(IM,N2,K))
      SGSJ1=FJUM*1./H(J )*(DY(JM)/2.*SGSVIS1(I,J)
     >	+DY(J )/2.*SGSVIS1(I,JM))
     >  +(1.-FJUM)*0.5*(SGSVIS(I,0,K)+SGSVIS(IM,0,K))
      
c	UU22=FJUP*(U(I,JP,K,1)-U(I,J,K,1))/H(JP)
c     >       +(1.-FJUP)*(U(I,N2,K,1)-U(I,N2M,K,1))/(DY(N2M)/2.)
c	UU11=FJUM*(U(I,J,K,1)-U(I,JM,K,1))/H(J)
c     >       +(1.-FJUM)*(U(I,1,K,1)-U(I,0,K,1))/(DY(1)/2.)
      UU22=(U(I,JP,K,1)-U(I,J,K,1))/H(JP)
      UU11=(U(I,J,K,1)-U(I,JM,K,1))/H(J)
      
      SGSDUY=1./DY(J)*((1.-IFWALL*(1.-FJUP)*(ALY-1.))*SGSJ2*UU22
     >            -(1.-IFWALL*(1.-FJUM))*SGSJ1*UU11) 
      
      HH(I,J,K,1)=HH(I,J,K,1)-SGSDUY
      
      ENDDO    
      
  64  CONTINUE

  62  CONTINUE
      
C-------------------D(SGS*(DV/DX))/DY-------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
      ENDDO
      

      DO 68 K=1,N3M

!$OMP  PARALLEL DO PRIVATE(D11,X11,D1,X1)      
	DO 69 J=1,N2M
          
	DO I=2,N1M
      D11(I)=24.*DX1*(U(I,J,K,2)-U(I-1,J,K,2))
	ENDDO
	D11(1)=24.*DX1*(U(1,J,K,2)-U(N1M,J,K,2))

	CALL CTDMA11(A11,B11,C11,D11,X11)

	DO I=1,N1M
	UU12(I,J)=X11(I)
	ENDDO

	DO I=2,N1M
	D1(I)=4.*(SGSVIS(I,J,K)+SGSVIS(I-1,J,K))
	ENDDO
	D1(1)=4.*(SGSVIS(1,J,K)+SGSVIS(N1M,J,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M
	SGSVIS12(I,J)=X1(I)
      ENDDO
      
  69  CONTINUE

!$OMP  PARALLEL DO PRIVATE(IP,IM,JP,JM,FJUM,FJUP,SGSJ2,SGSJ1,SGSDUXY)            
	DO 70 I=1,N1M
	IP=IPA(I)
	IM=IMA(I)
      
	DO J=1,N2M
	JP=J+1
	JM=J-1
	FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      SGSJ2=FJUP*1.0/H(JP)*(DY(J )/2.*SGSVIS12(I,JP)+DY(JP)/2.
     >   *SGSVIS12(I,J))+(1.-FJUP)*0.5*(SGSVIS(I,N2,K)+SGSVIS(IM,N2,K))
      SGSJ1=FJUM*1.0/H(J )*(DY(JM)/2.*SGSVIS12(I,J)+DY(J )/2.	
     >  *SGSVIS12(I,JM))+(1.-FJUM)*0.5*(SGSVIS(I,0,K)+SGSVIS(IM,0,K))
      
c	SGSDUXY=FJUP*1.0/DY(J)*(SGSJ2*UU12(I,JP)-SGSJ1*UU12(I,J))
c     >     +(1.-FJUP)*1.0/DY(N2M)*(-SGSJ1*UU12(I,N2M))
c      SGSDUXY=1.0/DY(J)*(SGSJ2*UU12(I,JP)-SGSJ1*UU12(I,J))
      SGSDUXY=1.0/DY(J)*
     >            ((1.-IFWALL*(1.-FJUP)*(ALY-1.))*SGSJ2*UU12(I,JP)
     >                -(1.-IFWALL*(1.-FJUM))*SGSJ1*UU12(I,J))
      
      HH(I,J,K,1)=HH(I,J,K,1)-SGSDUXY
      
      ENDDO    
      
  70  CONTINUE

  68  CONTINUE
      
C-----------------D(SGS*(DW/DX))/DZ-------------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
      DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
      ENDDO      


C    SGSVIS
      DO 71 J=1,N2M
      FJUM=FJMU(J)
      FJUP=FJPA(J)

!$OMP  PARALLEL DO  PRIVATE (D1,X1,D11,X11)      
	DO 72 K=1,N3M
          
	DO I=2,N1M
      D1(I)=4.*(SGSVIS(I,J,K)+SGSVIS(I-1,J,K))
	ENDDO
	D1(1)=4.*(SGSVIS(1,J,K)+SGSVIS(N1M,J,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M
	SGSVIS13(I,K)=X1(I)
	ENDDO

	DO I=2,N1M
	D11(I)=24.*DX1*(U(I,J,K,3)-U(I-1,J,K,3))
	ENDDO
	D11(1)=24.*DX1*(U(1,J,K,3)-U(N1M,J,K,3))

	CALL CTDMA11(A11,B11,C11,D11,X11)

      DO I=1,N1M
	WW13(I,K)=X11(I)
      ENDDO
      
  72  CONTINUE

!$OMP PARALLEL DO PRIVATE (D3,X3,D33,X33)            
      DO 73 I=1,N1M
          
	DO K=2,N3M
	D3(K)=4.*(SGSVIS13(I,K)+SGSVIS13(I,K-1))
	ENDDO
	D3(1)=4.*(SGSVIS13(I,1)+SGSVIS13(I,N3M))

	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M-1
      D33(K)=24.*DX3*(WW13(I,K+1)*X3(K+1)-WW13(I,K)*X3(K))
	ENDDO
	D33(N3M)=24.*DX3*(WW13(I,1)*X3(1)-WW13(I,N3M)*X3(N3M))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
      HH(I,J,K,1)=HH(I,J,K,1)-X33(K)
	ENDDO

      DO K=1,N3M
      ! WALL MODEL
      BC_DOWN=-TAUW(I,K,1)/DY(1)*IFWALL
      BC_UP  =-TAUW(I,K,2)/DY(N2M)*IFWALL*(ALY-1.)      
      BC=(1.-FJUM)*BC_DOWN+(1.-FJUP)*BC_UP 
      
      HH(I,J,K,1)=HH(I,J,K,1)
      
      RUH1(I,J,K)=RUH1(I,J,K)-3./2.*HH(I,J,K,1)+1./2.*HH2(I,J,K,1)+BC   
c      write(*,*) RUH1(I,J,K)
      ENDDO

   73 CONTINUE
      
   71 CONTINUE
      
      RETURN
      END
c***************** RHS2 ***********************     
      SUBROUTINE RHS2(U,P,HH,HH2,SS,SS2,RUH2,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC(0:M1,0:M3,3,2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH4/DYM(M2),DYC(M2),DYP(M2)
      COMMON/TSTEP/NTST,DT
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES


      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL HH2(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL SS2(0:M1,0:M2,0:M3,3)
      REAL RUH2(0:M1,0:M2,0:M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL NUTI,NUTJ,NUTK
      
      REAL X1(M1-1),X11(M1-1),X3(M3-1),X33(M3-1)
      REAL VV1(M1-1,M2-1),UU2(M1-1,M2-1),WW2(M2-1,M3-1),VV3(M2-1,M3-1)
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1)
      REAL A11(M1-1),B11(M1-1),C11(M1-1),D11(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1)
      REAL A33(M3-1),B33(M3-1),C33(M3-1),D33(M3-1)
      REAL VV21(M1-1,M2-1),SGSVISV21(M1-1,M2-1),SGSVISV1(M1-1,M2-1)
      REAL SGSVIS3(M2-1,M3-1),SGSVISV23(M2-1,M3-1),VV32(M2-1,M3-1)
      REAL SGSVISV3(M2-1,M3-1),VV23(M2-1,M3-1),SGSVISV32(M2-1,M3-1)
      REAL VV12(M1-1,M2-1),SGSVISV12(M1-1,M2-1),SGSVISV22(M1-1,M2-1)
      REAL VV22,VV11,SGSDVY
      
!$omp  parallel do private(KP,KM,JP,JM,IP,IM,FJUM,FJUP,
!$omp& VISCOS,PRESSG2,V1,V2,SGSJ2,SGSJ1,VIS22)
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO 10 J=2,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMV(J)
      FJUP=FJPA(J)
      
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
     
C     VISCOS TERM
      VISCOS= 
     >       (1./RE)*(DYP(J)*U(I,JP,K,2)
     >                         -DYC(J)*U(I,J ,K,2)
     >                         +DYM(J)*U(I,JM,K,2))
      
      PRESSG2=(P(I,J,K)-P(I,JM,K))/H(J)      
      
      RUH2(I,J,K)=-PRESSG2+VISCOS

 
C     CONVECTION TERM (HH)  
      V2=0.5*(U(I,JP,K,2)+U(I,J,K,2))      
      V1=0.5*(U(I,J,K,2)+U(I,JM,K,2))      
      HH(I,J,K,2)=HH(I,J,K,2)+1./H(J)*(V2*V2-V1*V1)      
c      HH(I,J,K,2)=HH(I,J,K,2)+U(I,J,K,2)*1./H(J)*(V2-V1)/2.
      
C     D(SGS*DV/DY)/DY     
      SGSJ2=SGSVIS(I,J ,K)
      SGSJ1=SGSVIS(I,JM,K)
    
      V2=1./DY(J)*(U(I,JP,K,2)-U(I,J,K,2))
	V1=1./DY(JM)*(U(I,J,K,2)-U(I,JM,K,2))     
      VIS22=1.0/H(J)*(SGSJ2*V2-SGSJ1*V1)
      HH(I,J,K,2)=HH(I,J,K,2)-2.*VIS22

 10   CONTINUE 
      
CCCCC FOURTH COMPACT DIFFERENCE
      
CCC   VISCOS TERM: (1/RE)*D2V/DX2
      
!$OMP PARALLEL DO
      DO 15 I=1,N1M
      A1(I)=1.0
      B1(I)=10.0
      C1(I)=1.0
 15   CONTINUE      

!$OMP PARALLEL DO PRIVATE (D1,X1)
      DO 20 K=1,N3M
          
      DO 20 J=2,N2M     
          
      DO I=2,N1M-1
      D1(I)=12.*DX1Q*(U(I-1,J,K,2)-2.*U(I,J,K,2)+U(I+1,J,K,2))
      ENDDO
      D1(1)=12.*DX1Q*(U(N1M,J,K,2)-2.*U(1,J,K,2)+U(2,J,K,2))
      D1(N1M)=12.*DX1Q*(U(N1M-1,J,K,2)-2.*U(N1M,J,K,2)+U(1,J,K,2))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      RUH2(I,J,K)=RUH2(I,J,K)+(1./RE)*X1(I)
      ENDDO      
      
  20  CONTINUE
      
      
CCC   VISCOS TERM: (1/RE)*D2V/DZ2      
!$OMP PARALLEL DO
      DO 25 K=1,N3M
      A3(K)=1.0
      B3(K)=10.0
      C3(K)=1.0
  25  CONTINUE      

!$OMP PARALLEL DO PRIVATE (D3,X3) 
      DO 30 I=1,N1M
          
      DO 30 J=2,N2M     
      
      DO K=2,N3M-1
      D3(K)=12.*DX3Q*(U(I,J,K-1,2)-2.*U(I,J,K,2)+U(I,J,K+1,2))
      ENDDO
      D3(1)=12.*DX3Q*(U(I,J,N3M,2)-2.*U(I,J,1,2)+U(I,J,2,2))
      D3(N3M)=12.*DX3Q*(U(I,J,N3M-1,2)-2.*U(I,J,N3M,2)+U(I,J,1,2))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      RUH2(I,J,K)=RUH2(I,J,K)+(1./RE)*X3(K)
      ENDDO      
      
  30  CONTINUE     

      
CCC   CONVECTION TERM HH  
!$OMP PARALLEL DO
      DO I=1,N1M
      A1(I)=1.0
      B1(I)=6.0
      C1(I)=1.0
      ENDDO  
      
!$OMP PARALLEL DO
      DO I=1,N1M
      A11(I)=1.0
      B11(I)=22.0
      C11(I)=1.0
      ENDDO
      

      DO 50 K=1,N3M          
      
! INTERPOLATION OF V IN X DIRECTION
!$OMP PARALLEL DO PRIVATE (D1,X1)      
      DO 35 J=2,N2M
      
      DO I=2,N1M
      D1(I)=4.*(U(I-1,J,K,2)+U(I,J,K,2))
      ENDDO      
      D1(1)=4.*(U(N1M,J,K,2)+U(1,J,K,2))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      VV1(I,J)=X1(I)      
      ENDDO  
      
  35  CONTINUE
      
! INTERPOLATION OF U IN Y DIRECTION  
!$OMP PARALLEL DO     
      DO 40 I=1,N1M
      DO 40 J=2,N2M
          
      UU2(I,J)=1./H(J)*(DY(J)/2.*U(I,J-1,K,1)+DY(J-1)/2.*U(I,J,K,1))
          
  40  CONTINUE
      
C      DUV/DX
!$OMP PARALLEL DO PRIVATE (D11,X11)      
      DO 45 J=2,N2M
      
      DO I=1,N1M-1
      D11(I)=24.*DX1*(-UU2(I,J)*VV1(I,J)+UU2(I+1,J)*VV1(I+1,J))
      ENDDO      
      D11(N1M)=24.*DX1*(-UU2(N1M,J)*VV1(N1M,J)+UU2(1,J)*VV1(1,J))
      
      CALL CTDMA11(A11,B11,C11,D11,X11)
      
      DO I=1,N1M
      HH(I,J,K,2)=HH(I,J,K,2)+X11(I)     
      ENDDO 

C     (1/2) U DV/DX      
C      DO I=1,N1M-1
C      D1(I)  =4.*(UU2(I,J)+UU2(I+1,J))
C      ENDDO      
C      D1(N1M)=4.*(UU2(N1M,J)+UU2(1,J))
      
C      CALL CTDMA11(A1,B1,C1,D1,X1)
      
C      DO I=1,N1M-1
C      D11(I)  =24.*DX1*(-VV1(I,J)+VV1(I+1,J))
C      ENDDO      
C      D11(N1M)=24.*DX1*(-VV1(N1M,J)+VV1(1,J))
      
C      CALL CTDMA11(A11,B11,C11,D11,X11)
      
C      DO I=1,N1M
C      HH(I,J,K,2)=HH(I,J,K,2)+X1(I)*X11(I)/2.     
C      ENDDO 
      
  45  CONTINUE
      
  50  CONTINUE   
      
!$OMP PARALLEL DO
      DO K=1,N3M
      A3(K)=1.0
      B3(K)=6.0
      C3(K)=1.0
      ENDDO   
      
!$OMP PARALLEL DO
      DO K=1,N3M
      A33(K)=1.0
      B33(K)=22.0
      C33(K)=1.0
      ENDDO
      

      DO 60 I=1,N1M     
      
      ! INTERPOLATION OF V IN Z DIRECTION
!$OMP PARALLEL DO PRIVATE (D3,X3)       
      DO 55 J=2,N2M
      
      DO K=2,N3M
      D3(K)=4.*(U(I,J,K-1,2)+U(I,J,K,2))
      ENDDO      
      D3(1)=4.*(U(I,J,N3M,2)+U(I,J,1,2))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      VV3(J,K)=X3(K)
      ENDDO 
      
  55  CONTINUE
      
      ! INTERPOLATION OF W IN Y DIRECTION
!$OMP PARALLEL DO        
      DO 57 K=1,N3M
      DO 57 J=2,N2M          
      WW2(J,K)=1./H(J)*(DY(J)/2.*U(I,J-1,K,3)+DY(J-1)/2.*U(I,J,K,3))    
  57  CONTINUE
      
C      DVW/DZ
!$OMP PARALLEL DO PRIVATE (D33,X33)       
      DO 59 J=2,N2M
      
      DO K=1,N3M-1
      D33(K)=24.*DX3*(-VV3(J,K)*WW2(J,K)+VV3(J,K+1)*WW2(J,K+1))
      ENDDO      
      D33(N3M)=24.*DX3*(-VV3(J,N3M)*WW2(J,N3M)+VV3(J,1)*WW2(J,1))
      
      CALL CTDMA33(A33,B33,C33,D33,X33)      

      DO K=1,N3M
      HH(I,J,K,2)=HH(I,J,K,2)+X33(K)  
      ENDDO
      
C     (1/2) W DV/DZ
C      DO K=1,N3M-1
C      D3(K)  =4.*(WW2(J,K)+WW2(J,K+1))
C      ENDDO      
C      D3(N3M)=4.*(WW2(J,N3M)+WW2(J,1))
      
C      CALL CTDMA33(A3,B3,C3,D3,X3)
      
C      DO K=1,N3M-1
C      D33(K)  =24.*DX3*(-VV3(J,K)+VV3(J,K+1))
C      ENDDO      
C      D33(N3M)=24.*DX3*(-VV3(J,N3M)+VV3(J,1))
      
C      CALL CTDMA33(A33,B33,C33,D33,X33)  
      
C      DO K=1,N3M
C      HH(I,J,K,2)=HH(I,J,K,2)+X3(K)*X33(K)/2.  
C      ENDDO      
  
  59  CONTINUE     
  60  CONTINUE       


      
C------------------D(SGS*(DU/DY))/DX----------------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
      ENDDO


      DO K=1,N3M
      
!$OMP PARALLEL DO PRIVATE (D1,X1)      
	DO J=1,N2M
ccc
	DO I=2,N1M
      D1(I)=4.*(SGSVIS(I,J,K)+SGSVIS(I-1,J,K))
	ENDDO
      D1(1)=4.*(SGSVIS(1,J,K)+SGSVIS(N1M,J,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M
	SGSVISV1(I,J)=X1(I)
      ENDDO
      
      ENDDO
ccc
      
ccc
!$OMP PARALLEL DO PRIVATE (IM)      
      DO I=1,N1M
	IM=IMA(I)
      
	DO J=2,N2M          
	SGSVISV21(I,J)=1./H(J)
     >             *(DY(J-1)/2.*SGSVISV1(I,J)+DY(J)/2.*SGSVISV1(I,J-1))
	VV21(I,J)=1./H(J)*(U(I,J,K,1)-U(I,J-1,K,1))
      ENDDO
c      SGSVISV21(I,2)=0.5*(SGSVIS(I,1,K)+SGSVIS(IM,1,K))
c      VV21(I,2)=1./H(2)*(U(I,2,K,1)-U(I,1,K,1))
      ENDDO
ccc
      
ccc    
!$OMP PARALLEL DO PRIVATE (D11,X11)
      DO J=2,N2M
          
	DO I=1,N1M-1
      D11(I)=24.*DX1
     >  	   *(SGSVISV21(I+1,J)*VV21(I+1,J)-SGSVISV21(I,J)*VV21(I,J))
	ENDDO
	D11(N1M)=24.*DX1*(SGSVISV21(1,J)*VV21(1,J)
     >         -SGSVISV21(N1M,J)*VV21(N1M,J))

	CALL CTDMA11(A11,B11,C11,D11,X11)

	DO I=1,N1M
          HH(I,J,K,2)=HH(I,J,K,2)-X11(I)
      ENDDO
      
      ENDDO
ccc
      
	ENDDO

C-----------------------D(SGS*(DW/DY))/DZ---------------------------------

!$OMP PARALLEL DO
      DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
      ENDDO
      

      DO 75 I=1,N1M
	IM=IMA(I)
      
!$OMP PARALLEL DO PRIVATE (D3,X3)      
	DO 76 J=1,N2M
          
	DO K=2,N3M
      D3(K)=4.*(SGSVIS(I,J,K-1)+SGSVIS(I,J,K))
	ENDDO
	D3(1)=4.*(SGSVIS(I,J,N3M)+SGSVIS(I,J,1))

	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M
	SGSVIS3(J,K)=X3(K)
	ENDDO

  76  CONTINUE
      
!$OMP PARALLEL DO      
      DO 77 K=1,N3M
          
	DO J=2,N2M
      SGSVISV23(J,K)=1./H(J)
     >       *(DY(J)/2.*SGSVIS3(J-1,K)+DY(J-1)/2.*SGSVIS3(J,K))
	VV32(J,K)=1./H(J)*(U(I,J,K,3)-U(I,J-1,K,3))
	ENDDO
c	SGSVISV23(2,K)=0.5*(SGSVIS(I,1,K)+SGSVIS(IM,1,K))
c      VV32(2,K)=1./H(2)*(U(I,2,K,3)-U(I,1,K,3))
  77  CONTINUE
     
!$OMP PARALLEL DO PRIVATE (D33,X33)
      DO 78 J=2,N2M
          
	DO K=1,N3M-1
      D33(K)=24.*DX3*
     >       (SGSVISV23(J,K+1)*VV32(J,K+1)-SGSVISV23(J,K)*VV32(J,K))
      ENDDO
	D33(N3M)=24.*DX3*
     >       (SGSVISV23(J,1)*VV32(J,1)-SGSVISV23(J,N3M)*VV32(J,N3M))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
      HH(I,J,K,2)=HH(I,J,K,2)-X33(K)
	ENDDO

  78  CONTINUE
      
  75  CONTINUE
 

C-------------------D(SGS*(DV/DZ))/DZ--------------------------------
!$OMP PARALLEL DO
      DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO

!$OMP PARALLEL DO
	DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
	ENDDO


      DO 67 I=1,N1M
	
!$OMP PARALLEL DO PRIVATE(D3,X3,D33,X33)      
      DO 68 J=1,N2M
	
      DO K=2,N3M
      D3(K)=4.*(SGSVIS(I,J,K-1)+SGSVIS(I,J,K))
	ENDDO
	D3(1)=4.*(SGSVIS(I,J,N3M)+SGSVIS(I,J,1))

	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M
	SGSVISV3(J,K)=X3(K)
	ENDDO

	DO K=2,N3M
	D33(K)=24.*DX3*(U(I,J,K,2)-U(I,J,K-1,2))
	ENDDO
	D33(1)=24.*DX3*(U(I,J,1,2)-U(I,J,N3M,2))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
	VV23(J,K)=X33(K)
	ENDDO

  68  CONTINUE
      
!$OMP PARALLEL DO PRIVATE(KM,KP)      
      DO 69 K=1,N3M
	KM=KMA(K)
	KP=KPA(K)
      
	DO J=2,N2M
      SGSVISV32(J,K)=1./H(J)
     >      *(DY(J)/2.*SGSVISV3(J-1,K)+DY(J-1)/2.*SGSVISV3(J,K))
	ENDDO
c	SGSVISV32(2,K)=0.5*(SGSVIS(I,1,K)+SGSVIS(I,1,KM))
  69  CONTINUE
      
!$OMP PARALLEL DO PRIVATE(D33,X33)      
      DO 70 J=2,N2M
          
	DO K=1,N3M-1
      D33(K)=24.*DX3*
     >	   (SGSVISV32(J,K+1)*VV23(J,K+1)-SGSVISV32(J,K)*VV23(J,K))
	ENDDO
	D33(N3M)=24.*DX3*
     >	   (SGSVISV32(J,1)*VV23(J,1)-SGSVISV32(J,N3M)*VV23(J,N3M))

      CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
      HH(I,J,K,2)=HH(I,J,K,2)-X33(K)
	ENDDO

  70  CONTINUE
      
  67  CONTINUE

C----------------------D(SGS*(DV/DX))/DX----------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
      ENDDO

!$OMP PARALLEL DO
	DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
	ENDDO
      


      DO 61 K=1,N3M

!$OMP  PARALLEL DO PRIVATE (D11,X11)       
	DO 62 J=2,N2M 
          
      DO I=2,N1M
      D11(I)=24.*DX1*(U(I,J,K,2)-U(I-1,J,K,2))
	ENDDO
	D11(1)=24.*DX1*(U(1,J,K,2)-U(N1M,J,K,2))

	CALL CTDMA11(A11,B11,C11,D11,X11)

      DO I=1,N1M
	VV12(I,J)=X11(I)
      ENDDO
      
   62 CONTINUE

!$OMP  PARALLEL DO PRIVATE (D1,X1)       
      DO 63 J=1,N2M
	
      DO I=2,N1M
      D1(I)=4.*(SGSVIS(I,J,K)+SGSVIS(I-1,J,K))
	ENDDO
	D1(1)=4.*(SGSVIS(1,J,K)+SGSVIS(N1M,J,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M
	SGSVISV12(I,J)=X1(I)
      ENDDO
      
  63  CONTINUE

!$OMP  PARALLEL DO PRIVATE (IM,IP) 
      DO 64 I=1,N1M
	IM=IMA(I)
	IP=IPA(I)
      
	DO J=2,N2M
      SGSVISV22(I,J)=1.0/H(J)*(DY(J-1)/2.*SGSVISV12(I,J)
     >  +DY(J)/2.*SGSVISV12(I,J-1))
      ENDDO
c	SGSVISV22(I,2)=0.5*(SGSVIS(I,1,K)+SGSVIS(IM,1,K))
      
  64  CONTINUE

!$OMP  PARALLEL DO PRIVATE (D11,X11)       
      DO 65 J=2,N2M
          
	DO I=1,N1M-1
	D11(I)=24.*DX1*(SGSVISV22(I+1,J)*VV12(I+1,J)
     >      -SGSVISV22(I,J)*VV12(I,J))
	ENDDO
	D11(N1M)=24.*DX1*(SGSVISV22(1,J)*VV12(1,J)
     >      -SGSVISV22(N1M,J)*VV12(N1M,J))

	CALL CTDMA11(A11,B11,C11,D11,X11)

	DO I=1,N1M
	HH(I,J,K,2)=HH(I,J,K,2)-X11(I)
	ENDDO

      DO I=1,N1M
	RUH2(I,J,K)=RUH2(I,J,K)-3./2.*HH(I,J,K,2)+1./2.*HH2(I,J,K,2)
c      write(*,*) RUH2(I,J,K)
      ENDDO

  65  CONTINUE
      
  61  CONTINUE
      
      RETURN
      END

c***************** RHS3 ***********************     
      SUBROUTINE RHS3(U,P,HH,HH2,SS,SS2,RUH3,PRESG3,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/BCON/UBC(0:M1,0:M3,3,2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/WALLCO/HPW(M2),HMW(M2),HCW(M2)
      COMMON/FINDX3/FUP(M2),FDOWN(M2)


      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL HH(0:M1,0:M2,0:M3,3)
      REAL HH2(0:M1,0:M2,0:M3,3)
      REAL SS(0:M1,0:M2,0:M3,3)
      REAL SS2(0:M1,0:M2,0:M3,3)
      REAL RUH3(0:M1,0:M2,0:M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL NUTI,NUTJ,NUTK
      
      REAL X1(M1-1), X11(M1-1),X3(M3-1),X33(M3-1)
      REAL WW3(M3-1),WW1(M1-1,M3-1),UU3(M1-1,M3-1),VV3(M2-1,M3-1)
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1)
      REAL A11(M1-1),B11(M1-1),C11(M1-1),D11(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1)
      REAL A33(M3-1),B33(M3-1),C33(M3-1),D33(M3-1)
      REAL SGSVISW3(M1-1,M3-1)
      REAL SGSVISW2(M2-1,M3-1),WW22,WW11,SGSDWY
	REAL SGSVIS31(M1-1,M3-1),UU31(M1-1,M3-1)
	REAL VV32(M2-1,M3-1),SGSVIS32(M2-1,M3-1),SGSDVY

!$omp  parallel do private(KP,KM,JP,JM,IP,IM,FJUM,FJUP,
!$omp& VISCOS,PRESSG3,BC)
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
    
C     VISCOS TERM
      VISCOS=
     >                (1./RE)*(HP(J)*U(I,JP,K,3)
     >                        -HC(J)*U(I,J ,K,3)
     >                        +HM(J)*U(I,JM,K,3)) 

    
      RUH3(I,J,K)=VISCOS

      
     
C      RUH3(I,J,K)=RUH3(I,J,K)+BC

 10   CONTINUE 
      
      
CCCCC FOURTH COMPACT DIFFERENCE
      
CCC   VISCOS TERM: (1/RE)*D2W/DX2
!$OMP PARALLEL DO
      DO 15 I=1,N1M
      A1(I)=1.0
      B1(I)=10.0
      C1(I)=1.0
 15   CONTINUE      

!$OMP PARALLEL DO PRIVATE (D1,X1)
      DO 20 K=1,N3M
      DO 20 J=1,N2M     
     
      DO I=2,N1M-1
      D1(I)=12.*DX1Q*(U(I-1,J,K,3)-2.*U(I,J,K,3)+U(I+1,J,K,3))
      ENDDO
      D1(1)=12.*DX1Q*(U(N1M,J,K,3)-2.*U(1,J,K,3)+U(2,J,K,3))
      D1(N1M)=12.*DX1Q*(U(N1M-1,J,K,3)-2.*U(N1M,J,K,3)+U(1,J,K,3))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      RUH3(I,J,K)=RUH3(I,J,K)+(1./RE)*X1(I)
      ENDDO      
      
  20  CONTINUE
      
      
CCC   VISCOS TERM: (1/RE)*D2W/DZ2      
!$OMP PARALLEL DO
      DO 25 K=1,N3M
      A3(K)=1.0
      B3(K)=10.0
      C3(K)=1.0
 25   CONTINUE      

!$OMP PARALLEL DO PRIVATE (D3,X3) 
      DO 30 I=1,N1M
      DO 30 J=1,N2M     
   
      DO K=2,N3M-1
      D3(K)=12.*DX3Q*(U(I,J,K-1,3)-2.*U(I,J,K,3)+U(I,J,K+1,3))
      ENDDO
      D3(1)=12.*DX3Q*(U(I,J,N3M,3)-2.*U(I,J,1,3)+U(I,J,2,3))
      D3(N3M)=12.*DX3Q*(U(I,J,N3M-1,3)-2.*U(I,J,N3M,3)+U(I,J,1,3))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      RUH3(I,J,K)=RUH3(I,J,K)+(1./RE)*X3(K)
      ENDDO      
      
  30  CONTINUE
      
CCC   PRESSURE GRADIENT TERM
!$OMP PARALLEL DO
      DO 35 K=1,N3M
      A3(K)=1.0
      B3(K)=22.0
      C3(K)=1.0
 35   CONTINUE
      
!$OMP PARALLEL DO PRIVATE (D3,X3)
      DO 40 I=1,N1M
      DO 40 J=1,N2M   
      
      DO K=2,N3M
      D3(K)=24.*DX3*(-P(I,J,K-1)+P(I,J,K))
      ENDDO      
      D3(1)=24.*DX3*(-P(I,J,N3M)+P(I,J,1))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      RUH3(I,J,K)=RUH3(I,J,K)-X3(K)-PRESG3
      ENDDO      
      
  40  CONTINUE   
      
CCC   CONVECTION TERM HH
      
C     D(WV)/DY
!$OMP PARALLEL DO
      DO K=1,N3M
      A3(K)=1.0
      B3(K)=6.0
      C3(K)=1.0
      ENDDO 


      DO I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      
!$OMP  PARALLEL DO PRIVATE(D3,X3)      
      DO J=1,N2M
      
      DO K=2,N3M
      D3(K) = 4.0*(U(I,J,K-1,2)+U(I,J,K,2))
      ENDDO        
      D3(1) = 4.0*(U(I,J,N3M,2)+U(I,J,1,2))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)
      
      DO K=1,N3M
      VV3(J,K)=X3(K)
      ENDDO       
      
      ENDDO

!$OMP  PARALLEL DO PRIVATE(KP,KM,JP,JM,FJUM,FJUP,W2,W1,V2,V1)
      DO K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      DO J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)

C      D(WV)/DY
      W2=FJUP*1.0/H(JP)*(DY(JP)/2.*U(I,J,K,3)+DY(J)/2.*U(I,JP,K,3))
     >  +(1.-FJUP)*U(I,N2,K,3)
      W1=FJUM*1.0/H(J)*(DY(J)/2.*U(I,JM,K,3)+DY(JM)/2.*U(I,J,K,3))
     >  +(1.-FJUM)*U(I,0,K,3)
      
      V2 = VV3(JP,K)
      V1 = VV3(J,K)
      IF (J.EQ.N2M) V2=0.0
      IF (J.EQ.1)   V1=0.0

      HH(I,J,K,3)=HH(I,J,K,3)+1./DY(J)*(V2*W2-V1*W1)
      
C     (1/2) V DW/DY
c      HH(I,J,K,3)=HH(I,J,K,3)+(V2+V1)/2.*1./DY(J)*(W2-W1)/2.
      
      ENDDO
      ENDDO
      
      ENDDO
  
      
C     D(WW)/DZ
!$OMP PARALLEL DO
      DO K=1,N3M
      A3(K)=1.0
c      B3(K)=4.0
      B3(K)=6.0
      C3(K)=1.0
      ENDDO 
      
!$OMP PARALLEL DO
      DO K=1,N3M
      A33(K)=1.0
      B33(K)=22.0
      C33(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO PRIVATE (D3,X3,WW3,D33,X33) 
      DO 50 I=1,N1M
      DO 50 J=1,N2M     
      
C      ! INTERPOLATION OF W IN Z DIRECTION
      DO K=1,N3M-1
      D3(K)=4.*(U(I,J,K,3)+U(I,J,K+1,3))
      ENDDO      
      D3(N3M)=4.*(U(I,J,N3M,3)+U(I,J,1,3))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      WW3(K)=X3(K)
      ENDDO   
      
C      DWW/DZ
      DO K=2,N3M
      D33(K)=24.*DX3*(WW3(K)**2-WW3(K-1)**2)
      ENDDO      
      D33(1)=24.*DX3*(WW3(1)**2-WW3(N3M)**2)      
      
      CALL CTDMA33(A33,B33,C33,D33,X33)
      
      DO K=1,N3M
      HH(I,J,K,3)=HH(I,J,K,3)+X33(K)     
      ENDDO     
      
C     (1/2) W DW/DZ
C      DO K=2,N3M
C      D33(K)=24.*DX3*(WW3(K)-WW3(K-1))
C      ENDDO      
C      D33(1)=24.*DX3*(WW3(1)-WW3(N3M))      
      
C      CALL CTDMA33(A33,B33,C33,D33,X33)
      
C      DO K=1,N3M
C      HH(I,J,K,3)=HH(I,J,K,3)+U(I,J,K,3)*X33(K)/2.      
C      ENDDO       
      
  50  CONTINUE
      
    
C     D(UW)/DX
!$OMP PARALLEL DO
      DO K=1,N3M
      A3(K)=1.0
      B3(K)=6.0
      C3(K)=1.0 
      ENDDO 
      
!$OMP PARALLEL DO
      DO I=1,N1M
      A1(I)=1.0
      B1(I)=6.0
      C1(I)=1.0
      ENDDO   
      
!$OMP PARALLEL DO
      DO I=1,N1M
      A11(I)=1.0
      B11(I)=22.0
      C11(I)=1.0   
      ENDDO
      

      DO 60 J=1,N2M     
      
      ! INTERPOLATION OF U IN Z DIRECTION
!$OMP PARALLEL DO PRIVATE (D3,X3)      
      DO 55 I=1,N1M
      
      DO K=2,N3M
      D3(K)=4.*(U(I,J,K-1,1)+U(I,J,K,1))
      ENDDO      
      D3(1)=4.*(U(I,J,N3M,1)+U(I,J,1,1))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      UU3(I,K)=X3(K)
      ENDDO 
      
  55  CONTINUE
      
      ! INTERPOLATION OF W IN X DIRECTION
!$OMP PARALLEL DO PRIVATE (D1,X1)      
      DO 57 K=1,N3M
      
      DO I=2,N1M
      D1(I)=4.*(U(I-1,J,K,3)+U(I,J,K,3))
      ENDDO      
      D1(1)=4.*(U(N1M,J,K,3)+U(1,J,K,3))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      WW1(I,K)=X1(I)
      ENDDO 
      
  57  CONTINUE
      
!$OMP PARALLEL DO PRIVATE (D11,X11)
      DO 59 K=1,N3M
      
C      D(UW)/DX
      DO I=1,N1M-1
      D11(I)=24.*DX1*(-UU3(I,K)*WW1(I,K)+UU3(I+1,K)*WW1(I+1,K))
      ENDDO      
      D11(N1M)=24.*DX1*(-UU3(N1M,K)*WW1(N1M,K)+UU3(1,K)*WW1(1,K))
      
      CALL CTDMA11(A11,B11,C11,D11,X11)      

      DO I=1,N1M
      HH(I,J,K,3)=HH(I,J,K,3)+X11(I) 
      ENDDO

C     (1/2) U D(W)/DX
C      DO I=1,N1M-1
C      D1(I)  =4.*(UU3(I,K)+UU3(I+1,K))
C      ENDDO      
C      D1(N1M)=4.*(UU3(N1M,K)+UU3(1,K))
      
C      CALL CTDMA11(A1,B1,C1,D1,X1)
      
C      DO I=1,N1M-1
C      D11(I)  =24.*DX1*(-WW1(I,K)+WW1(I+1,K))
C      ENDDO      
C      D11(N1M)=24.*DX1*(-WW1(N1M,K)+WW1(1,K))
      
C      CALL CTDMA11(A11,B11,C11,D11,X11) 
      
C      DO I=1,N1M
C      HH(I,J,K,3)=HH(I,J,K,3)+X1(I)*X11(I)/2.  
C      ENDDO
    
      
  59  CONTINUE
      
  60  CONTINUE
      
      
C-----------------D(SGS*(DW/DX))/DX-------------------------------------
!$OMP PARALLEL DO
      DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
      ENDDO
      

      DO 61 J=1,N2M
	
!$OMP PARALLEL DO PRIVATE (D3,X3)      
      DO 62 I=1,N1M
	
      DO K=2,N3M
	D3(K)=4.*(SGSVIS(I,J,K-1)+SGSVIS(I,J,K))
	ENDDO
	D3(1)=4.*(SGSVIS(I,J,N3M)+SGSVIS(I,J,1))  
      
	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M
	SGSVISW3(I,K)=X3(K)
      ENDDO
      
  62  CONTINUE
      
!$OMP PARALLEL DO PRIVATE (D1,X1,D11,X11)      
      DO 63 K=1,N3M
	
      DO I=2,N1M
	D1(I)=4.*(SGSVISW3(I-1,K)+SGSVISW3(I,K))
	ENDDO
	D1(1)=4.*(SGSVISW3(N1M,K)+SGSVISW3(1,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=2,N1M
	D11(I)=24.*DX1*(U(I,J,K,3)-U(I-1,J,K,3))
	ENDDO
	D11(1)=24.*DX1*(U(1,J,K,3)-U(N1M,J,K,3))

	CALL CTDMA11(A11,B11,C11,D11,X11)
	
	DO I=1,N1M-1
	D11(I)=24.*DX1*
     >	    (X1(I+1)*X11(I+1)-X1(I)*X11(I))
	ENDDO
      D11(N1M)=24.*DX1*
     >        (X1(1)*X11(1)-X1(N1M)*X11(N1M))

	CALL CTDMA11(A11,B11,C11,D11,X11)

	DO I=1,N1M
	HH(I,J,K,3)=HH(I,J,K,3)-X11(I)
      ENDDO  
      
  63  CONTINUE
      
  61  CONTINUE
      
C----------------D(SGS*(DW/DY))/DY-------------------------------------
!$OMP PARALLEL DO
      DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      

	DO 64 I=1,N1M
	IM=IMA(I)
	IP=IPA(I)
	
!$OMP  PARALLEL DO PRIVATE (D3,X3)      
      DO 65 J=1,N2M
	
      DO K=2,N3M
      D3(K)=4.*(SGSVIS(I,J,K)+SGSVIS(I,J,K-1))
	ENDDO
	D3(1)=4.*(SGSVIS(I,J,1)+SGSVIS(I,J,N3M))

	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M
	SGSVISW2(J,K)=X3(K)
      ENDDO
      
  65  CONTINUE

!$OMP  PARALLEL DO PRIVATE (JP,JM,
!$OMP& KM,KP,FJUM,FJUP,WW22,WW11,SGSJ1,SGSJ2,SGSDWY)      
      DO 66 K=1,N3M
	KM=KMA(K)
	KP=KPA(K)
      
	DO J=1,N2M
      JP=J+1
	JM=J-1
	FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      WW22=1./H(JP)*(U(I,JP,K,3)-U(I,J,K,3))
	WW11=1./H(J)*(U(I,J,K,3)-U(I,JM,K,3))
      
	SGSJ1=FJUM*1./H(J)
     >	  *(DY(JM)/2.*SGSVISW2(J,K)+DY(J)/2.*SGSVISW2(JM,K))
     >	  +(1.-FJUM)*0.5*(SGSVIS(I,0,K)+SGSVIS(I,0,KM))
      SGSJ2=FJUP*1.0/H(JP )*(DY(J)/2.*SGSVISW2(JP,K)+DY(JP )/2.	
     >   *SGSVISW2(J,K))+(1.-FJUP)*0.5*(SGSVIS(I,N2,K)+SGSVIS(I,N2,KM))
      
	SGSDWY=1./DY(J)*((1.-IFWALL*(1.-FJUP)*(ALY-1.))*SGSJ2*WW22
     >                    -(1.-IFWALL*(1.-FJUM))*SGSJ1*WW11)
      
      HH(I,J,K,3)=HH(I,J,K,3)-SGSDWY
      
      ENDDO
      
  66  CONTINUE
      
  64  CONTINUE
      
C-----------------------D(SGS*(DW/DZ))/DZ-----------------------------
!$OMP PARALLEL DO
	DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
	ENDDO
      
!$OMP PARALLEL DO PRIVATE (D33,X33)
	DO 67 I=1,N1M
	DO 67 J=1,N2M
	
      DO K=1,N3M-1
      D33(K)=24.*DX3*(U(I,J,K+1,3)-U(I,J,K,3))
	ENDDO
	D33(N3M)=24.*DX3*(U(I,J,1,3)-U(I,J,N3M,3))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=2,N3M
	D33(K)=24.*DX3*(X33(K)*SGSVIS(I,J,K)-X33(K-1)*SGSVIS(I,J,K-1))
	ENDDO
	D33(1)=24.*DX3*(X33(1)*SGSVIS(I,J,1)-X33(N3M)*SGSVIS(I,J,N3M))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
      HH(I,J,K,3)=HH(I,J,K,3)-2.*X33(K)
      ENDDO
      
  67  CONTINUE
      
C-------------------D(SGS*(DU/DZ))/DX--------------------------------------
!$OMP PARALLEL DO
      DO I=1,N1M
	A1(I)=1.0
	B1(I)=6.0
	C1(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO I=1,N1M
	A11(I)=1.0
	B11(I)=22.0
	C11(I)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
      ENDDO
      

	DO 68 J=1,N2M
	
!$OMP PARALLEL DO PRIVATE (D33,X33,D3,X3)      
      DO 69 I=1,N1M
          
	DO K=2,N3M
	D33(K)=24.*DX3*(U(I,J,K,1)-U(I,J,K-1,1))
	ENDDO
	D33(1)=24.*DX3*(U(I,J,1,1)-U(I,J,N3M,1))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
	UU31(I,K)=X33(K)
      ENDDO
      
	DO K=2,N3M
	D3(K)=4.*(SGSVIS(I,J,K)+SGSVIS(I,J,K-1))
	ENDDO
	D3(1)=4.*(SGSVIS(I,J,1)+SGSVIS(I,J,N3M))

	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M
	SGSVIS31(I,K)=X3(K)
      ENDDO
      
  69  CONTINUE

!$OMP PARALLEL DO PRIVATE (D1,X1,D11,X11)      
      DO 70 K=1,N3M
	
      DO I=2,N1M
	D1(I)=4.*(SGSVIS31(I,K)+SGSVIS31(I-1,K))
	ENDDO
	D1(1)=4.*(SGSVIS31(1,K)+SGSVIS31(N1M,K))

	CALL CTDMA11(A1,B1,C1,D1,X1)

	DO I=1,N1M-1
	D11(I)=24.*DX1*(X1(I+1)*UU31(I+1,K)-X1(I)*UU31(I,K))
	ENDDO
	D11(N1M)=24.*DX1*(X1(1)*UU31(1,K)-X1(N1M)*UU31(N1M,K))

	CALL CTDMA11(A11,B11,C11,D11,X11)

	DO I=1,N1M
      HH(I,J,K,3)=HH(I,J,K,3)-X11(I)
      ENDDO
      
  70  CONTINUE
      
  68  CONTINUE
      
C----------------D(SGS*(DV/DZ))/DY-------------------------------
!$OMP PARALLEL DO
      DO K=1,N3M
	A33(K)=1.0
	B33(K)=22.0
	C33(K)=1.0
      ENDDO
      
!$OMP PARALLEL DO
	DO K=1,N3M
	A3(K)=1.0
	B3(K)=6.0
	C3(K)=1.0
      ENDDO
      

      DO 71 I=1,N1M
	IM=IMA(I)
	IP=IPA(I)
      
!$OMP  PARALLEL DO  PRIVATE (D33,X33,D3,X3)      
	DO 72 J=1,N2M
          
	DO K=2,N3M
      D33(K)=24.*DX3*(U(I,J,K,2)-U(I,J,K-1,2))
	ENDDO
	D33(1)=24.*DX3*(U(I,J,1,2)-U(I,J,N3M,2))

	CALL CTDMA33(A33,B33,C33,D33,X33)

	DO K=1,N3M
	VV32(J,K)=X33(K)
      ENDDO
      
	DO K=2,N3M
	D3(K)=4.*(SGSVIS(I,J,K)+SGSVIS(I,J,K-1))
	ENDDO
	D3(1)=4.*(SGSVIS(I,J,1)+SGSVIS(I,J,N3M))

	CALL CTDMA33(A3,B3,C3,D3,X3)

	DO K=1,N3M
	SGSVIS32(J,K)=X3(K)
      ENDDO
      
  72  CONTINUE
      
!$OMP  PARALLEL DO  PRIVATE (KM,KP,JM,JP,FJUM,FJUP,SGSJ2,SGSJ1,SGSDVY)      
	DO 73 K=1,N3M
	KM=KMA(K)
	KP=KPA(K)
      
	DO J=1,N2M
	JP=J+1
	JM=J-1
	FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      SGSJ2=FJUP*1.0/H(JP)*(DY(J )/2.*SGSVIS32(JP,K)+DY(JP)/2.
     >  *SGSVIS32(J,K))+(1.-FJUP)*0.5*(SGSVIS(I,N2,K)+SGSVIS(I,N2,KM))
      SGSJ1=FJUM*1.0/H(J )*(DY(JM)/2.*SGSVIS32(J,K)+DY(J )/2.	
     >  *SGSVIS32(JM,K))+(1.-FJUM)*0.5*(SGSVIS(I,0,K)+SGSVIS(I,0,KM))
      
c	SGSDVY=FJUP*1.0/DY(J)*(SGSJ2*VV32(JP,K)-SGSJ1*VV32(J,K))
c     >     +(1.-FJUP)*1.0/DY(N2M)*(-SGSJ1*VV32(N2M,K))
      SGSDVY=1.0/DY(J)*((1.-IFWALL*(1.-FJUP)*(ALY-1.))*SGSJ2*VV32(JP,K)
     >                    -(1.-IFWALL*(1.-FJUM))*SGSJ1*VV32(J,K))
      
      HH(I,J,K,3)=HH(I,J,K,3)-SGSDVY
      
      ! WALL MODEL
      BC_DOWN=-TAUB(I,K,1)/DY(1)*IFWALL
      BC_UP  =-TAUB(I,K,2)/DY(N2M)*IFWALL*(ALY-1.)
      BC=(1.-FJUM)*BC_DOWN+(1.-FJUP)*BC_UP
      
      HH(I,J,K,3)=HH(I,J,K,3)
      
      RUH3(I,J,K)=RUH3(I,J,K)-3./2.*HH(I,J,K,3)+1./2.*HH2(I,J,K,3)+BC
      
c      write(*,*)  RUH3(I,J,K)
      
      ENDDO
      
  73  CONTINUE
      
  71  CONTINUE
        
      RETURN
      END

c***************** RHSDP ***********************     
      SUBROUTINE RHSDP(RDP,UH,U,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/BCON/UBC(0:M1,0:M3,3,2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)

      REAL U(0:M1,0:M2,0:M3,3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL RDP(M1,M2,M3)
      
      REAL X1(M1-1),X3(M3-1)
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1)
      
C     DIV(V)
!$omp parallel do private(KP,KM,JP,JM,IP,IM,FJUM,FJUP,DIVUH,CBC)
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      
      DIVUH=
C     > (UH(IP,J,K,1)-UH(I,J,K,1))*DX1+
     >   (FJUP*UH(I,JP,K,2)-FJUM*UH(I,J,K,2))/DY(J)
C     >   +(UH(I,J,KP,3)-UH(I,J,K,3))*DX3

      CBC=(1.-FJUM)*UBC(I,K,2,1)/DY(J)
     >       -(1.-FJUP)*UBC(I,K,2,2)/DY(J)

      RDP(I,J,K)=(DIVUH-CBC)/DT
      
  10  CONTINUE
      
!$OMP PARALLEL DO
      DO 15 I=1,N1M
      A1(I)=1.0
      B1(I)=22.0
      C1(I)=1.0
 15   CONTINUE      

C     DIV(U)
!$OMP PARALLEL DO PRIVATE (D1,X1)
      DO 20 K=1,N3M
      DO 20 J=1,N2M     

      DO I=1,N1M-1
      D1(I)=24.*DX1*(UH(I+1,J,K,1)-UH(I,J,K,1))
      ENDDO
      D1(N1M)=24.*DX1*(UH(1,J,K,1)-UH(N1M,J,K,1))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      RDP(I,J,K)=RDP(I,J,K)+X1(I)/DT      
      ENDDO      
      
  20  CONTINUE

C     DIV(W)
!$OMP PARALLEL DO
      DO 25 K=1,N3M
      A3(K)=1.0
      B3(K)=22.0
      C3(K)=1.0
  25  CONTINUE 
      
!$OMP PARALLEL DO PRIVATE (D3,X3)
      DO 30 I=1,N1M
      DO 30 J=1,N2M     
  
      DO K=1,N3M-1
      D3(K)=24.*DX3*(UH(I,J,K+1,3)-UH(I,J,K,3))
      ENDDO
      D3(N3M)=24.*DX3*(UH(I,J,1,3)-UH(I,J,N3M,3))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      RDP(I,J,K)=RDP(I,J,K)+X3(K)/DT
      ENDDO      
      
  30  CONTINUE
      
      RETURN 
      END


c***************** GETUH1 ***********************     
      Subroutine GETUH1(U,UH,RUH1,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/WALLCO/HPW(M2),HMW(M2),HCW(M2)
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL RUH1(0:M1,0:M2,0:M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL NUTI,NUTJ,NUTK
      REAL API(M1,M2,M3),ACI(M1,M2,M3),AMI(M1,M2,M3)
      REAL APJ(M1,M2,M3),ACJ(M1,M2,M3),AMJ(M1,M2,M3)
      REAL APK(M1,M3,M2),ACK(M1,M3,M2),AMK(M1,M3,M2)
      REAL R1(M1,M2,M3),R2(M1,M2,M3),R3(M1,M3,M2)
C      EQUIVALENCE(API,APJ)
C      EQUIVALENCE(ACI,ACJ)
C      EQUIVALENCE(AMI,AMJ)
C      EQUIVALENCE(R1,R2)


!$omp parallel do private(JP,JM,IP,IM,FJUM,FJUP)    
      DO 2 K=1,N3M    
  
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)      
      DO 20 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

      APJ(I,J,K)=FJUP*(
     >      -0.5*HP(J)*(1./RE))*DT    
      ACJ(I,J,K)=1.0+(
     >       0.5*HC(J)*(1./RE))*DT
       
      AMJ(I,J,K)=FJUM*(
     >      -0.5*HM(J)*(1./RE))*DT     

      R2(I,J,K)=RUH1(I,J,K)*DT
      
  20  CONTINUE
      
      CALL TDMAI(AMJ(1,1,K),ACJ(1,1,K),APJ(1,1,K),R2(1,1,K),
     >        R2(1,1,K),1,N2M,1,N1M)
      
      DO 21 J=1,N2M
      DO 21 I=1,N1M
      UH(I,J,K,1)=R2(I,J,K)     
  21  CONTINUE     
      
  
  2   CONTINUE


!$omp parallel do private(IP,IM)      
      DO 1 K=1,N3M      

      DO 10 J=1,N2M
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
    
   
      API(I,J,K)=1.0  -DT/2.0*(1./RE)*12.0*DX1Q*( 1.0)
      ACI(I,J,K)=10.0 -DT/2.0*(1./RE)*12.0*DX1Q*(-2.0)
      AMI(I,J,K)=1.0  -DT/2.0*(1./RE)*12.0*DX1Q*( 1.0)

      R1(I,J,K)=UH(I-1,J,K,1)+10.0*UH(I,J,K,1)+UH(I+1,J,K,1)      
      IF (I.EQ.1) R1(I,J,K)=UH(N1M,J,K,1)+10.0*UH(I,J,K,1)+UH(I+1,J,K,1)
      IF (I.EQ.N1M) R1(I,J,K)=UH(I-1,J,K,1)+10.0*UH(I,J,K,1)+UH(1,J,K,1)
      
  10  CONTINUE

      CALL CTDMA1J(AMI(1,1,K),ACI(1,1,K),API(1,1,K),R1(1,1,K),
     >          R1(1,1,K),1,N1M,1,N2M)

      DO 11 J=1,N2M 
      DO 11 I=1,N1M 
      UH(I,J,K,1)=R1(I,J,K)
  11  CONTINUE
  1   CONTINUE


!$omp parallel do private(IP,IM,KP,KM)      
      DO 3 J=1,N2M

      DO 30 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)      

      APK(I,K,J)=1.0  -DT/2.0*(1./RE)*12.0*DX3Q*( 1.0)
      ACK(I,K,J)=10.0 -DT/2.0*(1./RE)*12.0*DX3Q*(-2.0)
      AMK(I,K,J)=1.0  -DT/2.0*(1./RE)*12.0*DX3Q*( 1.0)

      R3(I,K,J)=UH(I,J,K-1,1)+10.0*UH(I,J,K,1)+UH(I,J,K+1,1)      
      IF (K.EQ.1) R3(I,K,J)=UH(I,J,N3M,1)+10.0*UH(I,J,K,1)+UH(I,J,K+1,1)
      IF (K.EQ.N3M) R3(I,K,J)=UH(I,J,K-1,1)+10.0*UH(I,J,K,1)+UH(I,J,1,1)

  30  CONTINUE
      
      CALL CTDMA3I(AMK(1,1,J),ACK(1,1,J),APK(1,1,J),
     >         R3(1,1,J),UH(0,0,0,1),J,N3M,1,N1M)
      
      
  3   CONTINUE

      RETURN
      END

c***************** GETUH2 ***********************     
      SUBROUTINE GETUH2(U,UH,RUH2,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH4/DYM(M2),DYC(M2),DYP(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES

      REAL U(0:M1,0:M2,0:M3,3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL RUH2(0:M1,0:M2,0:M3)
      REAL NUTI,NUTJ,NUTK
      REAL API(M1,M2,M3),ACI(M1,M2,M3),AMI(M1,M2,M3)
      REAL APJ(M1,M2,M3),ACJ(M1,M2,M3),AMJ(M1,M2,M3)
      REAL APK(M1,M3,M2),ACK(M1,M3,M2),AMK(M1,M3,M2)
      REAL R1(M1,M2,M3),R2(M1,M2,M3),R3(M1,M3,M2)
C      EQUIVALENCE(API,APJ)
C      EQUIVALENCE(ACI,ACJ)
C      EQUIVALENCE(AMI,AMJ)
C      EQUIVALENCE(R1,R2)


!$omp parallel do private(JP,JM,IP,IM,KP,KM,FJUM,FJUP)     
      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
 
      DO 20 J=2,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMV(J)
      FJUP=FJPA(J)
      DO 20 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

     
      APJ(I,J,K)=FJUP*(
     >      -0.5*DYP(J)*(1./RE)
     >      )*DT
      ACJ(I,J,K)=1.0+(
     >      +0.5*DYC(J)*(1./RE)
     >      )*DT
      AMJ(I,J,K)=FJUM*(
     >      -0.5*DYM(J)*(1./RE)
     >      )*DT


      R2(I,J,K)=DT*RUH2(I,J,K)

  20  CONTINUE

      CALL TDMAI(AMJ(1,1,K),ACJ(1,1,K),APJ(1,1,K),R2(1,1,K),
     >           R2(1,1,K),2,N2M,1,N1M)

      DO 21 J=2,N2M  
      DO 21 I=1,N1M
      UH(I,J,K,2)=R2(I,J,K)  
  21  CONTINUE  
  2   CONTINUE


!$omp parallel do private(JP,JM,IP,IM,KP,KM)   
      DO 1 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)      
   
      DO 10 J=2,N2M
      JP=J+1
      JM=J-1
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

      API(I,J,K)=1.0  -DT/2.0*(1./RE)*12.0*DX1Q*( 1.0)
      ACI(I,J,K)=10.0 -DT/2.0*(1./RE)*12.0*DX1Q*(-2.0)
      AMI(I,J,K)=1.0  -DT/2.0*(1./RE)*12.0*DX1Q*( 1.0)

      R1(I,J,K)=UH(I-1,J,K,2)+10.0*UH(I,J,K,2)+UH(I+1,J,K,2)      
      IF (I.EQ.1) R1(I,J,K)=UH(N1M,J,K,2)+10.0*UH(I,J,K,2)+UH(I+1,J,K,2)
      IF (I.EQ.N1M) R1(I,J,K)=UH(I-1,J,K,2)+10.0*UH(I,J,K,2)+UH(1,J,K,2)

  10  CONTINUE

      CALL CTDMA1J(AMI(1,1,K),ACI(1,1,K),API(1,1,K),
     >          R1(1,1,K),R1(1,1,K),1,N1M,2,N2M)

      DO 11 J=2,N2M 
      DO 11 I=1,N1M 
      UH(I,J,K,2)=R1(I,J,K)
  11  CONTINUE
  1   CONTINUE
  

!$omp parallel do private(KP,KM,IP,IM,JP,JM)    
      DO 3 J=2,N2M
      JP=J+1
      JM=J-1      
  
      DO 30 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      
      APK(I,K,J)=1.0  -DT/2.0*(1./RE)*12.0*DX3Q*( 1.0)
      ACK(I,K,J)=10.0 -DT/2.0*(1./RE)*12.0*DX3Q*(-2.0)
      AMK(I,K,J)=1.0  -DT/2.0*(1./RE)*12.0*DX3Q*( 1.0)

      R3(I,K,J)=UH(I,J,K-1,2)+10.0*UH(I,J,K,2)+UH(I,J,K+1,2)      
      IF (K.EQ.1) R3(I,K,J)=UH(I,J,N3M,2)+10.0*UH(I,J,K,2)+UH(I,J,K+1,2)
      IF (K.EQ.N3M) R3(I,K,J)=UH(I,J,K-1,2)+10.0*UH(I,J,K,2)+UH(I,J,1,2)

  30  CONTINUE
      
      CALL CTDMA3I(AMK(1,1,J),ACK(1,1,J),APK(1,1,J),
     >          R3(1,1,J),UH(0,0,0,2),J,N3M,1,N1M)
      
     
  3   CONTINUE

      RETURN
      END

c***************** GETUH3 ***********************     
      Subroutine GETUH3(U,UH,RUH3,SGSVIS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL       
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS      
      COMMON/WALLCO/HPW(M2),HMW(M2),HCW(M2)
      
      REAL U(0:M1,0:M2,0:M3,3)
      REAL UH(0:M1,0:M2,0:M3,3)
      REAL RUH3(0:M1,0:M2,0:M3)
      REAL SGSVIS(0:M1,0:M2,0:M3)
      REAL NUTI,NUTJ,NUTK
      REAL API(M1,M2,M3),ACI(M1,M2,M3),AMI(M1,M2,M3)
      REAL APJ(M1,M2,M3),ACJ(M1,M2,M3),AMJ(M1,M2,M3)
      REAL APK(M1,M3,M2),ACK(M1,M3,M2),AMK(M1,M3,M2)
      REAL R1(M1,M2,M3),R2(M1,M2,M3),R3(M1,M3,M2)
C      EQUIVALENCE(API,APJ)
C      EQUIVALENCE(ACI,ACJ)
C      EQUIVALENCE(AMI,AMJ)
C      EQUIVALENCE(R1,R2)

!$omp parallel do private(JP,JM,IP,IM,KP,KM,FJUM,FJUP)     
      DO 2 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)      
 
      DO 20 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      
      DO 20 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

    
      APJ(I,J,K)=FJUP*(
     >      -0.5*HP(J)*(1./RE))*DT  
      ACJ(I,J,K)=1.0+(
     >      +0.5*HC(J)*(1./RE))*DT   
      AMJ(I,J,K)=FJUM*(
     >      -0.5*HM(J)*(1./RE))*DT
    

      R2(I,J,K)=DT*RUH3(I,J,K)

  20  CONTINUE

      CALL TDMAI(AMJ(1,1,K),ACJ(1,1,K),APJ(1,1,K),
     >         R2(1,1,K),R2(1,1,K),1,N2M,1,N1M)

      DO 21 J=1,N2M
      DO 21 I=1,N1M
      UH(I,J,K,3)=R2(I,J,K)
  21  CONTINUE    
  2   CONTINUE


!$omp parallel do private(KP,KM,IP,IM,JP,JM,FJUM,FJUP)    
      DO 3 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)      
  
      DO 30 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 30 I=1,N1M
      IP=IPA(I)
      IM=IMA(I) 
 
      APK(I,K,J)=1.0  -DT/2.0*(1./RE)*12.0*DX3Q*( 1.0)
      ACK(I,K,J)=10.0 -DT/2.0*(1./RE)*12.0*DX3Q*(-2.0)
      AMK(I,K,J)=1.0  -DT/2.0*(1./RE)*12.0*DX3Q*( 1.0)

      R3(I,K,J)=UH(I,J,K-1,3)+10.0*UH(I,J,K,3)+UH(I,J,K+1,3)      
      IF (K.EQ.1) R3(I,K,J)=UH(I,J,N3M,3)+10.0*UH(I,J,1,3)+UH(I,J,2,3)
      IF (K.EQ.N3M) R3(I,K,J)=UH(I,J,N3M-1,3)+10.0*UH(I,J,N3M,3)
     >                        +UH(I,J,1,3)

  30  CONTINUE

      CALL CTDMA3I(AMK(1,1,J),ACK(1,1,J),APK(1,1,J),
     >         R3(1,1,J),UH(0,0,0,3),J,N3M,1,N1M)

   
  3   CONTINUE


!$omp parallel do private(IP,IM,KP,KM)     
      DO 1 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
 
      DO 10 J=1,N2M
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

      API(I,J,K)=1.0  -DT/2.0*(1./RE)*12.0*DX1Q*( 1.0)
      ACI(I,J,K)=10.0 -DT/2.0*(1./RE)*12.0*DX1Q*(-2.0)
      AMI(I,J,K)=1.0  -DT/2.0*(1./RE)*12.0*DX1Q*( 1.0)

      R1(I,J,K)=UH(I-1,J,K,3)+10.0*UH(I,J,K,3)+UH(I+1,J,K,3)      
      IF (I.EQ.1) R1(I,J,K)=UH(N1M,J,K,3)+10.0*UH(I,J,K,3)+UH(I+1,J,K,3)
      IF (I.EQ.N1M) R1(I,J,K)=UH(I-1,J,K,3)+10.0*UH(I,J,K,3)+UH(1,J,K,3)

  10  CONTINUE

      CALL CTDMA1J(AMI(1,1,K),ACI(1,1,K),API(1,1,K),
     >         R1(1,1,K),R1(1,1,K),1,N1M,1,N2M)

      DO 11 J=1,N2M 
      DO 11 I=1,N1M 
      UH(I,J,K,3)=R1(I,J,K) 
  11  CONTINUE
      
     
  1   CONTINUE


C     INTERMEDIATE VELOCITY UPDATE
!$omp parallel do 
      DO 300 K=1,N3M
      DO 300 J=1,N2M
      DO 300 I=1,N1M
      UH(I,J,K,1)=U(I,J,K,1)+UH(I,J,K,1)
      UH(I,J,K,3)=U(I,J,K,3)+UH(I,J,K,3)
300   CONTINUE

!$omp parallel do 
      DO 310 K=1,N3M
      DO 310 J=2,N2M
      DO 310 I=1,N1M
      UH(I,J,K,2)=U(I,J,K,2)+UH(I,J,K,2)
310   CONTINUE

      RETURN
      END
      
c***************** TDMAI ***********************     
      SUBROUTINE TDMAI(A,B,C,D,X,NS,NF,NIS,NIF)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M1,M2),B(M1,M2),C(M1,M2),D(M1,M2),X(M1,M2)
      REAL BET,GAM(M1,0:M2)
     
      J=NS
      DO I=NIS,NIF
      BET=1.0/B(I,J)
      GAM(I,J)=C(I,J)*BET
      X(I,J)=D(I,J)*BET
      ENDDO
         
      DO 10 J=NS+1,NF
      DO I=NIS,NIF
      BET=1.0/(B(I,J)-A(I,J)*GAM(I,J-1))
      GAM(I,J)=C(I,J)*BET
      X(I,J)=(D(I,J)-A(I,J)*X(I,J-1))*BET
      ENDDO
 10   CONTINUE

      DO 20 J=NF-1,NS,-1
      DO I=NIS,NIF
      X(I,J)=X(I,J)-X(I,J+1)*GAM(I,J)
      ENDDO
 20   CONTINUE

      RETURN
      END

c***************** CTDMA1J ***********************    
      SUBROUTINE CTDMA1J(A,B,C,D,X,NS,NF,NJS,NJF)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M1,*),B(M1,*),C(M1,*),D(M1,*),X(M1,*)
      REAL BET1,BET2,GAM(0:M1,M2)
      REAL P(0:M1,M2),Q(0:M1,M2)

      DO J=NJS,NJF-1,2
      I=NS
      BET1=1.0/B(I,J)
      BET2=1.0/B(I,J+1)
      GAM(I,J)=C(I,J)*BET1
      GAM(I,J+1)=C(I,J+1)*BET2
      X(I,J)=D(I,J)*BET1 
      X(I,J+1)=D(I,J+1)*BET2 
      P(I,J)=C(NF,J)
      P(I,J+1)=C(NF,J+1)
      Q(I,J)=A(I,J)*BET1
      Q(I,J+1)=A(I,J+1)*BET2

      DO 10 I=NS+1,NF-1
      BET1=1.0/(B(I,J)-A(I,J)*GAM(I-1,J))
      BET2=1.0/(B(I,J+1)-A(I,J+1)*GAM(I-1,J+1))
      GAM(I,J)=C(I,J)*BET1
      GAM(I,J+1)=C(I,J+1)*BET2
      P(I,J)=-P(I-1,J)*GAM(I-1,J)
      P(I,J+1)=-P(I-1,J+1)*GAM(I-1,J+1)
      Q(I,J)=-A(I,J)*BET1*Q(I-1,J)
      Q(I,J+1)=-A(I,J+1)*BET2*Q(I-1,J+1)
      X(I,J)=(D(I,J)-A(I,J)*X(I-1,J))*BET1
 10   X(I,J+1)=(D(I,J+1)-A(I,J+1)*X(I-1,J+1))*BET2
      
      BET1=1.0/(B(NF-1,J)-A(NF-1,J)*GAM(NF-2,J))
      BET2=1.0/(B(NF-1,J+1)-A(NF-1,J+1)*GAM(NF-2,J+1))
      P(NF-1,J)=A(NF,J)-P(NF-2,J)*GAM(NF-2,J)
      P(NF-1,J+1)=A(NF,J+1)-P(NF-2,J+1)*GAM(NF-2,J+1)
      Q(NF-1,J)=(C(NF-1,J)-A(NF-1,J)*Q(NF-2,J))*BET1
      Q(NF-1,J+1)=(C(NF-1,J+1)-A(NF-1,J+1)*Q(NF-2,J+1))*BET2
      X(NF,J)=D(NF,J)
      X(NF,J+1)=D(NF,J+1)
      P(NF,J)=B(NF,J)
      P(NF,J+1)=B(NF,J+1)

      DO 20 I=NS,NF-1
      X(NF,J)=X(NF,J)-P(I,J)*X(I,J)
      X(NF,J+1)=X(NF,J+1)-P(I,J+1)*X(I,J+1)
      P(NF,J)=P(NF,J)-P(I,J)*Q(I,J)
 20   P(NF,J+1)=P(NF,J+1)-P(I,J+1)*Q(I,J+1)
      X(NF,J)=X(NF,J)/P(NF,J)
      X(NF,J+1)=X(NF,J+1)/P(NF,J+1)

      GAM(NF-1,J)=0.0
      GAM(NF-1,J+1)=0.0
      DO 30 I=NF-1,NS,-1
      X(I,J)=X(I,J)-GAM(I,J)*X(I+1,J)-Q(I,J)*X(NF,J)
 30   X(I,J+1)=X(I,J+1)-GAM(I,J+1)*X(I+1,J+1)-Q(I,J+1)*X(NF,J+1)
      ENDDO

      IF(MOD(NJF-NJS+1,2).EQ.1) THEN
      J=NJF
      I=NS
      BET1=1.0/B(I,J)
      GAM(I,J)=C(I,J)*BET1
      X(I,J)=D(I,J)*BET1 
      P(I,J)=C(NF,J)
      Q(I,J)=A(I,J)*BET1

      DO 40 I=NS+1,NF-1
      BET1=1.0/(B(I,J)-A(I,J)*GAM(I-1,J))
      GAM(I,J)=C(I,J)*BET1
      P(I,J)=-P(I-1,J)*GAM(I-1,J)
      Q(I,J)=-A(I,J)*BET1*Q(I-1,J)
 40   X(I,J)=(D(I,J)-A(I,J)*X(I-1,J))*BET1
      
      BET1=1.0/(B(NF-1,J)-A(NF-1,J)*GAM(NF-2,J))
      P(NF-1,J)=A(NF,J)-P(NF-2,J)*GAM(NF-2,J)
      Q(NF-1,J)=(C(NF-1,J)-A(NF-1,J)*Q(NF-2,J))*BET1
      X(NF,J)=D(NF,J)
      P(NF,J)=B(NF,J)

      DO 50 I=NS,NF-1
      X(NF,J)=X(NF,J)-P(I,J)*X(I,J)
 50   P(NF,J)=P(NF,J)-P(I,J)*Q(I,J)
      X(NF,J)=X(NF,J)/P(NF,J)

      GAM(NF-1,J)=0.0
      DO 60 I=NF-1,NS,-1
 60   X(I,J)=X(I,J)-GAM(I,J)*X(I+1,J)-Q(I,J)*X(NF,J)
      ENDIF

      RETURN
      END

c***************** CTDMA3I ***********************    
      SUBROUTINE CTDMA3I(A,B,C,D,X,J,N,NIS,NIF)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M1,M3),B(M1,M3),C(M1,M3),D(M1,M3)
      REAL X(0:M1,0:M2,0:M3)
      REAL BET,GAM(M1,M3)
      REAL P(M1,M3),Q(M1,M3)

      DO I=NIS,NIF
      BET=1./B(I,1)
      GAM(I,1)=C(I,1)*BET
      X(I,J,1)=D(I,1)*BET 
      P(I,1)=C(I,N)
      Q(I,1)=A(I,1)*BET
      ENDDO

      DO 10 K=2,N-1
      DO I=NIS,NIF
      BET=1./(B(I,K)-A(I,K)*GAM(I,K-1))
      GAM(I,K)=C(I,K)*BET
      P(I,K)=-P(I,K-1)*GAM(I,K-1)
      Q(I,K)=-A(I,K)*BET*Q(I,K-1)
      X(I,J,K)=(D(I,K)-A(I,K)*X(I,J,K-1))*BET
      ENDDO
 10   CONTINUE
      
      DO I=NIS,NIF
      BET=1./(B(I,N-1)-A(I,N-1)*GAM(I,N-2))
      P(I,N-1)=A(I,N)-P(I,N-2)*GAM(I,N-2)
      Q(I,N-1)=(C(I,N-1)-A(I,N-1)*Q(I,N-2))*BET
      
      X(I,J,N)=D(I,N)
      P(I,N)=B(I,N)
      ENDDO
      DO 20 K=1,N-1
      DO I=NIS,NIF
      X(I,J,N)=X(I,J,N)-P(I,K)*X(I,J,K)
      P(I,N)=P(I,N)-P(I,K)*Q(I,K)
      ENDDO
 20   CONTINUE
      DO I=NIS,NIF
      X(I,J,N)=X(I,J,N)/P(I,N)
      GAM(I,N-1)=0.0
      ENDDO

      DO 30 K=N-1,1,-1
      DO I=NIS,NIF
      X(I,J,K)=X(I,J,K)-GAM(I,K)*X(I,J,K+1)-Q(I,K)*X(I,J,N)
      ENDDO
 30   CONTINUE
      RETURN
      END
      
C     *************** TDMA1 *************** 
      SUBROUTINE TDMA1(A,B,C,D,X)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M1-1),B(M1-1),C(M1-1),D(M1-1),X(M1-1)
      REAL BB(M1-1),CC(M1-1),DD(M1-1)
      
      BB(1)=B(1)
      CC(1)=C(1)
      DD(1)=D(1)
      
      DO 10 I=2,M1-1
      BB(I)=B(I)-CC(I-1)*A(I)/BB(I-1)
      CC(I)=C(I)
      DD(I)=D(I)-DD(I-1)*A(I)/BB(I-1)          
  10  CONTINUE
      
      X(M1-1)=DD(M1-1)/BB(M1-1)
      DO 20 I=M1-2,1,-1
      X(I)=(DD(I)-CC(I)*X(I+1))/BB(I)
  20  CONTINUE           
      
      RETURN
      END      
      
C     *************** CTDMA11 ***************
      SUBROUTINE CTDMA1(A,B,C,D,X)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M1-1),B(M1-1),C(M1-1),D(M1-1),X(M1-1)
      REAL U(M1-1),V(M1-1),Y(M1-1),Q(M1-1)
      REAL AP(M1-1),BP(M1-1),CP(M1-1)
      
      
      DO 10 I=1,M1-1
      U(I)=0.0
      V(I)=0.0
      AP(I)=A(I)
      BP(I)=B(I)
      CP(I)=C(I)
  10  CONTINUE
      U(1) =-B(1)
      U(M1-1)= C(M1-1)
      V(1) = 1.
      V(M1-1)=-A(1)/B(1)
      AP(1)=0.0
      BP(1)= 2.*B(1)
      BP(M1-1)=B(M1-1)+A(1)*C(M1-1)/B(1)
      CP(M1-1)=0.0
      
      CALL TDMA1(AP,BP,CP,D,Y)
      CALL TDMA1(AP,BP,CP,U,Q)
      
      COE1=0.0
      COE2=0.0
      
      DO 20 I=1,M1-1
      COE1=COE1+V(I)*Y(I)
      COE2=COE2+V(I)*Q(I)
  20  CONTINUE
      
      
      DO 30 I=1,M1-1
      X(I)=Y(I)-COE1/(1.+COE2)*Q(I)
  30  CONTINUE          
      
      RETURN
      END
      


C     *************** TDMA3 *************** 
      SUBROUTINE TDMA3(A,B,C,D,X)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M3-1),B(M3-1),C(M3-1),D(M3-1),X(M3-1)
      REAL BB(M3-1),CC(M3-1),DD(M3-1)
      
      BB(1)=B(1)
      CC(1)=C(1)
      DD(1)=D(1)
      
      DO 10 K=2,M3-1
      BB(K)=B(K)-CC(K-1)*A(K)/BB(K-1)
      CC(K)=C(K)
      DD(K)=D(K)-DD(K-1)*A(K)/BB(K-1)          
  10  CONTINUE
      
      X(M3-1)=DD(M3-1)/BB(M3-1)
      DO 20 K=M3-2,1,-1
      X(K)=(DD(K)-CC(K)*X(K+1))/BB(K)
  20  CONTINUE           
      
      RETURN
      END
      
C     *************** CTDMA33 ***************
      SUBROUTINE CTDMA3(A,B,C,D,X)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M3-1),B(M3-1),C(M3-1),D(M3-1),X(M3-1)
      REAL U(M3-1),V(M3-1),Y(M3-1),Q(M3-1)
      REAL AP(M3-1),BP(M3-1),CP(M3-1)
      
      
      DO 10 K=1,M3-1
      U(K)=0.0
      V(K)=0.0
      AP(K)=A(K)
      BP(K)=B(K)
      CP(K)=C(K)
  10  CONTINUE
      U(1) =-B(1)
      U(M3-1)= C(M3-1)
      V(1) = 1.
      V(M3-1)=-A(1)/B(1)
      AP(1)=0.0
      BP(1)= 2.*B(1)
      BP(M3-1)=B(M3-1)+A(1)*C(M3-1)/B(1)
      CP(M3-1)=0.0
      
      CALL TDMA3(AP,BP,CP,D,Y)
      CALL TDMA3(AP,BP,CP,U,Q)
      
      COE1=0.0
      COE2=0.0
      
      DO 20 K=1,M3-1
      COE1=COE1+V(K)*Y(K)
      COE2=COE2+V(K)*Q(K)
  20  CONTINUE
      
      
      DO 30 K=1,M3-1
      X(K)=Y(K)-COE1/(1.+COE2)*Q(K)
  30  CONTINUE          
      
      RETURN
      END
      
      SUBROUTINE CTDMA11(A,B,C,D,X)
      INCLUDE 'PARAM.H'
      PARAMETER(M1M=M1-1)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M      
      REAL A(M1M),B(M1M),C(M1M),D(M1M),X(M1M),Z(M1M)
      REAL BETA(M1M),P(M1M),GAMA(M1M),Q(M1M)
      
      BETA(1)=B(1)
      GAMA(1)=C(1)/BETA(1)      
      DO I=2,M1M-1
      BETA(I)=B(I)-A(I)*GAMA(I-1)
      GAMA(I)=C(I)/BETA(I)      
      ENDDO
      
      P(1)=C(M1M)
      Q(1)=A(1)/BETA(1)
      DO I=2,M1M-2
      P(I)=-P(I-1)*GAMA(I-1)
      Q(I)=-A(I)*Q(I-1)/BETA(I)
      ENDDO
      P(M1M-1)=A(M1M)-P(M1M-2)*GAMA(M1M-2)
      Q(M1M-1)=(C(M1M-1)-A(M1M-1)*Q(M1M-2))/BETA(M1M-1)
      P(M1M)=B(M1M)
      DO I=1,M1M-1
      P(M1M)=P(M1M)-P(I)*Q(I)
      ENDDO
      
C     FORWARD SUBSTITUTION
      Z(1)=D(1)/BETA(1)
      DO I=2,M1M-1
      Z(I)=(D(I)-A(I)*Z(I-1))/BETA(I)
      ENDDO
      Z(M1M)=D(M1M)
      DO I=1,M1M-1
      Z(M1M)=Z(M1M)-P(I)*Z(I)
      ENDDO
      Z(M1M)=Z(M1M)/P(M1M)
      
C     BACKWARD SUBSTITUTION
      X(M1M)=Z(M1M)
      X(M1M-1)=Z(M1M-1)-Q(M1M-1)*X(M1M)
      DO I=M1M-2,1,-1
      X(I)=Z(I)-GAMA(I)*X(I+1)-Q(I)*X(M1M)
      ENDDO              
      
      RETURN
      END
      
      SUBROUTINE CTDMA33(A,B,C,D,X)
      INCLUDE 'PARAM.H'
      PARAMETER(M3M=M3-1)
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      REAL A(M3M),B(M3M),C(M3M),D(M3M),X(M3M),Z(M3M)
      REAL BETA(M3M),P(M3M),GAMA(M3M),Q(M3M)
      
      BETA(1)=B(1)
      GAMA(1)=C(1)/BETA(1)      
      DO I=2,M3M-1
      BETA(I)=B(I)-A(I)*GAMA(I-1)
      GAMA(I)=C(I)/BETA(I)      
      ENDDO
      
      P(1)=C(M3M)
      Q(1)=A(1)/BETA(1)
      DO I=2,M3M-2
      P(I)=-P(I-1)*GAMA(I-1)
      Q(I)=-A(I)*Q(I-1)/BETA(I)
      ENDDO
      P(M3M-1)=A(M3M)-P(M3M-2)*GAMA(M3M-2)
      Q(M3M-1)=(C(M3M-1)-A(M3M-1)*Q(M3M-2))/BETA(M3M-1)
      P(M3M)=B(M3M)
      DO I=1,M3M-1
      P(M3M)=P(M3M)-P(I)*Q(I)
      ENDDO
      
C     FORWARD SUBSTITUTION
      Z(1)=D(1)/BETA(1)
      DO I=2,M3M-1
      Z(I)=(D(I)-A(I)*Z(I-1))/BETA(I)
      ENDDO
      Z(M3M)=D(M3M)
      DO I=1,M3M-1
      Z(M3M)=Z(M3M)-P(I)*Z(I)
      ENDDO
      Z(M3M)=Z(M3M)/P(M3M)
      
C     BACKWARD SUBSTITUTION
      X(M3M)=Z(M3M)
      X(M3M-1)=Z(M3M-1)-Q(M3M-1)*X(M3M)
      DO I=M3M-2,1,-1
      X(I)=Z(I)-GAMA(I)*X(I+1)-Q(I)*X(M3M)
      ENDDO              
      
      RETURN
      END  

c***************** DIVCHECK ***********************    
      SUBROUTINE DIVCHECK(U,TIME,DIVMAX,NTIME) 
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX

      REAL U(0:M1,0:M2,0:M3,3)
      REAL A1(M1-1),B1(M1-1),C1(M1-1),D1(M1-1),X1(M1-1)
      REAL A3(M3-1),B3(M3-1),C3(M3-1),D3(M3-1),X3(M3-1)
      REAL DIV(M1-1,M2-1,M3-1)
      
      DIVMAX=0.0

!$omp parallel do private(KP,JP,IP) 
      DO 20 K=1,N3M
      KP=KPA(K)
      DO 20 J=1,N2M
      JP=J+1
      DO 20 I=1,N1M
      IP=IPA(I)
      DIV(I,J,K)= (U(I,JP,K,2)-U(I,J,K,2))/DY(J)
C     >            +(U(IP,J,K,1)-U(I,J,K,1))*DX1+
C     >            (U(I,J,KP,3)-U(I,J,K,3))*DX3
C      DIVMAX=AMAX1(ABS(DIV),DIVMAX)
  20  CONTINUE
    
!$OMP PARALLEL DO
      DO 15 I=1,N1M
      A1(I)=1.0
      B1(I)=22.0
      C1(I)=1.0
 15   CONTINUE      

C     DIV(U)
!$OMP PARALLEL DO PRIVATE (D1,X1)
      DO 21 K=1,N3M
      DO 21 J=1,N2M     
 
      DO I=1,N1M-1
      D1(I)=24.*DX1*(U(I+1,J,K,1)-U(I,J,K,1))
      ENDDO
      D1(N1M)=24.*DX1*(U(1,J,K,1)-U(N1M,J,K,1))
      
      CALL CTDMA11(A1,B1,C1,D1,X1)

      DO I=1,N1M
      DIV(I,J,K)=DIV(I,J,K)+X1(I)
      ENDDO      
      
  21  CONTINUE

C     DIV(W)
!$OMP PARALLEL DO
      DO 25 K=1,N3M
      A3(K)=1.0
      B3(K)=22.0
      C3(K)=1.0
  25  CONTINUE 
      
!$OMP PARALLEL DO reduction(max:DIVMAX) PRIVATE (D3,X3)
      DO 30 I=1,N1M
      DO 30 J=1,N2M     
      
      DO K=1,N3M-1
      D3(K)=24.*DX3*(U(I,J,K+1,3)-U(I,J,K,3))
      ENDDO
      D3(N3M)=24.*DX3*(U(I,J,1,3)-U(I,J,N3M,3))
      
      CALL CTDMA33(A3,B3,C3,D3,X3)

      DO K=1,N3M
      DIV(I,J,K)=DIV(I,J,K)+X3(K) 
      DIVMAX=AMAX1(ABS(DIV(I,J,K)),DIVMAX)
      ENDDO      
      
  30  CONTINUE
     
      RETURN
      END

c*********************** CFL ***********************
c     This subroutine calculate the maximum local CFL number
c     devided by DT
c     AT THE CELL CENTER
c
      SUBROUTINE CFL(U,CFLM)
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX

      REAL U(0:M1,0:M2,0:M3,3)
      CFLM=0.0
!$omp parallel do private(KP,JP,IP,CFLL) reduction(max:CFLM)
      DO 10 K=1,N3M
      KP=KPA(K)
      DO 10 J=1,N2M
      JP=J+1
      DO 10 I=1,N1M
      IP=IPA(I)
      CFLL=ABS(U(I,J,K,1)+U(IP,J,K,1))*0.5*DX1
     >     +ABS(U(I,J,K,2)+U(I,JP,K,2))*0.5/DY(J)
     >     +ABS(U(I,J,K,3)+U(I,J,KP,3))*0.5*DX3
      CFLM=AMAX1(CFLM,CFLL)
 10   CONTINUE
      
      

C      WRITE(*,*) CFLM
      
      RETURN
      END

C-JIChoi 99/10/30 X1/X3 FOURIER COMPUTATIONS FOR PERIODIC CHANNEL FLOW
C**********************************************************************
C
C     POISSON EQUATION IS SOLVED BY
C     1. FFT IN X1 & X3-DIRECTION
C     3. TRIDIAGONAL SOLVER IN X2-DIRECTION
C
C
C   ********************* INIWAVE ******************************
C
C     GET THE MODIFIED WAVENUMBERS AND THE LOCATION
      SUBROUTINE INIWAVE

      INCLUDE 'PARAM.H'
      PARAMETER (M3M=M3-1,M3MH=M3M/2+1)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/DAK/N3MH,N3MD,N3MDU,N3MU
      COMMON/WAVK13/AK1(M1),AK3(M3MH)

      REAL ANI(M1),ANK(M3)

      N3MH=N3M/2+1
      PI=ACOS(-1.0)

C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X1 DIRECTION
C     IT MUST BE CONSIDERED ABOUT FREQUENCY SHIFTING IN X1 DIRECTION
      DO 10 I=1,N1M
      II=I
      IF(I.GT.N1M/2+1) II=I-N1M
   10 ANI(I)=(II-1)*2.0*PI        ! 2*PI*m
!      WRITE(6,764) (ANI(I),I=1,N1M)
c  764 FORMAT(1X,'ANI',1X,11(F8.3,1X))

      DO 11 II=1,N1M
c   11 AK1(II)=12.*DX1Q*(1.-COS(ANI(II)/N1M))/(5.+COS(ANI(II)/N1M))
   11 AK1(II)=288.*DX1Q*(1.-COS(ANI(II)/N1M))/(11.+COS(ANI(II)/N1M))**2
!      WRITE(6,764) (AK1(I),I=1,N1M)

C     MODIFIED WAVE NUMBER DEFINITION NECESSARY FOR THE X3 DIRECTION
      DO 20 K=1,N3MH
   20 ANK(K)=(K-1)*2.*PI
!      WRITE(6,765) (ANK(K),K=1,N3MH)
c  765 FORMAT(1X,'ANK',1X,11(F8.3,1X))

      DO 21 KK=1,N3MH
c  21  AK3(KK)=12.*DX3Q*(1.-COS(ANK(KK)/N3M))/(5.+COS(ANK(KK)/N3M))
   21  AK3(KK)=288.*DX3Q*(1.-COS(ANK(KK)/N3M))/(11.+COS(ANK(KK)/N3M))**2
!      WRITE(6,765) (AK3(K),K=1,N3MH)

      CALL METRICPOISSON

      RETURN
      END
      
  
      

C  ************************ METRICPOISSON  **********************
C
C     CALCULATE THE COEFFICIENTS OF THE POISSON EQUATION FOR DPH.
C     COEFFICIENTS FOR THE THREE POINTS STENCIL OF THE POISSON EQ.
C

      SUBROUTINE METRICPOISSON
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/METPOI/PMJ(M2),PCJ(M2),PPJ(M2)

      DO 100 JC=1,N2M
      JP=JC+1
      JM=JC-1
      FJUM=FJMU(JC)
      FJUP=FJPA(JC)
      PMJ(JC)=FJUM*FJUP*(1.0/DY(JC)/H(JC))
     >       +(1.-FJUP)*(1.0/DY(N2M)/H(N2M))
      PCJ(JC)=FJUM*FJUP*(-1.0/DY(JC)/H(JC)-1.0/DY(JC)/H(JP))
     >       +(1.-FJUP)*(-1.0/DY(N2M)/H(N2M))
     >       +(1.-FJUM)*(-1.0/DY(1)/H(2))
      PPJ(JC)=FJUM*FJUP*(1.0/DY(JC)/H(JP))+(1.-FJUM)*(1.0/DY(1)/H(2))
  100 CONTINUE

      RETURN
      END


C   ********************* GETDP ********************************
C     CALCULATE DP.
C     PERIODIC DIRECTION ALONG X1 & X3 ALLOWS
C     TO USE THE REAL FOURIER TRANSFORM.

      SUBROUTINE GETDP(DP,RDP,NTIME)

      INCLUDE 'PARAM.H'
!       include 'fftw3.f'
      PARAMETER (M1M=M1-1,M3M=M3-1,M3MH=M3M/2+1,M3MD=M3M+2)
      PARAMETER (MM=M1M*M3M)

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/WAVK13/AK1(M1),AK3(M3MH)
      COMMON/METPOI/PMJ(M2),PCJ(M2),PPJ(M2)

      integer*8 fwd , bwd 
      REAL RDP(M1,M2,M3), DP(M1,M2,M3)
      REAL FRDP(M1,M2,M3MD), FDP(M1,M2,M3MD)
      REAL CMJ(M2),CCJ(M2),CPJ(M2),CFJ(M2)
      COMPLEX*16 COEF(M1M,M3M,M2)

      INTEGER fftw_forward
      PARAMETER (fftw_forward=-1)
      INTEGER fftw_backward
      PARAMETER (fftw_backward=+1)

      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)

      INTEGER fftw_estimate
      PARAMETER (fftw_estimate=64)
      
      integer iret
      
      call dfftw_init_threads(iret)
      
c      NTHREADS = 16

      nthds = omp_get_max_threads()

      call dfftw_plan_with_nthreads(nthds)
      
      call dfftw_plan_dft_2d(fwd, n1m, n3m, COEF, COEF, 
     >                       FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(bwd, n1m, n3m, COEF, COEF, 
     >                      FFTW_BACKWARD, FFTW_ESTIMATE)

!     MKL_DFT IS APPLIED TO RDP!!

      N1MH=N1M/2+1
      N3MH=N3M/2+1
      
!$omp parallel do private(KR,KI) 
      DO 10 J=1,N2M      
    
      DO 1 K=1,N3M
      DO 1 I=1,N1M
      COEF(I,K,J)=CMPLX(RDP(I,J,K),0.0)         ! MAKE COMPLEX VARIABLE
    1 CONTINUE

      call dfftw_execute_dft(fwd,coef(1,1,J),coef(1,1,J))      

!      RATIO=DBLE(N1M*N3M) ! DEC ONLY USE
!     REDUCE FOURIER COMPONENT N3M->N3MH
     
      DO 11 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 11 I=1,N1M
      FRDP(I,J,KR)=REAL(COEF(I,K,J))            
      FRDP(I,J,KI)=AIMAG(COEF(I,K,J))
   11 CONTINUE
      
   10 CONTINUE

!     SOLVE FDP BY USING TDMA IN X2-DIRECTION
!     SPLIT REAL AND IMAGINARY PART
!---- REAL PART CALCULATION

!$omp parallel do private(KR,CMJ,CCJ,CPJ,CFJ)
      DO 110 K=1,N3MH
      KR=2*K-1

      DO 110 I=1,N1M
     
      DO 120 J=1,N2M
      CMJ(J)=PMJ(J)
      CCJ(J)=PCJ(J)-AK3(K)-AK1(I)
      CPJ(J)=PPJ(J)
      CFJ(J)=FRDP(I,J,KR)
      IF(I.EQ.1.AND.K.EQ.1) THEN
      CMJ(1)=0.0
      CCJ(1)=1.0
      CPJ(1)=0.0
      CFJ(1)=0.0
      ENDIF
  120 CONTINUE
      
      CALL TDMAP(CMJ,CCJ,CPJ,CFJ,N2M,CFJ)
      
      DO 130 J=1,N2M
      FDP(I,J,KR)=CFJ(J)
  130 CONTINUE
      
  110 CONTINUE

!---- IMAG PART CALCULATION
!$omp parallel do private(KI,CMJ,CCJ,CPJ,CFJ)
      DO 140 K=1,N3MH
      KI=2*K
      DO 140 I=1,N1M
      
      DO 150 J=1,N2M
      CMJ(J)=PMJ(J)
      CCJ(J)=PCJ(J)-AK3(K)-AK1(I)
      CPJ(J)=PPJ(J)
      CFJ(J)=FRDP(I,J,KI)
!      WRITE(*,*) I,K,J,CMJ(J),CCJ(J),CPJ(J),CFJ(J)
      IF(((I.EQ.1).AND.(K.EQ.1)).OR.((I.EQ.1).AND.(K.EQ.N3MH)).OR.
     >((I.EQ.N1MH).AND.(K.EQ.1)).OR.((I.EQ.N1MH).AND.(K.EQ.N3MH))) THEN
      CMJ(J)=0.0
      CCJ(J)=1.0
      CPJ(J)=0.0
      CFJ(J)=0.0
      ENDIF
  150 CONTINUE

      CALL TDMAP(CMJ,CCJ,CPJ,CFJ,N2M,CFJ)
      
      DO 160 J=1,N2M
      FDP(I,J,KI)=CFJ(J)
  160 CONTINUE
      
  140 CONTINUE

!     DO THE INVERSE FFT.
!$omp parallel do private(KR,KI,KK,II)  
      DO 20 J=1,N2M      
    
      DO 21 K=1,N3MH
      KR=2*K-1
      KI=2*K
      DO 21 I=1,N1M
      COEF(I,K,J)=CMPLX(FDP(I,J,KR),FDP(I,J,KI))
   21 CONTINUE

      DO 22 K=N3MH+1,N3M
      KK=N3M-K+2
      DO 22 I=1,N1M
      II=N1M-I+2
      IF(I.EQ.1) II=1
      COEF(I,K,J)=CONJG(COEF(II,KK,J))
   22 CONTINUE

      call dfftw_execute_dft(bwd, COEF(1,1,J), COEF(1,1,J))

      DO 23 K=1,N3M
      DO 23 I=1,N1M
      DP(I,J,K)=REAL(COEF(I,K,J))/FLOAT(N1M*N3M)
      
   23 CONTINUE
   
   20 CONTINUE

      call dfftw_destroy_plan(FWD)
      call dfftw_destroy_plan(BWD)
      
      RETURN
      END

C  ****************************** TDMAP **********************
      SUBROUTINE TDMAP(A,B,C,R,N,X)
      INCLUDE 'PARAM.H'
      REAL GAM(M2),A(M2),B(M2),C(M2),R(M2),X(M2)

      BET=B(1)
      X(1)=R(1)/BET
      DO 11 J=2,N
      GAM(J)=C(J-1)/BET
      BET=B(J)-A(J)*GAM(J)
      X(J)=(R(J)-A(J)*X(J-1))/BET
   11 CONTINUE
      DO 12 J=N-1,1,-1
      X(J)=X(J)-GAM(J+1)*X(J+1)
   12 CONTINUE
      RETURN
      END
      

c***************** CHKMF ***********************    
      SUBROUTINE CHKMF(U,TIME)  
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/FLOWRATE/FLOW1,FLOW3

      REAL U(0:M1,0:M2,0:M3,3)

C     Calculate the mass flow rate in x1 direction      
      FLOW1=0.0 
!$omp parallel do reduction(+:FLOW1)
      DO 1 K=1,N3M
      DO 1 J=1,N2M
      DO 1 I=1,N1M
      FLOW1=FLOW1
     >       +U(1,J,K,1)*DY(J)/DX3/DX1
    1 CONTINUE    
      FLOW1=FLOW1/ALX/ALZ    ! NORMALIZED LENGTH IS Y

C     Calculate the mass flow rate in x3 direction      
      FLOW3=0.0 
!$omp parallel do reduction(+:FLOW3)
      DO 2 K=1,N3M
      DO 2 J=1,N2M
      DO 2 I=1,N1M
      FLOW3=FLOW3
     >       +U(I,J,1,3)*DY(J)/DX1/DX3
    2 CONTINUE
      FLOW3=FLOW3/ALX/ALZ
   
      WRITE(32,100) TIME,FLOW1,FLOW3
 100  FORMAT(3(E12.5,2X))

      RETURN
      END

C***************** PROFILE ***********************    
C     SPATIALLY AVERAGED FLOW FILED
C     AT A CERTAIN INSTANTANEOUS TIME

      SUBROUTINE PROFILE(U,P,NTIME)  
      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/PARA/RE
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)

      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL NXZ

      REAL PM_STEP(0:M2)
      REAL UMS_STEP(0:M2),VMS_STEP(M2),WMS_STEP(0:M2),UVMS_STEP(0:M2)
      REAL PMS_STEP(0:M2)

      REAL KAPA,B
            
      KAPA=0.4
      B=5.2
      
      OPEN(57,FILE='KAPAB.PLT')
      
      DO I=1,13500,100
      YH=I*1.0
      WRITE(57,*)I, (1.0/KAPA)*ALOG(YH)+B
      ENDDO
      
      NXZ=1.0/DBLE(N1M*N3M)

!$omp parallel do 
      DO 1 J=1,N2M
      UM_STEP(J)=0.0
      WM_STEP(J)=0.0
      PM_STEP(J)=0.0
      DO 10 K=1,N3M
      DO 10 I=1,N1M
      UM_STEP(J)=UM_STEP(J)+U(I,J,K,1)
      WM_STEP(J)=WM_STEP(J)+U(I,J,K,3)
      PM_STEP(J)=PM_STEP(J)+P(I,J,K)
   10 CONTINUE
      UM_STEP(J)=UM_STEP(J)*NXZ
      WM_STEP(J)=WM_STEP(J)*NXZ
      PM_STEP(J)=PM_STEP(J)*NXZ
    1 CONTINUE    

!$omp parallel do 
      DO 2 J=1,N2
      VM_STEP(J)=0.0
      DO 20 K=1,N3M
      DO 20 I=1,N1M
      VM_STEP(J)=VM_STEP(J)+U(I,J,K,2)
   20 CONTINUE
      VM_STEP(J)=VM_STEP(J)*NXZ
    2 CONTINUE

      UM_STEP(0)=0.0
      WM_STEP(0)=0.0
      PM_STEP(0)=PM_STEP(1)
      PM_STEP(N2)=PM_STEP(N2M)
      
!      print*,'um1 =',um_step(1),'umn2m =',um_step(n2m)
!      pause
      
      IF(NTIME.EQ.NTST)THEN
      
!$omp parallel do 
      DO 3 J=1,N2M
      UMS_STEP(J)=0.0
      WMS_STEP(J)=0.0
      PMS_STEP(J)=0.0
      DO 30 K=1,N3M
      DO 30 I=1,N1M
      UMS_STEP(J)=UMS_STEP(J)+(U(I,J,K,1)-UM_STEP(J))**2
      WMS_STEP(J)=WMS_STEP(J)+(U(I,J,K,3)-WM_STEP(J))**2
      PMS_STEP(J)=PMS_STEP(J)+(P(I,J,K)  -PM_STEP(J))**2
   30 CONTINUE
      UMS_STEP(J)=UMS_STEP(J)*NXZ
      WMS_STEP(J)=WMS_STEP(J)*NXZ
      PMS_STEP(J)=PMS_STEP(J)*NXZ
    3 CONTINUE    

!$omp parallel do 
      DO 4 J=1,N2
      VMS_STEP(J)=0.0
      DO 40 K=1,N3M
      DO 40 I=1,N1M
      VMS_STEP(J)=VMS_STEP(J)+(U(I,J,K,2)-VM_STEP(J))**2
   40 CONTINUE
      VMS_STEP(J)=VMS_STEP(J)*NXZ
    4 CONTINUE

      UMS_STEP(0)=0.0
      WMS_STEP(0)=0.0
      PMS_STEP(0)=PM_STEP(1)
      PMS_STEP(N2)=PMS_STEP(N2M)
      
      OPEN(36,FILE='U.plt',STATUS='UNKNOWN')
      OPEN(37,FILE='V.plt',STATUS='UNKNOWN')
      OPEN(38,FILE='W.plt',STATUS='UNKNOWN')
      OPEN(39,FILE='P.plt',STATUS='UNKNOWN')

      X2=0.0
      WRITE(36,200) X2,UM_STEP(0),UMS_STEP(0),U(1,0,1,1)
      WRITE(38,200) X2,WM_STEP(0),WMS_STEP(0),U(1,0,1,3)
      DO 50 J=1,N2
      X2=X2+H(J)
      WRITE(36,200) X2,UM_STEP(J),UMS_STEP(J),U(1,J,1,1)
      WRITE(38,200) X2,WM_STEP(J),WMS_STEP(J),U(1,J,1,3)
      IF (J.NE.N2) WRITE(39,200) X2,PM_STEP(J),PMS_STEP(J),P(1,J,1)
   50 CONTINUE

      DO 60 J=1,N2
      WRITE(37,200) Y(J),VM_STEP(J),VMS_STEP(J),U(1,J,1,2)
   60 CONTINUE


  100 FORMAT(2(E12.5,2X))
  200 FORMAT(4(E12.5,2X))
      
      CLOSE(36)
      CLOSE(37)
      CLOSE(38)
      CLOSE(39)

      ENDIF

      RETURN
      END

C********************** SGSTRESS **************************
C     THIS SUBROUTINE IS TO OBTAIN X-Z PLANE AVERAGED SUBGRID-SCALE STRESS AND DISSIPATION. 
C     VALUES ARE EVALUATED AT THE CELL CENTER.
 
      SUBROUTINE SGSTRESS(U,SGSVIS)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/SGSMEAN/STRM(M2,6),SGSNUT(M2),SGSDS(M2),SMAGC(M2)

      REAL U(0:M1,0:M2,0:M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL SR(M1,M3,6)
      REAL NXZ
      
      NXZ=1.0/DBLE(N1M*N3M)
      
!$omp parallel do private(STRM_1,STRM_2,STRM_3,STRM_4,
!$omp& STRM_5,STRM_6,STRMMS,SGS_NUT,COEFF,SUMSR)
      DO 10 J=1,N2M
      CALL STRAIN(U,SR,J) 
      STRM_1=0.0
      STRM_2=0.0
      STRM_3=0.0
      STRM_4=0.0
      STRM_5=0.0
      STRM_6=0.0
      STRMMS=0.0
      SGS_NUT=0.0
      DO 20 K=1,N3M
      DO 20 I=1,N1M
      COEFF=-2.0*SGSVIS(I,J,K)
      SUMSR= 2.*SR(I,K,2)**2+SR(I,K,1)**2
     >      +2.*SR(I,K,3)**2+SR(I,K,4)**2
     >      +2.*SR(I,K,5)**2+SR(I,K,6)**2
      SGS_NUT=SGS_NUT+SGSVIS(I,J,K)
      STRM_1=STRM_1+COEFF*SR(I,K,1)
      STRM_2=STRM_2+COEFF*SR(I,K,2)
      STRM_3=STRM_3+COEFF*SR(I,K,3)
      STRM_4=STRM_4+COEFF*SR(I,K,4)
      STRM_5=STRM_5+COEFF*SR(I,K,5)
      STRM_6=STRM_6+COEFF*SR(I,K,6)
      STRMMS=STRMMS+COEFF*SUMSR
   20 CONTINUE
      SGSNUT(J)=SGS_NUT*NXZ   ! x-z plane averaged sgs viscosity
      STRM(J,1)=STRM_1*NXZ    ! x-z plane averaged sgs stress tao_11
      STRM(J,2)=STRM_2*NXZ    ! x-z plane averaged sgs stress tao_12
      STRM(J,3)=STRM_3*NXZ    ! x-z plane averaged sgs stress tao_13
      STRM(J,4)=STRM_4*NXZ    ! x-z plane averaged sgs stress tao_22
      STRM(J,5)=STRM_5*NXZ    ! x-z plane averaged sgs stress tao_23
      STRM(J,6)=STRM_6*NXZ    ! x-z plane averaged sgs stress tao_33
      SGSDS(J)=STRMMS*NXZ     ! x-z plane averaged sgs dissipation rate
   10 CONTINUE

      RETURN
      END

C  ***************************** INSFIELD **********************

      SUBROUTINE INSFIELD(U,P,NAVG,SGSVIS)

      INCLUDE 'PARAM.H'
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WMEANO/WMO(M2,3),WRMSO(M2,3)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/SIZE/ALX,ALY,ALZ,VOL

      REAL U(0:M1,0:M2,0:M3,3),P(M1,M2,M3)
      REAL VOR(M1,M2,M3,3),SGSVIS(0:M1,0:M2,0:M3)
      REAL WS(M1,M3,2)

      OPEN(500,FILE='INSU_YZ.plt',STATUS='UNKNOWN')
      OPEN(501,FILE='INSOME_YZ.plt',STATUS='UNKNOWN')
      OPEN(502,FILE='INSU_XZ.plt',STATUS='UNKNOWN')
      OPEN(503,FILE='INSOME_XZ.plt',STATUS='UNKNOWN')
      OPEN(504,FILE='INSXZW.plt',STATUS='UNKNOWN')
      OPEN(505,FILE='INSNUT_YZ.plt',STATUS='UNKNOWN')
C     TO WRITE SPANWISE VELOCITY FIELD (JIN LEE)
      OPEN(506,FILE='INSV_XZ.plt',STATUS='UNKNOWN') ! JIN

      RNAVG=DBLE(NAVG)

C     CALCULATE VORTICITY at Cell Center Points
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

      U1KP=0.5*(U(I,J,KP,1)+U(IP,J,KP,1))
      U1KM=0.5*(U(I,J,KM,1)+U(IP,J,KM,1))
      U2KP=0.5*(U(I,J,KP,2)+U(I,JP,KP,2))
      U2KM=0.5*(U(I,J,KM,2)+U(I,JP,KM,2))
      U2IP=0.5*(U(IP,J,K,2)+U(IP,JP,K,2))
      U2IM=0.5*(U(IM,J,K,2)+U(IM,JP,K,2))
      U3IP=0.5*(U(IP,J,K,3)+U(IP,J,KP,3))
      U3IM=0.5*(U(IM,J,K,3)+U(IM,J,KP,3))

      U1JP=0.5*(U(I,JP,K,1)+U(IP,JP,K,1))
      U1JC=0.5*(U(I,J ,K,1)+U(IP,J ,K,1))
      U1JM=0.5*(U(I,JM,K,1)+U(IP,JM,K,1))
      U12=0.5/H(JP)*(DY(JP)*U1JC+DY(J)*U1JP)
      U11=0.5/H(J )*(DY(J)*U1JM+DY(JM)*U1JC)
      U12=U12*FJUP+(1.-FJUP)*(U(I,N2,K,1)+U(IP,N2,K,1))*0.5
      U11=U11*FJUM+(1.-FJUM)*(U(I,0 ,K,1)+U(IP,0 ,K,1))*0.5

      U3JP=0.5*(U(I,JP,K,3)+U(I,JP,KP,3))
      U3JC=0.5*(U(I,J ,K,3)+U(I,J ,KP,3))
      U3JM=0.5*(U(I,JM,K,3)+U(I,JM,KP,3))
      U32=0.5/H(JP)*(DY(JP)*U3JC+DY(J)*U3JP)
      U31=0.5/H(J )*(DY(J)*U3JM+DY(JM)*U3JC)
      U32=U32*FJUP+(1.-FJUP)*(U(I,N2,K,3)+U(I,N2,KP,3))*0.5
      U31=U31*FJUM+(1.-FJUM)*(U(I,0 ,K,3)+U(I,0 ,KP,3))*0.5

      DV3DX2=(U32-U31)/DY(J) 
      DV2DX3=(U2KP-U2KM)*0.5*DX3
      DV3DX1=(U3IP-U3IM)*0.5*DX1
      DV1DX3=(U1KP-U1KM)*0.5*DX3
      DV2DX1=(U2IP-U2IM)*0.5*DX1
      DV1DX2=(U12-U11)/DY(J) 

      VOR(I,J,K,1)=DV3DX2-DV2DX3
      VOR(I,J,K,2)=DV1DX3-DV3DX1
      VOR(I,J,K,3)=DV2DX1-DV1DX2
   10 CONTINUE

      DO 15 K=1,N3M
      DO 15 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)
      UHC=0.5*(U(I,1,K,1)+U(IP,1,K,1))     ! Bottom Wall
      UHW=0.5*(U(I,0,K,1)+U(IP,0,K,1))
      WS(I,K,1)=(UHC-UHW)/(DY(1)*0.5)

      UHC=0.5*(U(I,N2M,K,1)+U(IP,N2M,K,1)) ! Top Wall
      UHW=0.5*(U(I,N2 ,K,1)+U(IP,N2 ,K,1))
      WS(I,K,2)=(UHC-UHW)/(DY(N2M)*0.5)
   15 CONTINUE

C     YZ Plane view 
      WRITE(500,201) N2M,N3M
      WRITE(501,201) N2M,N3M
      WRITE(505,201) N2M,N3M
      NX1=N1M/4*0+1
      NX2=N1M/4*1+1
      NX3=N1M/4*2+1
      NX4=N1M/4*3+1
      DO 20 K=1,N3M
      X3=(REAL(K-1)+0.5)/DX3
      DO 20 J=1,N2M
      X2=Y(J)+DY(J)*0.5
      U1=U(NX1,J,K,1)
      U2=U(NX2,J,K,1)
      U3=U(NX3,J,K,1)
      U4=U(NX4,J,K,1)
      VM1=VMO(J,1)/RNAVG
      WM1=WMO(J,1)/RNAVG
      WRITE(500,202) X3,X2,U1-VM1,U2-VM1,U3-VM1,U4-VM1
      WRITE(501,202) X3,X2,VOR(NX1,J,K,1)-WM1,VOR(NX2,J,K,1)-WM1
     >                   ,VOR(NX3,J,K,1)-WM1,VOR(NX4,J,K,1)-WM1
      WRITE(505,202) X3,X2,SGSVIS(NX1,J,K),SGSVIS(NX2,J,K)
     >                    ,SGSVIS(NX3,J,K),SGSVIS(NX4,J,K)
   20 CONTINUE   

C     XZ Plane view 
      WRITE(502,211) N1M,N3M
      WRITE(503,211) N1M,N3M
      WRITE(504,211) N1M,N3M
      WRITE(506,211) N1M,N3M ! JIN
      NY1=10   ! y^+=5
      NY2=21   ! y^+=
      NY3=40   ! y^+=
      NY4=65  ! y^+=5
      VM1=VMO(NY1,1)/RNAVG
      VM2=VMO(NY2,1)/RNAVG
      VM3=VMO(NY3,1)/RNAVG
      VM4=VMO(NY4,1)/RNAVG
      WM1=WMO(NY1,1)/RNAVG
      WM2=WMO(NY2,1)/RNAVG
      WM3=WMO(NY3,1)/RNAVG
      WM4=WMO(NY4,1)/RNAVG
      PM1=PMO(1)/RNAVG
      PM2=PMO(N2)/RNAVG
      VVM1=VMO(NY1,2)/RNAVG ! JIN
      VVM2=VMO(NY2,2)/RNAVG ! JIN
      VVM3=VMO(NY3,2)/RNAVG ! JIN
      VVM4=VMO(NY4,2)/RNAVG ! JIN
      DO 30 K=1,N3M
      X3=(REAL(K-1)+0.5)/DX3
      DO 30 I=1,N1M
      IP=IPA(I)
      X1=(REAL(I-1)+0.5)/DX1
      U1=0.5*(U(I,NY1,K,1)+U(IP,NY1,K,1))
      U2=0.5*(U(I,NY2,K,1)+U(IP,NY2,K,1))
      U3=0.5*(U(I,NY3,K,1)+U(IP,NY3,K,1))
      U4=0.5*(U(I,NY4,K,1)+U(IP,NY4,K,1))

      V1=U(I,NY1,K,2) ! JIN
      V2=U(I,NY2,K,2) ! JIN
      V3=U(I,NY3,K,2) ! JIN
      V4=U(I,NY4,K,2) ! JIN

      WRITE(502,212) X1,X3,U1-VM1,U2-VM2
     >                    ,U3-VM3,U4-VM4
      WRITE(503,212) X1,X3,VOR(I,NY1,K,1)-WM1,VOR(I,NY2,K,1)-WM2
     >                   ,VOR(I,NY3,K,1)-WM3,VOR(I,NY4,K,1)-WM4
      WRITE(504,232) X1,X3,WS(I,K,1),WS(I,K,2)
     >                   ,P(I,1,K)-PM1,P(I,N2M,K)-PM2
C     JIN
      WRITE(506,212) X1,X3,V1-VVM1,V2-VVM2
     >                    ,V3-VVM3,V4-VVM4
   30 CONTINUE   

  201 FORMAT('ZONE I=',I4,'J=',I4,'F=POINT')
  202 FORMAT(6(E12.4,2X))
  211 FORMAT('ZONE I=',I4,'J=',I4,'F=POINT')
  212 FORMAT(6(E12.4,2X))
  221 FORMAT('ZONE I=',I4,'J=',I4,'F=POINT')
  222 FORMAT(6(E12.4,2X))
  231 FORMAT('ZONE I=',I4,'J=',I4,'F=POINT')
  232 FORMAT(6(E12.4,2X))

      CLOSE(500)
      CLOSE(501)
      CLOSE(502)
      CLOSE(503)
      CLOSE(504)
      CLOSE(505)
      CLOSE(506)

      RETURN
      END

C  ****************************** ENERGY **********************
C     CALCULATE TOTAL ENERGY = (V1**2+V2**2+V3**2)*0.5.
C     CALCULATE MEAN AND RMS VALUES AT THE CENTER OF CELL.

      SUBROUTINE ENERGY(U,P,NAVG)

      INCLUDE 'PARAM.H'
      
      PARAMETER (M1M=M1-1,M3M=M3-1)
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/PARA/RE
      COMMON/SIZE/ALX,ALY,ALZ,VOL
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH2/Y(0:M2)
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/TSTEP/NTST,DT,CFLMAX
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/FINDX2/FJPA(M2),FJMU(M2),FJMV(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/WMEAN/WM(M2,3),WRMS(M2,3)
      COMMON/PMEAN/PM(M2),PRMS(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      
      COMMON/SPEC/XSPCTR(M1,M2,3),ZSPCTR(M3,M2,3)
      COMMON/SPECO/XSPCTRO(M1,M2,3),ZSPCTRO(M3,M2,3)
      COMMON/SPECTAU/SPECTAUW(M1,M3,2),SPECTAUB(M1,M3,2)
      COMMON/SPECTAUO/SPECTAUWO(M1,M3,2),SPECTAUBO(M1,M3,2)
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/CORR/RXX(M1,M2,3),RZZ(M3,M2,3)
      COMMON/CORRO/RXXO(M1,M2,3),RZZO(M3,M2,3)
      COMMON/CORR2/RXY(M1,M2,3)
      COMMON/CORR2O/RXYO(M1,M2,3)
      COMMON/CORR3/RXYZ(M1,M2,M3,3)
      COMMON/CORR3O/RXYZO(M1,M2,M3,3)
      COMMON/QUAD/Q(M2,4),QP(M2,4)     
      COMMON/QUADO/QO(M2,4),QPO(M2,4)      
      
c-----for turbulent kinetic budgets 
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PVMEAN/PV(M2)
      COMMON/VHIGH/VSKEW(M2,3),VFLAT(M2,3),U2V(M2),VW2(M2)
      COMMON/PSTR/PDUDX(M2),PDVDY(M2),PDWDZ(M2)
      COMMON/DISSU/DUDX2(M2),DUDY2(M2),DUDZ2(M2)
      COMMON/DISSV/DVDX2(M2),DVDY2(M2),DVDZ2(M2)
      COMMON/DISSW/DWDX2(M2),DWDY2(M2),DWDZ2(M2) 
      COMMON/WALLTAO/TAOW
      COMMON/WALLTAOAVE/TAOWO


      REAL U(0:M1,0:M2,0:M3,3)
      REAL P(M1,M2,M3)
      REAL VOR(M1,M2,M3,3)
      REAL VP(M1,M2,M3,3)
      REAL VMP(M2,3)
      REAL NXZ,PI,FS,FS1,TAUDOWN
      REAL UC1(M1M),VC1(M1M),WC1(M1M)
      REAL UC3(M3M),VC3(M3M),WC3(M3M)
      REAL WAVE1(M1M),WAVE3(M3M),TAUC(M1M,M3M)
      REAL PYY1(M1M),PYY3(M3M),PYY(M1M,M3M)
      INTEGER iret
      COMPLEX*16 COEF1(M1M),COEF3(M3M),COEF(M1M,M3M)
      INTEGER*8 FWD1,FWD3,FWD
      
      INTEGER fftw_forward
      PARAMETER (fftw_forward=-1)
      INTEGER fftw_backward
      PARAMETER (fftw_backward=+1)

      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)

      INTEGER fftw_estimate
      PARAMETER (fftw_estimate=64)
      
      nthds = omp_get_max_threads()
      
      call dfftw_init_threads(iret)      
      call dfftw_plan_with_nthreads(nthds) 
      
      call dfftw_plan_dft_1d(FWD1, N1M, COEF1, COEF1, 
     >                       FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_1d(FWD3, N3M, COEF3, COEF3, 
     >                       FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_plan_dft_2d(FWD, N1M, N3M, COEF, COEF, 
     >                       FFTW_FORWARD, FFTW_ESTIMATE)

      PI=ACOS(-1.)      
      NXZ=1.0/DBLE(N1M*N3M)
      RNAVG=DBLE(NAVG)

      IF(WMODEL.EQ.0)THEN
      UTAU=SQRT(WSMO(1))
      ELSE
      UTAU=SQRT(TAOWO)
      ENDIF      

C     CALCULATE VORTICITY at Cell Center Points
!$omp parallel do private(KP,KM,JP,JM,FJUM,FJUP,IP,IM,
!$omp& U1KP,U1KM,U2KP,U2KM,U2IP,U2IM,U3IP,U3IM,U1JP,U1JC,U1JM,
!$omp& U12,U11,U3JP,U3JC,U3JM,U31,U32,DV3DX2,DV2DX3,DV3DX1,
!$omp& DV1DX3,DV2DX1,DV1DX2)
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

      U1KP=0.5*(U(I,J,KP,1)+U(IP,J,KP,1))
      U1KM=0.5*(U(I,J,KM,1)+U(IP,J,KM,1))
      U2KP=0.5*(U(I,J,KP,2)+U(I,JP,KP,2))
      U2KM=0.5*(U(I,J,KM,2)+U(I,JP,KM,2))
      U2IP=0.5*(U(IP,J,K,2)+U(IP,JP,K,2))
      U2IM=0.5*(U(IM,J,K,2)+U(IM,JP,K,2))
      U3IP=0.5*(U(IP,J,K,3)+U(IP,J,KP,3))
      U3IM=0.5*(U(IM,J,K,3)+U(IM,J,KP,3))

      U1JP=0.5*(U(I,JP,K,1)+U(IP,JP,K,1))
      U1JC=0.5*(U(I,J ,K,1)+U(IP,J ,K,1))
      U1JM=0.5*(U(I,JM,K,1)+U(IP,JM,K,1))
      U12=0.5/H(JP)*(DY(JP)*U1JC+DY(J )*U1JP)
      U11=0.5/H(J )*(DY(J )*U1JM+DY(JM)*U1JC)
      U12=U12*FJUP+(1.-FJUP)*(U(I,N2,K,1)+U(IP,N2,K,1))*0.5
      U11=U11*FJUM+(1.-FJUM)*(U(I,0 ,K,1)+U(IP,0 ,K,1))*0.5

      U3JP=0.5*(U(I,JP,K,3)+U(I,JP,KP,3))
      U3JC=0.5*(U(I,J ,K,3)+U(I,J ,KP,3))
      U3JM=0.5*(U(I,JM,K,3)+U(I,JM,KP,3))
      U32=0.5/H(JP)*(DY(JP)*U3JC+DY(J )*U3JP)
      U31=0.5/H(J )*(DY(J )*U3JM+DY(JM)*U3JC)
      U32=U32*FJUP+(1.-FJUP)*(U(I,N2,K,3)+U(I,N2,KP,3))*0.5
      U31=U31*FJUM+(1.-FJUM)*(U(I,0 ,K,3)+U(I,0 ,KP,3))*0.5

      DV3DX2=(U32-U31)/DY(J) 
      DV2DX3=(U2KP-U2KM)*0.5*DX3
      DV3DX1=(U3IP-U3IM)*0.5*DX1
      DV1DX3=(U1KP-U1KM)*0.5*DX3
      DV2DX1=(U2IP-U2IM)*0.5*DX1
      DV1DX2=(U12-U11)/DY(J) 

      VOR(I,J,K,1)=DV3DX2-DV2DX3
      VOR(I,J,K,2)=DV1DX3-DV3DX1
      VOR(I,J,K,3)=DV2DX1-DV1DX2
      
   10 CONTINUE

C     VORTICITY Tilting and Stretching must be considered!!!

C     CALCULATE MEAN,Mean Square and Higher-order values at Cell Center Points.

!$omp parallel do private(JP,KP,IP,V1,V2,V3,
!$omp& V1M,V2M,V3M,PMM,VM1M,VM2M,VM3M,U11MS,U22MS,U33MS,U12MS,
!$omp& U13MS,U23MS,PMMMS,PVMS,VOR1MS,VOR2MS,VOR3MS,U1SK,U2SK,U3SK,
!$omp& U1FL,U2FL,U3FL,U11U2,U2U33)      
      DO 20 J=1,N2M
      JP=J+1
      
      V1M=0.0
      V2M=0.0
      V3M=0.0
      PMM=0.0
      VM1M=0.0
      VM2M=0.0
      VM3M=0.0
      U11MS=0.0
      U22MS=0.0
      U33MS=0.0
      U12MS=0.0
      U13MS=0.0
      U23MS=0.0
      PMMMS=0.0
      PVMS=0.0
      VOR1MS=0.0
      VOR2MS=0.0
      VOR3MS=0.0
c      U1SK=0.0
c      U2SK=0.0
c      U3SK=0.0
c      U1FL=0.0
c      U2FL=0.0
c      U3FL=0.0
      U11U2=0.0
      U2U33=0.0
      
      DO 30 K=1,N3M
      KP=KPA(K)
      DO 30 I=1,N1M
      IP=IPA(I)
      V1=(U(I,J,K,1)+U(IP,J,K,1))*0.5
      V2=(U(I,J,K,2)+U(I,JP,K,2))*0.5
      V3=(U(I,J,K,3)+U(I,J,KP,3))*0.5
      V1M=V1M+V1  
      V2M=V2M+V2  
      V3M=V3M+V3  
      PMM=PMM+P(I,J,K)
      VM1M=VM1M+VOR(I,J,K,1)  
      VM2M=VM2M+VOR(I,J,K,2)  
      VM3M=VM3M+VOR(I,J,K,3)  
      U11MS=U11MS+V1**2  
      U22MS=U22MS+V2**2  
      U33MS=U33MS+V3**2  
      U12MS=U12MS+V1*V2  
      U13MS=U13MS+V1*V3
      U23MS=U23MS+V2*V3  
      PMMMS=PMMMS+P(I,J,K)**2
      PVMS=PVMS+P(I,J,K)*V2
      VOR1MS=VOR1MS+VOR(I,J,K,1)**2
      VOR2MS=VOR2MS+VOR(I,J,K,2)**2
      VOR3MS=VOR3MS+VOR(I,J,K,3)**2
c      U1SK=U1SK+V1**3
c      U2SK=U2SK+V2**3
c      U3SK=U3SK+V3**3
c      U1FL=U1FL+V1**4
c      U2FL=U2FL+V2**4
c      U3FL=U3FL+V3**4
      U11U2=U11U2+V1*V1*V2
      U2U33=U2U33+V2*V3**2
   30 CONTINUE
      
      VM(J,1)=V1M*NXZ         ! average u
      VM(J,2)=V2M*NXZ         ! average v
      VM(J,3)=V3M*NXZ         ! average w
      PM(J)=PMM*NXZ           ! average p   
      WM(J,1)=VM1M*NXZ        ! average vort1
      WM(J,2)=VM2M*NXZ        ! average vort2
      WM(J,3)=VM3M*NXZ        ! average vort3
      VRMS(J,1)=U11MS*NXZ     ! average u*u
      VRMS(J,2)=U22MS*NXZ     ! average v*v
      VRMS(J,3)=U33MS*NXZ     ! average w*w
      VRMS(J,4)=U12MS*NXZ     ! average u*v
      VRMS(J,5)=U13MS*NXZ     ! average u*w
      VRMS(J,6)=U23MS*NXZ     ! average v*w
      PRMS(J)=PMMMS*NXZ       ! average p*p
      PV(J)=PVMS*NXZ          ! average p*v
      WRMS(J,1)=VOR1MS*NXZ    ! average vort1*vort1
      WRMS(J,2)=VOR2MS*NXZ    ! average vort2*vort2
      WRMS(J,3)=VOR3MS*NXZ    ! average vort3*vort3
c      VSKEW(J,1)=U1SK*NXZ
c      VSKEW(J,2)=U2SK*NXZ
c      VSKEW(J,3)=U3SK*NXZ
c      VFLAT(J,1)=U1FL*NXZ
c      VFLAT(J,2)=U2FL*NXZ
c      VFLAT(J,3)=U3FL*NXZ
      U2V(J)=U11U2*NXZ
      VW2(J)=U2U33*NXZ
      
   20 CONTINUE
      
cccc  skewness and flatness      
!$omp parallel do private(JP,KP,IP,V1,V2,V3,
!$omp& U1SK,U2SK,U3SK,U1FL,U2FL,U3FL)         
      DO 21 J=1,N2M
      JP=J+1      
      
      U1SK=0.0
      U2SK=0.0
      U3SK=0.0
      U1FL=0.0
      U2FL=0.0
      U3FL=0.0
      
      DO 31 K=1,N3M
      KP=KPA(K)
      DO 31 I=1,N1M
      IP=IPA(I)
      
      V1=(U(I,J,K,1)+U(IP,J,K,1))*0.5-VM(J,1)
      V2=(U(I,J,K,2)+U(I,JP,K,2))*0.5-VM(J,2)
      V3=(U(I,J,K,3)+U(I,J,KP,3))*0.5-VM(J,3)      
      
      U1SK=U1SK+V1**3
      U2SK=U2SK+V2**3
      U3SK=U3SK+V3**3
      U1FL=U1FL+V1**4
      U2FL=U2FL+V2**4
      U3FL=U3FL+V3**4

   31 CONTINUE      
      
      VSKEW(J,1)=U1SK*NXZ
      VSKEW(J,2)=U2SK*NXZ
      VSKEW(J,3)=U3SK*NXZ
      VFLAT(J,1)=U1FL*NXZ
      VFLAT(J,2)=U2FL*NXZ
      VFLAT(J,3)=U3FL*NXZ
      
   21 CONTINUE

C     CALCULATE DU/DY,DV/DY,DW/DY
!$omp parallel do private(JP,JM,FJUM,FJUP,KP,IP,U_JP,U_JC,U_JM,
!$omp& W_JP,W_JC,W_JM,U1,U2,V1,V2,W1,W2,DU_DY,DV_DY,DW_DY)
      DO 40 J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J)
      DU_DY=0.0
      DV_DY=0.0
      DW_DY=0.0
      DO 50 K=1,N3M
      KP=KPA(K)
      DO 50 I=1,N1M
      IP=IPA(I)
      
      U_JP=(U(I,JP,K,1)+U(IP,JP,K,1))*0.5
      U_JC=(U(I,J ,K,1)+U(IP,J ,K,1))*0.5
      U_JM=(U(I,JM,K,1)+U(IP,JM,K,1))*0.5
      W_JP=(U(I,JP,K,3)+U(I,JP,KP,3))*0.5
      W_JC=(U(I,J ,K,3)+U(I,J ,KP,3))*0.5
      W_JM=(U(I,JM,K,3)+U(I,JM,KP,3))*0.5
      U1 = FJUM*(1.0/H(J)*(DY(J)*U_JM/2.0+DY(JM)*U_JC/2.0))
     >    +(1.-FJUM)*(U(I,0,K,1)+U(IP,0,K,1))*0.5
      U2 = FJUP*(1.0/H(JP)*(DY(J)*U_JP/2.0+DY(JP)*U_JC/2.0))
     >    +(1.-FJUP)*(U(I,N2,K,1)+U(IP,N2,K,1))*0.5
      V1=U(I,J,K,2)
      V2=U(I,JP,K,2)
      W1 = FJUM*(1.0/H(J)*(DY(J)*W_JM/2.0+DY(JM)*W_JC/2.0))
     >    +(1.-FJUM)*(U(I,0,K,3)+U(I,0,KP,3))*0.5
      W2 = FJUP*(1.0/H(JP)*(DY(J)*W_JP/2.0+DY(JP)*W_JC/2.0))
     >    +(1.-FJUP)*(U(I,N2,K,3)+U(I,N2,KP,3))*0.5
      DU_DY=DU_DY+(U2-U1)/DY(J)
      DV_DY=DV_DY+(V2-V1)/DY(J)
      DW_DY=DW_DY+(W2-W1)/DY(J)
   50 CONTINUE 
      
      DUDY(J)=DU_DY*NXZ   ! average du/dy
      DVDY(J)=DV_DY*NXZ   ! average dv/dy
      DWDY(J)=DW_DY*NXZ   ! average dw/dy
   40 CONTINUE
      
      
CCCCCCCCCCCCCCCCC 1D SEPCTRA CCCCCCCCCCCCCCCC
      
C     CALCULATE X-SPECTRUM 
      FS=PI*DX1
      
!$omp parallel do      
      DO J=1,N2M
      DO I=1,N1M
      DO L=1,3
      XSPCTR(I,J,L)=0.0
      ENDDO
      ENDDO
      ENDDO
      
!$omp parallel do private(JP,JM,FJUM,FJUP,
!$omp& KP,IP,UC1,VC1,WC1,COEF1,PYY1)
      DO J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J) 
      
      DO K=1,N3M
      KP=KPA(K)  
      
     
C     U SPECTRUM IN X DIRECTION
      DO I=1,N1M
      IP=IPA(I)      
      UC1(I)=(U(I,J,K,1)+U(IP,J,K,1))*0.5-VMO(J,1)
      COEF1(I)=CMPLX(UC1(I),0.0)      
      ENDDO
      
      call dfftw_execute_dft(FWD1,COEF1,COEF1)      
      
      DO I=1,N1M
      PYY1(I)=COEF1(I)*CONJG(COEF1(I))/N1M/FS
      XSPCTR(I,J,1)=XSPCTR(I,J,1)+PYY1(I)
      ENDDO
      
C     V SPECTRUM IN X DIRECTION      
      DO I=1,N1M
      VC1(I)=(U(I,J,K,2)+U(I,JP,K,2))*0.5-VMO(J,2)
      COEF1(I)=CMPLX(VC1(I),0.0)      
      ENDDO
      
      call dfftw_execute_dft(FWD1,COEF1,COEF1)      
      
      DO I=1,N1M
      PYY1(I)=COEF1(I)*CONJG(COEF1(I))/N1M/FS
      XSPCTR(I,J,2)=XSPCTR(I,J,2)+PYY1(I)
      ENDDO
      
C     W SPECTRUM IN X DIRECTION      
      DO I=1,N1M
      WC1(I)=(U(I,J,K,3)+U(I,J,KP,3))*0.5-VMO(J,3)
      COEF1(I)=CMPLX(WC1(I),0.0)      
      ENDDO
      
      call dfftw_execute_dft(FWD1,COEF1,COEF1)      
      
      DO I=1,N1M
      PYY1(I)=COEF1(I)*CONJG(COEF1(I))/N1M/FS
      XSPCTR(I,J,3)=XSPCTR(I,J,3)+PYY1(I)
      ENDDO      
      
      ENDDO
     
C     AVERAGE
      DO I=1,N1M
      XSPCTR(I,J,1)=XSPCTR(I,J,1)/REAL(N3M)
      XSPCTR(I,J,2)=XSPCTR(I,J,2)/REAL(N3M)
      XSPCTR(I,J,3)=XSPCTR(I,J,3)/REAL(N3M)
      ENDDO
          
      ENDDO
      
      
C-----CALCULATE Z-SPECTRUM 
      FS=PI*DX3
      
!$omp parallel do      
      DO J=1,N2M
      DO K=1,N3M
      DO L=1,3
      ZSPCTR(K,J,L)=0.0
      ENDDO
      ENDDO
      ENDDO
      
!$omp parallel do private(JP,JM,FJUM,FJUP,KP,IP,
!$omp& UC3,VC3,WC3,COEF3,PYY3)
      DO J=1,N2M
      JP=J+1
      JM=J-1
      FJUM=FJMU(J)
      FJUP=FJPA(J) 
      
      DO I=1,N1M
      IP=IPA(I)  
      
C     U SPECTRUM IN Z DIRECTION
      DO K=1,N3M            
      UC3(K)=(U(I,J,K,1)+U(IP,J,K,1))*0.5-VMO(J,1)
      COEF3(K)=CMPLX(UC3(K),0.0)      
      ENDDO
      
      call dfftw_execute_dft(FWD3,COEF3,COEF3)      
      
      DO K=1,N3M
      PYY3(K)=COEF3(K)*CONJG(COEF3(K))/N3M/FS
      ZSPCTR(K,J,1)=ZSPCTR(K,J,1)+PYY3(K)
      ENDDO
      
C     V SPECTRUM IN Z DIRECTION      
      DO K=1,N3M
      VC3(K)=(U(I,J,K,2)+U(I,JP,K,2))*0.5-VMO(J,2)
      COEF3(K)=CMPLX(VC3(K),0.0)      
      ENDDO
      
      call dfftw_execute_dft(FWD3,COEF3,COEF3)      
      
      DO K=1,N3M
      PYY3(K)=COEF3(K)*CONJG(COEF3(K))/N3M/FS
      ZSPCTR(K,J,2)=ZSPCTR(K,J,2)+PYY3(K)
      ENDDO
      
C     W SPECTRUM IN Z DIRECTION      
      DO K=1,N3M
      KP=KPA(K)
      WC3(K)=(U(I,J,K,3)+U(I,J,KP,3))*0.5-VMO(J,3)
      COEF3(K)=CMPLX(WC3(K),0.0)      
      ENDDO
      
      call dfftw_execute_dft(FWD3,COEF3,COEF3)      
      
      DO K=1,N3M
      PYY3(K)=COEF3(K)*CONJG(COEF3(K))/N3M/FS
      ZSPCTR(K,J,3)=ZSPCTR(K,J,3)+PYY3(K)
      ENDDO      
      
      ENDDO
      
      DO K=1,N3M
      ZSPCTR(K,J,1)=ZSPCTR(K,J,1)/REAL(N1M)
      ZSPCTR(K,J,2)=ZSPCTR(K,J,2)/REAL(N1M)
      ZSPCTR(K,J,3)=ZSPCTR(K,J,3)/REAL(N1M)
      ENDDO
          
      ENDDO
      
CCCCCCCCCCCCCCCC 2D energy spectra of the streamwise wall shearstress CCCCCCCCCCCCCCC  
      FS  = PI*DX1
      FS1 = PI*DX3      
      
      TAUDOWN=TAUDOWN/DBLE(N1M*N3M)
      DO K=1,N3M
      DO I=1,N1M
      TAUC(I,K)=TAUW(I,K,1)-UTAU**2
      COEF(I,K)=CMPLX(TAUC(I,K),0.0)
      ENDDO
      ENDDO
      
      call dfftw_execute_dft(FWD,COEF,COEF)
      
      DO K=1,N3M
      DO I=1,N1M
      PYY(I,K)=COEF(I,K)*CONJG(COEF(I,K))/N1M/N3M/FS1/FS
      SPECTAUW(I,K,1)=PYY(I,K)
      ENDDO    
      ENDDO      
      
      
      
CCCCCCCCCCCCCCCCC 1D TWO POINT CORRELATION CCCCCCCCCCCCCCCC

      N1H=N1M/2+1      
      N3H=N3M/2+1
      
      IF (ALY.EQ.1.0) THEN
      N2H=N2M/2+1
      ELSE
      N2H=N2M/4+1
      ENDIF

!$omp parallel do PRIVATE (JP,KP,IP,V1,V2,V3)
      DO 130 J=1,N2M
      JP=JPA(J)
      DO 130 K=1,N3M
      KP=KPA(K)
      DO 130 I=1,N1M
      IP=IPA(I)
      V1=(U(I,J,K,1)+U(IP,J, K,1))*0.5
      V2=(U(I,J,K,2)+U(I ,JP,K,2))*0.5
      V3=(U(I,J,K,3)+U(I ,J,KP,3))*0.5
      VP(I,J,K,1)=V1-VM(J,1)
      VP(I,J,K,2)=V2-VM(J,2)
      VP(I,J,K,3)=V3-VM(J,3)
  130 CONTINUE

C-----STREAMWISE TWO-POINT CORRELATION-----
!$omp parallel do
      DO J=1,N2M
      DO I=1,N1M
      DO L=1,3
      RXX(I,J,L)=0.0
      ENDDO
      ENDDO
      ENDDO      

      II=N1H
!$omp parallel do 
      DO J=1,N2M      
      DO I=1,N1M
      DO K=1,N3M
      RXX(I,J,1)=RXX(I,J,1)+VP(II,J,K,1)*VP(I,J,K,1)
      RXX(I,J,2)=RXX(I,J,2)+VP(II,J,K,2)*VP(I,J,K,2)
      RXX(I,J,3)=RXX(I,J,3)+VP(II,J,K,3)*VP(I,J,K,3)
      ENDDO
      RXX(I,J,1)=RXX(I,J,1)/DBLE(N3M)/(VRMS(J,1)-VM(J,1)*VM(J,1))
      RXX(I,J,2)=RXX(I,J,2)/DBLE(N3M)/(VRMS(J,2)-VM(J,2)*VM(J,2))
      RXX(I,J,3)=RXX(I,J,3)/DBLE(N3M)/(VRMS(J,3)-VM(J,3)*VM(J,3))
      ENDDO
      ENDDO      
      

C-----SPANWISE TWO-POINT CORRELATION-----
!$omp parallel do
      DO J=1,N2M
      DO K=1,N3M
      DO L=1,3
      RZZ(K,J,L)=0.0
      ENDDO
      ENDDO
      ENDDO      
 
      KK=N3H
!$omp parallel do
      DO J=1,N2M      
      DO K=1,N3M
      DO I=1,N1M      
      RZZ(K,J,1)=RZZ(K,J,1)+VP(I,J,KK,1)*VP(I,J,K,1)
      RZZ(K,J,2)=RZZ(K,J,2)+VP(I,J,KK,2)*VP(I,J,K,2)
      RZZ(K,J,3)=RZZ(K,J,3)+VP(I,J,KK,3)*VP(I,J,K,3)
      ENDDO
      RZZ(K,J,1)=RZZ(K,J,1)/DBLE(N1M)/(VRMS(J,1)-VM(J,1)*VM(J,1))
      RZZ(K,J,2)=RZZ(K,J,2)/DBLE(N1M)/(VRMS(J,2)-VM(J,2)*VM(J,2))
      RZZ(K,J,3)=RZZ(K,J,3)/DBLE(N1M)/(VRMS(J,3)-VM(J,3)*VM(J,3))
      ENDDO
      ENDDO    
      
CCCCCCCCCC TWO-POINT CORRELATION IN X-Y PLANE CCCCCCCCCC
      
!$omp parallel do
      DO J=1,N2M
      DO I=1,N1M
      DO L=1,3
      RXY(I,J,L)=0.0
      ENDDO
      ENDDO
      ENDDO
      
      JJ=1
      II=N1H
!$omp parallel do
      DO J=1,N2M      
      DO I=1,N1M
      DO K=1,N3M
      RXY(I,J,1)=RXY(I,J,1)+VP(II,JJ,K,1)*VP(I,J,K,1)
      RXY(I,J,2)=RXY(I,J,2)+VP(II,JJ,K,2)*VP(I,J,K,2)
      RXY(I,J,3)=RXY(I,J,3)+VP(II,JJ,K,3)*VP(I,J,K,3)
      ENDDO
      RXY(I,J,1)=RXY(I,J,1)/DBLE(N3M)/SQRT(VRMS(J,1)-VM(J,1)*VM(J,1))
     >            /SQRT(VRMS(JJ,1)-VM(JJ,1)*VM(JJ,1))
      RXY(I,J,2)=RXY(I,J,2)/DBLE(N3M)/SQRT(VRMS(J,2)-VM(J,2)*VM(J,2))
     >            /SQRT(VRMS(JJ,2)-VM(JJ,2)*VM(JJ,2))
      RXY(I,J,3)=RXY(I,J,3)/DBLE(N3M)/SQRT(VRMS(J,3)-VM(J,3)*VM(J,3))
     >            /SQRT(VRMS(JJ,3)-VM(JJ,3)*VM(JJ,3))
      ENDDO
      ENDDO  
      
CCCCCCCCCC 3-D TWO-POINT CORRELATION CCCCCCCCCC
!$omp parallel do
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      DO L=1,3
      RXYZ(I,J,K,L)=0.0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      JJ=N2H
      II=N1H
      KK=N3H
!$omp parallel do
      DO J=1,N2M      
      DO I=1,N1M
      DO K=1,N3M
      RXYZ(I,J,K,1)=VP(II,JJ,KK,1)*VP(I,J,K,1)
     >            /SQRT(VRMS(J,1)-VM(J,1)*VM(J,1))
     >            /SQRT(VRMS(JJ,1)-VM(JJ,1)*VM(JJ,1))
      RXYZ(I,J,K,2)=VP(II,JJ,KK,2)*VP(I,J,K,2)
     >            /SQRT(VRMS(J,2)-VM(J,2)*VM(J,2))
     >            /SQRT(VRMS(JJ,2)-VM(JJ,2)*VM(JJ,2)) 
      RXYZ(I,J,K,3)=VP(II,JJ,KK,3)*VP(I,J,K,3)
     >            /SQRT(VRMS(J,3)-VM(J,3)*VM(J,3))
     >            /SQRT(VRMS(JJ,3)-VM(JJ,3)*VM(JJ,3))
      ENDDO
      ENDDO
      ENDDO    
      
      
CCCCCCCCCCCCCCCCC QUADRANT ANALYSIS CCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     CALCULATE REYNOLDS SHEAR STRESS FOR EACH QUADRANT
C     NORMALIZED BY THE MEAN REYNOLDS SHEAR STRESS and TURBULENT HEAT FLUX.

      DO 230 J=1,N2M
      JP=J+1
      UV1=0.0
      UV2=0.0
      UV3=0.0
      UV4=0.0
      UV1P=0.0
      UV2P=0.0
      UV3P=0.0
      UV4P=0.0
      DO 240 K=1,N3M
      DO 240 I=1,N1M
      IP=IPA(I)
      V1=(U(I,J,K,1)+U(IP,J,K,1))*0.5
      V2=(U(I,J,K,2)+U(I,JP,K,2))*0.5
      V1P=V1-VMP(J,1)
      V2P=V2-VMP(J,2)
      IF(V1P .GT. 0.0 .AND. V2P .GT. 0.0) THEN
      UV1=UV1+V1P*V2P 
      UV1P=UV1P+1.0 
      ENDIF 
      IF(V1P .LT. 0.0 .AND. V2P .GT. 0.0) THEN
      UV2=UV2+V1P*V2P 
      UV2P=UV2P+1.0 
      ENDIF 
      IF(V1P .LT. 0.0 .AND. V2P .LT. 0.0) THEN
      UV3=UV3+V1P*V2P 
      UV3P=UV3P+1.0 
      ENDIF 
      IF(V1P .GT. 0.0 .AND. V2P .LT. 0.0) THEN
      UV4=UV4+V1P*V2P 
      UV4P=UV4P+1.0 
      ENDIF 
  240 CONTINUE
      Q(J,1)=-UV1*NXZ
      Q(J,2)=-UV2*NXZ
      Q(J,3)=-UV3*NXZ
      Q(J,4)=-UV4*NXZ
      QP(J,1)=UV1P*NXZ
      QP(J,2)=UV2P*NXZ
      QP(J,3)=UV3P*NXZ
      QP(J,4)=UV4P*NXZ
  230 CONTINUE
      
      
      
CCCCCCCCCCCCCCCCC BUDGET ANALYSIS CCCCCCCCCCCCCCCCCCC   
      
C     CALCULATE p*du/dx, p*dw/dz, (du/dx)^2, (dv/dx)^2, (dw/dz)^2
!c!$omp parallel do private(KP,IP,P_DUDX,P_DWDZ,DU_DX2,DV_DX2,DW_DZ2,
!c!$omp& U1,U2,V1,V2,W1,W2,V_IP,V_IC,V_IM)
!      DO 60 J=1,N2M
!      P_DUDX=0.0
!      P_DWDZ=0.0
!      DU_DX2=0.0
!      DV_DX2=0.0
!      DW_DZ2=0.0
!      DO 70 K=1,N3M
!      KP=KPA(K)
!      DO 70 I=1,N1M
!      IP=IPA(I)
!      IM=IMA(I)
!      U1=U(I ,J,K,1)
!      U2=U(IP,J,K,1)
!      V_IP=(U(IP,J,K,2)+U(IP,JP,K,2))*0.5
!      V_IC=(U(I ,J,K,2)+U(I ,JP,K,2))*0.5
!      V_IM=(U(IM,J,K,2)+U(IM,JP,K,2))*0.5
!      V2=(V_IP+V_IC)*0.5
!      V1=(V_IM+V_IC)*0.5
!      W1=U(I,J,K ,3)
!      W2=U(I,J,KP,3)
!      P_DUDX=P_DUDX+P(I,J,K)*(U2-U1)*DX1
!      P_DWDZ=P_DWDZ+P(I,J,K)*(W2-W1)*DX3
!      DU_DX2=DU_DX2+((U2-U1)*DX1)**2
!      DV_DX2=DV_DX2+((V2-V1)*DX1)**2
!      DW_DZ2=DW_DZ2+((W2-W1)*DX3)**2
!   70 CONTINUE
!      PDUDX(J)=P_DUDX*NXZ
!      PDWDZ(J)=P_DWDZ*NXZ
!      DUDX2(J)=DU_DX2*NXZ
!      DVDX2(J)=DV_DX2*NXZ
!      DWDZ2(J)=DW_DZ2*NXZ
!   60 CONTINUE

C     CALCULATE p*dv/dy,(du/dy)^2,(dv/dy)^2,(dw/dy)^2
!c!$omp parallel do private(JP,JM,FJUM,FJUP,P_DVDY,KP,IP,
!c!$omp& DU_DY2,DV_DY2,DW_DY2,U_JP,U_JC,U_JM,W_JP,W_JC,W_JM,
!c!$omp& U1,U2,V1,V2,W1,W2)
!      DO 80 J=1,N2M
!      JP=J+1
!      JM=J-1
!      FJUM=FJMU(J)
!      FJUP=FJPA(J)
!      P_DVDY=0.0
!      DU_DY2=0.0
!      DV_DY2=0.0
!      DW_DY2=0.0
!      DO 90 K=1,N3M
!      KP=KPA(K)
!      DO 90 I=1,N1M
!      IP=IPA(I)

!      U_JP=(U(I,JP,K,1)+U(IP,JP,K,1))*0.5
!      U_JC=(U(I,J ,K,1)+U(IP,J ,K,1))*0.5
!      U_JM=(U(I,JM,K,1)+U(IP,JM,K,1))*0.5
!      U1=FJUM*(1.0/H(J )*(DY(J)*U_JM/2.0+DY(JM)*U_JC/2.0))
!     >  +(1.-FJUM)*(U(I,0,K,1)+U(IP,0,K,1))*0.5
!      U2=FJUP*(1.0/H(JP)*(DY(J)*U_JP/2.0+DY(JP)*U_JC/2.0))
!     >  +(1.-FJUP)*(U(I,N2,K,1)+U(IP,N2,K,1))*0.5
!      V1=U(I,J ,K,2)
!      V2=U(I,JP,K,2)
!      W_JP=(U(I,JP,K,3)+U(I,JP,KP,3))*0.5
!      W_JC=(U(I,J ,K,3)+U(I,J ,KP,3))*0.5
!      W_JM=(U(I,JM,K,3)+U(I,JM,KP,3))*0.5
!      W1=FJUM*(1.0/H(J )*(DY(J)*W_JM/2.0+DY(JM)*W_JC/2.0))
!     >  +(1.-FJUM)*(U(I,0,K,3)+U(I,0,KP,3))*0.5
!      W2=FJUP*(1.0/H(JP)*(DY(J)*W_JP/2.0+DY(JP)*W_JC/2.0))
!     >  +(1.-FJUP)*(U(I,N2,K,3)+U(I,N2,KP,3))*0.5

!      P_DVDY=P_DVDY+P(I,J,K)*(V2-V1)/DY(J)
!      DU_DY2=DU_DY2+((U2-U1)/DY(J))**2
!      DV_DY2=DV_DY2+((V2-V1)/DY(J))**2
!      DW_DY2=DW_DY2+((W2-W1)/DY(J))**2
!   90 CONTINUE
!      PDVDY(J)=P_DVDY*NXZ
!      DUDY2(J)=DU_DY2*NXZ
!      DVDY2(J)=DV_DY2*NXZ
!      DWDY2(J)=DW_DY2*NXZ
!   80 CONTINUE

C     CALCULATE (dw/dx)^2,(du/dz)^2,(dv/dz)^2
!c!$omp parallel do private(JP,KP,KM,IP,IM,U_KP,U_KC,U_KM,
!c!$omp& V_KP,V_KC,V_KM,W_IP,W_IC,W_IM,DW_DX2,DU_DZ2,DV_DZ2,
!c!$omp& U1,U2,V1,V2,W1,W2)

!      DO 100 J=1,N2M
!      DW_DX2=0.0
!      DU_DZ2=0.0
!      DV_DZ2=0.0
!      JP=J+1
!      DO 110 K=1,N3M
!      KP=KPA(K)
!      KM=KMA(K)
!      DO 110 I=1,N1M
!      IP=IPA(I)
!      IM=IMA(I)
!      U_KP=(U(I,J,KP,1)+U(IP,J,KP,1))*0.5
!      U_KC=(U(I,J,K ,1)+U(IP,J,K ,1))*0.5
!      U_KM=(U(I,J,KM,1)+U(IP,J,KM,1))*0.5
!      V_KP=(U(I,J,KP,2)+U(I,JP,KP,2))*0.5
!      V_KC=(U(I,J,K ,2)+U(I,JP,K ,2))*0.5
!      V_KM=(U(I,J,KM,2)+U(I,JP,KM,2))*0.5
!      W_IP=(U(IP,J,K,3)+U(IP,J,KP,3))*0.5
!      W_IC=(U(I ,J,K,3)+U(I ,J,KP,3))*0.5
!      W_IM=(U(IM,J,K,3)+U(IM,J,KP,3))*0.5
!      U2=(U_KP+U_KC)*0.5
!      U1=(U_KM+U_KC)*0.5
!      V2=(V_KP+V_KC)*0.5
!      V1=(V_KM+V_KC)*0.5
!      W2=(W_IP+W_IC)*0.5
!      W1=(W_IM+W_IC)*0.5
!      DW_DX2=DW_DX2+((W2-W1)*DX1)**2
!      DU_DZ2=DU_DZ2+((U2-U1)*DX3)**2
!      DV_DZ2=DV_DZ2+((V2-V1)*DX3)**2
!  110 CONTINUE
!      DWDX2(J)=DW_DX2*NXZ
!      DUDZ2(J)=DU_DZ2*NXZ
!      DVDZ2(J)=DV_DZ2*NXZ
!  100 CONTINUE

C     CALCULATE Time and spatially averaged mean values
!      DO 120 J=1,N2M
!      IF(NAVG.EQ.1) THEN
!      DO 125 M=1,3
!  125 VMP(J,M)=VM(J,M)
!      ENDIF
!      IF(NAVG.NE.1) THEN
!      DO 127 M=1,3
!  127 VMP(J,M)=(VMO(J,M)+VM(J,M))/RNAVG
!      ENDIF
!  120 CONTINUE

      call dfftw_destroy_plan(FWD1)
      call dfftw_destroy_plan(FWD3)
      call dfftw_destroy_plan(FWD)

      RETURN
      END


C  ****************************** TAVER **********************
C     CALCULATE TIME AVERAGED MEAN QUANTITIES :

      SUBROUTINE TAVER(NAVG)
      INCLUDE 'PARAM.H'
      
C-----FOR MEAN VALUES      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/WMEAN/WM(M2,3),WRMS(M2,3)
      COMMON/PMEAN/PM(M2),PRMS(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)
      
      COMMON/SPEC/XSPCTR(M1,M2,3),ZSPCTR(M3,M2,3)
      COMMON/SPECO/XSPCTRO(M1,M2,3),ZSPCTRO(M3,M2,3)
      COMMON/SPECTAU/SPECTAUW(M1,M3,2),SPECTAUB(M1,M3,2)
      COMMON/SPECTAUO/SPECTAUWO(M1,M3,2),SPECTAUBO(M1,M3,2)
      COMMON/CORR/RXX(M1,M2,3),RZZ(M3,M2,3)
      COMMON/CORRO/RXXO(M1,M2,3),RZZO(M3,M2,3)
      COMMON/CORR2/RXY(M1,M2,3)
      COMMON/CORR2O/RXYO(M1,M2,3)
      COMMON/CORR3/RXYZ(M1,M2,M3,3)
      COMMON/CORR3O/RXYZO(M1,M2,M3,3)
      COMMON/QUAD/Q(M2,4),QP(M2,4)     
      COMMON/QUADO/QO(M2,4),QPO(M2,4)   
      
c-----for turbulent kinetic budgets 
      COMMON/VMEAN2/DUDY(M2),DVDY(M2),DWDY(M2)
      COMMON/PVMEAN/PV(M2)
      COMMON/VHIGH/VSKEW(M2,3),VFLAT(M2,3),U2V(M2),VW2(M2)
      COMMON/PSTR/PDUDX(M2),PDVDY(M2),PDWDZ(M2)
      COMMON/DISSU/DUDX2(M2),DUDY2(M2),DUDZ2(M2)
      COMMON/DISSV/DVDX2(M2),DVDY2(M2),DVDZ2(M2)
      COMMON/DISSW/DWDX2(M2),DWDY2(M2),DWDZ2(M2)      
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEAN/STRM(M2,6),SGSNUT(M2),SGSDS(M2),SMAGC(M2)
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WMEANO/WMO(M2,3),WRMSO(M2,3)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      
c-----for turbulent kinetic budgets 
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/PVMEANO/PVO(M2)
      COMMON/VHIGHO/VSKEWO(M2,3),VFLATO(M2,3),U2VO(M2),VW2O(M2)
      COMMON/PSTRO/PDUDXO(M2),PDVDYO(M2),PDWDZO(M2)
      COMMON/DISSUO/DUDX2O(M2),DUDY2O(M2),DUDZ2O(M2)
      COMMON/DISSVO/DVDX2O(M2),DVDY2O(M2),DVDZ2O(M2)
      COMMON/DISSWO/DWDX2O(M2),DWDY2O(M2),DWDZ2O(M2)      
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLTAO/TAOW
      COMMON/WALLTAOAVE/TAOWO
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/TAUXYZO/TAUWO(0:M1,0:M3,2),TAUBO(0:M1,0:M3,2)
      
      IF(NAVG.EQ.1) THEN
          
!$omp parallel do 
      DO 10 J=1,N2M
      DO 10 L=1,3      
      VMO(J,L)=VM(J,L)            ! velocity
      WMO(J,L)=WM(J,L)            ! vorticity
      WRMSO(J,L)=WRMS(J,L)        ! enstrophy
      VSKEWO(J,L)=VSKEW(J,L)      ! skewness
      VFLATO(J,L)=VFLAT(J,L)      ! flatness
   10 CONTINUE
      
!$omp parallel do     
      DO J=1,N2M
      DO I=1,N1M
      DO L=1,3
      XSPCTRO(I,J,L)=XSPCTR(I,J,L)
      ENDDO
      ENDDO
      ENDDO
      
!$omp parallel do       
      DO J=1,N2M
      DO K=1,N3M
      DO L=1,3
      ZSPCTRO(K,J,L)=ZSPCTR(K,J,L)
      ENDDO
      ENDDO
      ENDDO 
      
!$omp parallel do       
      DO I=1,N1M
      DO K=1,N3M
      DO L=1,2
      SPECTAUWO(I,K,L)=SPECTAUW(I,K,L)
      ENDDO
      ENDDO
      ENDDO        

!$omp parallel do 
      DO 20 L=1,3
      WSMO(L)=WSM(L)          ! wall shear stress
   20 CONTINUE
   
      TAOWO=TAOW

!$omp parallel do 
      DO 30 J=1,N2M
      DO 30 L=1,6
      VRMSO(J,L)=VRMS(J,L)
      STRMO(J,L)=STRM(J,L)
   30 CONTINUE
      
!$omp parallel do 
      DO  K=1,N3M
      DO  I=1,N1M
      DO  L=1,2
      TAUWO(I,K,L)=TAUW(I,K,L)
      TAUBO(I,K,L)=TAUB(I,K,L)
      ENDDO
      ENDDO
      ENDDO
      
!$omp parallel do 
      DO 40 J=1,N2M
          PMO(J)=PM(J)
          PRMSO(J)=PRMS(J)
          PVO(J)=PV(J)
          DUDYO(J)=DUDY(J)
          DVDYO(J)=DVDY(J)
          DWDYO(J)=DWDY(J)
          U2VO(J)=U2V(J)
          VW2O(J)=VW2(J)
          PDUDXO(J)=PDUDX(J)
          PDVDYO(J)=PDVDY(J)
          PDWDZO(J)=PDWDZ(J)
          DUDX2O(J)=DUDX2(J)
          DUDY2O(J)=DUDY2(J)
          DUDZ2O(J)=DUDZ2(J)
          DVDX2O(J)=DVDX2(J)
          DVDY2O(J)=DVDY2(J)
          DVDZ2O(J)=DVDZ2(J)
          DWDX2O(J)=DWDX2(J)
          DWDY2O(J)=DWDY2(J)
          DWDZ2O(J)=DWDZ2(J)
          SGSNUTO(J)=SGSNUT(J)
          SGSDSO(J)=SGSDS(J)
          SMAGCO(J)=SMAGC(J)
   40 CONTINUE

!$omp parallel do 
      DO 50 J=1,N2M
      DO 50 L=1,4      
          QO(J,L)=Q(J,L)
          QPO(J,L)=QP(J,L)
   50 CONTINUE

!$omp parallel do 
      DO 60 I=1,N1M
      DO 60 J=1,N2M
      DO 60 L=1,3
          RXXO(I,J,L)=RXX(I,J,L)
   60 CONTINUE

!$omp parallel do 
      DO 70 K=1,N3M
      DO 70 J=1,N2M
      DO 70 L=1,3      
          RZZO(K,J,L)=RZZ(K,J,L)
   70 CONTINUE
       
      
!$omp parallel do 
      DO 71 I=1,N1M
      DO 71 J=1,N2M
      DO 71 L=1,3      
          RXYO(I,J,L)=RXY(I,J,L)
   71 CONTINUE
      
      
!$omp parallel do 
      DO 72 K=1,N3M
      DO 72 I=1,N1M
      DO 72 J=1,N2M
      DO 72 L=1,3      
          RXYZO(I,J,K,L)=RXYZ(I,J,K,L)
   72 CONTINUE      
      
      ENDIF 

      IF(NAVG.NE.1) THEN
!$omp parallel do 
      DO 80 J=1,N2M
      DO 80 L=1,3      
          VMO(J,L)=(VMO(J,L)*REAL(NAVG-1)+VM(J,L))/REAL(NAVG)
          WMO(J,L)=(WMO(J,L)*REAL(NAVG-1)+WM(J,L))/REAL(NAVG)
          WRMSO(J,L)=(WRMSO(J,L)*REAL(NAVG-1)+WRMS(J,L))/REAL(NAVG)
          VSKEWO(J,L)=(VSKEWO(J,L)*REAL(NAVG-1)+VSKEW(J,L))/REAL(NAVG)
          VFLATO(J,L)=(VFLATO(J,L)*REAL(NAVG-1)+VFLAT(J,L))/REAL(NAVG)
   80 CONTINUE
      
!$omp parallel do      
      DO J=1,N2M
      DO I=1,N1M
      DO L=1,3
      XSPCTRO(I,J,L)=(XSPCTRO(I,J,L)*REAL(NAVG-1)+XSPCTR(I,J,L))
     >                    /REAL(NAVG)
      ENDDO
      ENDDO
      ENDDO
      
!$omp parallel do      
      DO J=1,N2M
      DO K=1,N3M
      DO L=1,3
      ZSPCTRO(K,J,L)=(ZSPCTRO(K,J,L)*REAL(NAVG-1)+ZSPCTR(K,J,L))
     >                    /REAL(NAVG)
      ENDDO
      ENDDO
      ENDDO 
!$omp parallel do
      DO I=1,N1M
      DO K=1,N3M
      DO L=1,2
      SPECTAUWO(I,K,L)=(SPECTAUWO(I,K,L)*REAL(NAVG-1)+SPECTAUW(I,K,L))
     >                    /REAL(NAVG)
      ENDDO
      ENDDO
      ENDDO
      
!$omp parallel do      
      DO I=1,N1M
      DO K=1,N3M
      DO L=1,2
      TAUWO(I,K,L)=(TAUWO(I,K,L)*REAL(NAVG-1)+TAUW(I,K,L))
     >                    /REAL(NAVG)
      ENDDO
      ENDDO
      ENDDO 

!$omp parallel do 
      DO 90 L=1,3
          WSMO(L)=(WSMO(L)*REAL(NAVG-1)+WSM(L))/REAL(NAVG)
   90 CONTINUE

      TAOWO=(TAOWO*REAL(NAVG-1)+TAOW)/REAL(NAVG)


!$omp parallel do      
      DO 100 J=1,N2M
      DO 100 L=1,6
          VRMSO(J,L)=(VRMSO(J,L)*REAL(NAVG-1)+VRMS(J,L))/REAL(NAVG)
          STRMO(J,L)=(STRMO(J,L)*REAL(NAVG-1)+STRM(J,L))/REAL(NAVG)
  100 CONTINUE

!$omp parallel do 
      DO 110 J=1,N2M
          PMO(J)=(PMO(J)*REAL(NAVG-1)+PM(J))/REAL(NAVG)
          PRMSO(J)=(PRMSO(J)*REAL(NAVG-1)+PRMS(J))/REAL(NAVG)
          PVO(J)=(PVO(J)*REAL(NAVG-1)+PV(J))/REAL(NAVG)
          DUDYO(J)=(DUDYO(J)*REAL(NAVG-1)+DUDY(J))/REAL(NAVG)
          DVDYO(J)=(DVDYO(J)*REAL(NAVG-1)+DVDY(J))/REAL(NAVG)
          DWDYO(J)=(DWDYO(J)*REAL(NAVG-1)+DWDY(J))/REAL(NAVG)
          U2VO(J)=(U2VO(J)*REAL(NAVG-1)+U2V(J))/REAL(NAVG)
          VW2O(J)=(VW2O(J)*REAL(NAVG-1)+VW2(J))/REAL(NAVG)
          PDUDXO(J)=(PDUDXO(J)*REAL(NAVG-1)+PDUDX(J))/REAL(NAVG)
          PDVDYO(J)=(PDVDYO(J)*REAL(NAVG-1)+PDVDY(J))/REAL(NAVG)
          PDWDZO(J)=(PDWDZO(J)*REAL(NAVG-1)+PDWDZ(J))/REAL(NAVG)
          DUDX2O(J)=(DUDX2O(J)*REAL(NAVG-1)+DUDX2(J))/REAL(NAVG)
          DUDY2O(J)=(DUDY2O(J)*REAL(NAVG-1)+DUDY2(J))/REAL(NAVG)
          DUDZ2O(J)=(DUDZ2O(J)*REAL(NAVG-1)+DUDZ2(J))/REAL(NAVG)
          DVDX2O(J)=(DVDX2O(J)*REAL(NAVG-1)+DVDX2(J))/REAL(NAVG)
          DVDY2O(J)=(DVDY2O(J)*REAL(NAVG-1)+DVDY2(J))/REAL(NAVG)
          DVDZ2O(J)=(DVDZ2O(J)*REAL(NAVG-1)+DVDZ2(J))/REAL(NAVG)
          DWDX2O(J)=(DWDX2O(J)*REAL(NAVG-1)+DWDX2(J))/REAL(NAVG)
          DWDY2O(J)=(DWDY2O(J)*REAL(NAVG-1)+DWDY2(J))/REAL(NAVG)
          DWDZ2O(J)=(DWDZ2O(J)*REAL(NAVG-1)+DWDZ2(J))/REAL(NAVG)
          SGSNUTO(J)=(SGSNUTO(J)*REAL(NAVG-1)+SGSNUT(J))/REAL(NAVG)
          SGSDSO(J)=(SGSDSO(J)*REAL(NAVG-1)+SGSDS(J))/REAL(NAVG)
          SMAGCO(J)=(SMAGCO(J)*REAL(NAVG-1)+SMAGC(J))/REAL(NAVG)
  110 CONTINUE

!$omp parallel do      
      DO 120 J=1,N2M
      DO 120 L=1,4
          QO(J,L)=(QO(J,L)*REAL(NAVG-1)+Q(J,L))/REAL(NAVG)
          QPO(J,L)=(QPO(J,L)*REAL(NAVG-1)+QP(J,L))/REAL(NAVG)
  120 CONTINUE

!$omp parallel do 
      DO 130 I=1,N1M
      DO 130 J=1,N2M
      DO 130 L=1,3      
          RXXO(I,J,L)=(RXXO(I,J,L)*REAL(NAVG-1)+RXX(I,J,L))/REAL(NAVG)
  130 CONTINUE

!$omp parallel do 
      DO 140 K=1,N3M
      DO 140 J=1,N2M
      DO 140 L=1,3      
          RZZO(K,J,L)=(RZZO(K,J,L)*REAL(NAVG-1)+RZZ(K,J,L))/REAL(NAVG)
  140 CONTINUE
       
      
!$omp parallel do 
      DO 131 I=1,N1M
      DO 131 J=1,N2M
      DO 131 L=1,3      
          RXYO(I,J,L)=(RXYO(I,J,L)*REAL(NAVG-1)+RXY(I,J,L))/REAL(NAVG)
  131 CONTINUE
      
!$omp parallel do 
      DO 132 K=1,N3M
      DO 132 I=1,N1M
      DO 132 J=1,N2M
      DO 132 L=1,3      
          RXYZO(I,J,K,L)=(RXYZO(I,J,K,L)*REAL(NAVG-1)
     >     +RXYZ(I,J,K,L))/REAL(NAVG)
  132 CONTINUE     
      
      ENDIF

      RETURN
      END

C  ****************************** WALLSS **********************
C     CALCULATE WALL SHEAR STRESS and Pressure at Top and Bottom WALL

      SUBROUTINE WALLSS(U,P,TIME)
      INCLUDE 'PARAM.H'

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/PARA/RE
      COMMON/WSMEM/WSM(3),WSMO(3)
      COMMON/SIZE/ALX,ALY,ALZ,VOL

      REAL U(0:M1,0:M2,0:M3,3),P(M1,M2,M3)
      REAL WSB(M1,M3,3),WST(M1,M3,3)
      REAL NXZ

      NXZ=1.0/DBLE(N1M*N3M)

C     Bottom and Top WALL at I+1/2,J,K+1/2 Points
      WSM(1)=0.0
      WSM(2)=0.0
      WSM(3)=0.0

!$omp parallel do private(KP,KM,IP,IM,UHC,UHW,VHC,VHW,WHC,WHW,
!$omp& WSB,WST)
!$omp& reduction(+:WSM)
      DO 10 K=1,N3M
      KP=KPA(K)
      KM=KMA(K)
      DO 10 I=1,N1M
      IP=IPA(I)
      IM=IMA(I)

      UHC=0.5*(U(I,1,K,1)+U(IP,1,K,1))     ! Bottom Wall
      UHW=0.5*(U(I,0,K,1)+U(IP,0,K,1))
      VHC=U(I,2,K,2)
      VHW=U(I,1,K,2)
      WHC=0.5*(U(I,1,K,3)+U(I,1,KP,3))
      WHW=0.5*(U(I,0,K,3)+U(I,0,KP,3))
      WSB(I,K,1)=(UHC-UHW)/(DY(1)*0.5)
      WSB(I,K,2)=(VHC-VHW)/DY(1)
      WSB(I,K,3)=(WHC-WHW)/(DY(1)*0.5)

      UHC=0.5*(U(I,N2M,K,1)+U(IP,N2M,K,1)) ! Top Wall
      UHW=0.5*(U(I,N2 ,K,1)+U(IP,N2 ,K,1))
      VHC=U(I,N2M,K,2)
      VHW=U(I,N2 ,K,2)
      WHC=0.5*(U(I,N2M,K,3)+U(I,N2M,KP,3))
      WHW=0.5*(U(I,N2 ,K,3)+U(I,N2 ,KP,3))
      WST(I,K,1)=(UHC-UHW)/(DY(N2M)*0.5)
      WST(I,K,2)=(VHC-VHW)/DY(N2M)
      WST(I,K,3)=(WHC-WHW)/(DY(N2M)*0.5)

      WSM(1)=WSM(1)+WSB(I,K,1)+WST(I,K,1)*(ALY-1.0) !)+WST(I,K,1)) half-ch 0828
      WSM(2)=WSM(2)+WSB(I,K,2)+WST(I,K,1)*(ALY-1.0) !)+WST(I,K,2))
      WSM(3)=WSM(3)+WSB(I,K,3)+WST(I,K,1)*(ALY-1.0) !)+WST(I,K,3))
   10 CONTINUE
      WSM(1)=WSM(1)*NXZ/RE/ALY
      WSM(2)=WSM(2)*NXZ/RE/ALY
      WSM(3)=WSM(3)*NXZ/RE/ALY
     
!         WRITE(*,100) TIME,SQRT(ABS(WSM(1)))*RE,SQRT(ABS(WSM(1))),RE
  100 FORMAT(4(E12.5,2X))
  
!     CALL WALLSS_WRITE(U,P,PRESG,TIME)
      RETURN
      END

      
      
      SUBROUTINE WALLSS_WRITE(U,P,PRESG,TIME)
      INCLUDE 'PARAM.H'

      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)
      COMMON/WSMEM/WSM(3),WSMO(3)
        COMMON/PARA/RE
        REAL PRESG
      WRITE(31,100) TIME,WSM(1),-PRESG 
C      WRITE(31,100) TIME,WSM(1),WSM(2),WSM(3)
      WRITE(315,100) TIME,SQRT(ABS(WSM(1)))*RE,SQRT(ABS(WSM(1)))
  100 FORMAT(4(E12.5,2X))

      RETURN
      END


C  **********************  WRITEAVG **********************
C     WRITE AVERAGED VALUES FOR TURBULENT STATISTICS

      SUBROUTINE WRITEAVG(U,P,TIME,NAVG,IPLUS,NTIME)
      INCLUDE 'PARAM.H'
      
      COMMON/DIM/N1,N2,N3,N1M,N2M,N3M
      COMMON/MESH1/DX1,DX1Q,DX3,DX3Q
      COMMON/MESH3/DY(0:M2),H(M2),HM(M2),HC(M2),HP(M2)
      COMMON/PARA/RE      
      
      COMMON/VMEANO/VMO(M2,3),VRMSO(M2,6)
      COMMON/WMEANO/WMO(M2,3),WRMSO(M2,3)
      COMMON/PMEANO/PMO(M2),PRMSO(M2)
      COMMON/VMEAN/VM(M2,3),VRMS(M2,6)
      COMMON/WSMEM/WSM(3),WSMO(3)
      
      COMMON/SPEC/XSPCTR(M1,M2,3),ZSPCTR(M3,M2,3)
      COMMON/SPECO/XSPCTRO(M1,M2,3),ZSPCTRO(M3,M2,3)
      COMMON/CORR/RXX(M1,M2,3),RZZ(M3,M2,3)
      COMMON/CORRO/RXXO(M1,M2,3),RZZO(M3,M2,3)
      COMMON/CORR2/RXY(M1,M2,3)
      COMMON/CORR2O/RXYO(M1,M2,3)
      COMMON/CORR3/RXYZ(M1,M2,M3,3)
      COMMON/CORR3O/RXYZO(M1,M2,M3,3)   
      COMMON/QUAD/Q(M2,4),QP(M2,4)     
      COMMON/QUADO/QO(M2,4),QPO(M2,4)   
      COMMON/TAUXYZO/TAUWO(0:M1,0:M3,2),TAUBO(0:M1,0:M3,2)
      COMMON/SPECTAU/SPECTAUW(M1,M3,2),SPECTAUB(M1,M3,2)
      COMMON/SPECTAUO/SPECTAUWO(M1,M3,2),SPECTAUBO(M1,M3,2)
      
c-----for turbulent kinetic budgets 
      COMMON/VMEAN2O/DUDYO(M2),DVDYO(M2),DWDYO(M2)
      COMMON/PVMEANO/PVO(M2)
      COMMON/VHIGHO/VSKEWO(M2,3),VFLATO(M2,3),U2VO(M2),VW2O(M2)
      COMMON/PSTRO/PDUDXO(M2),PDVDYO(M2),PDWDZO(M2)
      COMMON/DISSUO/DUDX2O(M2),DUDY2O(M2),DUDZ2O(M2)
      COMMON/DISSVO/DVDX2O(M2),DVDY2O(M2),DVDZ2O(M2)
      COMMON/DISSWO/DWDX2O(M2),DWDY2O(M2),DWDZ2O(M2)
      COMMON/PROFIL/UM_STEP(0:M2),VM_STEP(M2),WM_STEP(0:M2)
      
c-----for SGS stress and dissipation 
      COMMON/SGSMEANO/STRMO(M2,6),SGSNUTO(M2),SGSDSO(M2),SMAGCO(M2)
      COMMON/WALLMODEL/WMODEL,IFWALL,TAUW1,TAUB1,IS
      COMMON/TAUXYZ/TAUW(0:M1,0:M3,2),TAUB(0:M1,0:M3,2)
      COMMON/FINOUT/INCODE,IDTOPT,NWRITE,NREAD,IAVG,NPRN,INSF,NINS,FILES

      COMMON/WALLTAO/TAOW
      COMMON/WALLTAOAVE/TAOWO

      REAL U(3,0:M1,0:M2,0:M3)
      REAL P(M1,M2,M3)
      REAL UTAU,TAUW0(N1M,N3M)
      REAL A,B,DXN,T1,T2,TT(M1)
      REAL T(M1),PDF(M1)
      
       CHARACTER*80  FILEW3,FILEW1,FILEW2,FILEW4,FILEW5,
     >                FILEW6,FILEW7,FILEW8,FILEW9,FILEW10,
     >                FILEW11,FILEW12,FILEW13,FILEW14,FILEW15,
     >                FILEW16,FILEW17,FILEW18
      
C      IF(1.EQ.1)THEN

!      UTAU=SQRT((UM_STEP(1)+UM_STEP(N2M))/RE/dy(1))
      IF(WMODEL.EQ.0)THEN
      UTAU=SQRT(WSMO(1))
      ELSE
      UTAU=SQRT(TAOWO)
      ENDIF      
!    UTAU=SQRT(WSM(1))
!      print*,utau
!      pause 123
!      WRITE(*,*)TAOWO,UTAU
      
      FILEW1='AVE/AVEUVW.'
      N=INDEX(FILEW1,'.')
      WRITE(UNIT=FILEW1(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(41,FILE=FILEW1,STATUS='UNKNOWN')

      FILEW2='AVE/AVERMS.'
      N=INDEX(FILEW2,'.')
      WRITE(UNIT=FILEW2(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(42,FILE=FILEW2,STATUS='UNKNOWN')

!      FILEW3='VOR_AVE_RMS.'
!      N=INDEX(FILEW3,'.')
!      WRITE(UNIT=FILEW3(N+1:),FMT='(BN,I6.6)') NTIME
!      OPEN(43,FILE=FILEW3,STATUS='UNKNOWN')


!      FILEW4='P_RMS.'
!      N=INDEX(FILEW4,'.')
!      WRITE(UNIT=FILEW4(N+1:),FMT='(BN,I6.6)') NTIME
!      OPEN(44,FILE=FILEW4,STATUS='UNKNOWN')

      FILEW5='AVE/SGS_AVE.'
      N=INDEX(FILEW5,'.')
      WRITE(UNIT=FILEW5(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(45,FILE=FILEW5,STATUS='UNKNOWN')
      
      FILEW6='AVE/MEANSU.'
      N=INDEX(FILEW6,'.')
      WRITE(UNIT=FILEW6(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(46,FILE=FILEW6,STATUS='UNKNOWN')
      
      FILEW7='AVE/XSPCTR.'
      N=INDEX(FILEW7,'.')
      WRITE(UNIT=FILEW7(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(47,FILE=FILEW7,STATUS='UNKNOWN')
      
      FILEW10='AVE/ZSPCTR.'
      N=INDEX(FILEW10,'.')
      WRITE(UNIT=FILEW10(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(50,FILE=FILEW10,STATUS='UNKNOWN')
      
      FILEW11='AVE/XCORR.'
      N=INDEX(FILEW11,'.')
      WRITE(UNIT=FILEW11(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(51,FILE=FILEW11,STATUS='UNKNOWN')
      
      FILEW12='AVE/ZCORR.'
      N=INDEX(FILEW12,'.')
      WRITE(UNIT=FILEW12(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(52,FILE=FILEW12,STATUS='UNKNOWN')
      
      FILEW13='AVE/XYCORR.'
      N=INDEX(FILEW13,'.')
      WRITE(UNIT=FILEW13(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(53,FILE=FILEW13,STATUS='UNKNOWN')
      
      FILEW14='AVE/XYZCORR.'
      N=INDEX(FILEW14,'.')
      WRITE(UNIT=FILEW14(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(54,FILE=FILEW14,STATUS='UNKNOWN')
      WRITE(54,*)'VARIABLES="X","Y","Z","RUU","RVV","RWW"'
      WRITE(54,*)'ZONE T="',1,'" I=',N1M ,' J=',N2M,' K=',N3M
      
      FILEW15='AVE/QUAD.'
      N=INDEX(FILEW15,'.')
      WRITE(UNIT=FILEW15(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(55,FILE=FILEW15,STATUS='UNKNOWN')
      
      FILEW16='AVE/SKEWNESS.'
      N=INDEX(FILEW16,'.')
      WRITE(UNIT=FILEW16(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(56,FILE=FILEW16,STATUS='UNKNOWN')
      
      FILEW17='AVE/STRAINPDF.'
      N=INDEX(FILEW17,'.')
      WRITE(UNIT=FILEW17(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(57,FILE=FILEW17,STATUS='UNKNOWN')
      
      FILEW18='AVE/SPCTRTAUW.'
      N=INDEX(FILEW18,'.')
      WRITE(UNIT=FILEW18(N+1:),FMT='(BN,I6.6)') NTIME
      OPEN(58,FILE=FILEW18,STATUS='UNKNOWN')
      
       
      X2=0.0

      DO 50 J=1,N2M
          
      X2=X2+H(J)
      
      WRITE(41,100) REAL(X2*UTAU*RE),REAL(VMO(J,1)/UTAU),
     >    REAL(VMO(J,2)/UTAU),REAL(VMO(J,3)/UTAU)
      
      WRITE(46,100) REAL(X2),REAL(VMO(J,1)),
     >    REAL(VMO(J,2)),REAL(VMO(J,3))
      
      IF(FILES.NE.0.0)THEN
      WRITE(42,100) REAL(X2*UTAU*RE)
     >            ,REAL(SQRT(VRMSO(J,1)-VMO(J,1)**2)/UTAU)
     >             ,REAL(SQRT(VRMSO(J,2)-VMO(J,2)**2)/UTAU)
     >            ,REAL(SQRT(VRMSO(J,3)-VMO(J,3)**2)/UTAU)
     >             ,-REAL((VRMSO(J,4)-VMO(J,1)*VMO(J,2))/UTAU**2)
      ELSE
c      WRITE(42,100) REAL(X2*UTAU*RE)
c     >                        ,REAL(SQRT(VRMSO(J,1))/UTAU)
c     >                        ,REAL(SQRT(VRMSO(J,2))/UTAU)
c     >                        ,SQRT(VRMSO(J,3))/UTAU
c     >                 ,-(VRMSO(J,4))/UTAU**2          
      WRITE(42,100) REAL(X2*UTAU*RE)
     >                        ,REAL(SQRT(VRMSO(J,1)-VMO(J,1)**2)/UTAU)
     >                        ,REAL(SQRT(VRMSO(J,2)-VMO(J,2)**2)/UTAU)
     >                        ,SQRT(VRMSO(J,3)-VMO(J,3)**2)/UTAU
     >                 ,-(VRMSO(J,4)-VMO(J,1)*VMO(J,2))/UTAU**2
      ENDIF
!      WRITE(43,*) X2*UTAU*RE,(SQRT(WRMSO(J,1)))/UTAU**2/RE,
!     >                        (SQRT(WRMSO(J,2)))/UTAU**2/RE,
!     >                        (SQRT(WRMSO(J,3)))/UTAU**2/RE
!      WRITE(44,*) X2*UTAU*RE,SQRT(PRMSO(J)) !-PMO(J))

      UV=-REAL((VRMSO(J,4)-VMO(J,1)*VMO(J,2)))
       WRITE(55,100) REAL(X2),REAL(QO(J,1))/UV,
     >    REAL(QO(J,2))/UV,REAL(QO(J,3))/UV,REAL(QO(J,4))/UV,
     >    REAL(QPO(J,1)),REAL(QPO(J,2)),REAL(QPO(J,3)),
     >    REAL(QPO(J,4)) 
      
      DO I=1,N1M
      WRITE(47,100) REAL(X2*UTAU*RE),REAL(X2),XSPCTRO(I,J,1)/UTAU**2,
     >                 XSPCTRO(I,J,2)/UTAU**2,XSPCTRO(I,J,3)/UTAU**2
      ENDDO
      
      DO K=1,N3M
      WRITE(50,100) REAL(X2*UTAU*RE),REAL(X2),ZSPCTRO(K,J,1)/UTAU**2,
     >                 ZSPCTRO(K,J,2)/UTAU**2,ZSPCTRO(K,J,3)/UTAU**2
      ENDDO
      
      DO I=1,N1M
      X1=REAL(I-N1M/2-1)/DX1
      WRITE(51,100) REAL(X1),REAL(X2),RXXO(I,J,1),
     >                 RXXO(I,J,2),RXXO(I,J,3)
      ENDDO
      
      DO K=1,N3M
      Z1=REAL(K-N3M/2-1)/DX3
      WRITE(52,100) REAL(Z1),REAL(X2),RZZO(K,J,1),
     >                 RZZO(K,J,2),RZZO(K,J,3)
      ENDDO
      
      DO I=1,N1M
      X1=REAL(I-N1M/2-1)/DX1
      WRITE(53,100) REAL(X1),REAL(X2),RXYO(I,J,1),
     >                 RXYO(I,J,2),RXYO(I,J,3)
      ENDDO  
      
      WRITE(55,100) REAL(X2*UTAU*RE)
     >                        ,QO(J,1)/UTAU**2,QO(J,2)/UTAU**2,
     >                         QO(J,3)/UTAU**2,QO(J,4)/UTAU**2
      
      URMS1=REAL(SQRT(VRMSO(J,1)-VMO(J,1)**2)/UTAU)
      VRMS1=REAL(SQRT(VRMSO(J,2)-VMO(J,2)**2)/UTAU)
      WRMS1=REAL(SQRT(VRMSO(J,3)-VMO(J,3)**2)/UTAU)
      
      WRITE(56,100) REAL(X2*UTAU*RE)
     >            ,VSKEWO(J,1)/URMS1**3
     >            ,VSKEWO(J,2)/VRMS1**3
     >            ,VSKEWO(J,3)/WRMS1**3
      
C      ENDIF
      
   50 CONTINUE
      
      X2=0.0
      DO K=1,N3M
      X3=REAL(K-N3M/2-1)/DX3
      DO J=1,N2M
      X2=X2+H(J)
      DO I=1,N1M
      X1=REAL(I-N1M/2-1)/DX1      
      WRITE(54,100) REAL(X1),REAL(X2),REAL(X3),RXYZO(I,J,K,1),
     >                 RXYZO(I,J,K,2),RXYZO(I,J,K,3)
      ENDDO
      ENDDO
      ENDDO
      
         
      X2=0.0
      DO J=1,N2
      X2=X2+H(J)
      WRITE(45,*) REAL(X2),REAL(X2*UTAU*RE),REAL(SGSNUTO(J)/UTAU),
     >                REAL(SMAGCO(J))
      ENDDO
      
      DO I=1,N1M
      DO K=1,N3M
      WRITE(57,*) SPECTAUWO(I,K,1)/UTAU**2
      ENDDO
      ENDDO


      DO K=1,N3M
      DO I=1,N1M
      TAUW0(I,K) = (TAUWO(I,K,1) - UTAU**2)/(UTAU**2)
      ENDDO
      ENDDO
      
      DO M=1,N1M
      T(M)=0.
      ENDDO
      
      A=MINVAL(TAUW0)
      B=MAXVAL(TAUW0)
      DXN=(B-A)/N1M
      DO K=1,N3M
      DO I=1,N1M
      DO M=1,N1M
          T1=A+DXN*(M-1)
          T2=A+DXN*M
          IF(TAUW0(I,K).GE.T1.AND.TAUW0(I,K).LE.T2) THEN
          T(M) = T(M) + 1.
          ENDIF
      ENDDO
      ENDDO
      ENDDO
     
      DO M=1,N1M
          TT(M)=A+DXN*(2*M-1)/2.
          PDF(M)=T(M)/DBLE(N1M*N3M)/dxn                
          WRITE(58,*) TT(M),PDF(M)         
      ENDDO
      
  100 FORMAT(8(E12.5,2X))  
      
      RETURN
      END
      
      
!-----------------------linearly fitting----------------------------------
!     ndata is the number of experimental data
!     x is the array for independent variable
!     y is the array for dependent variable
!     y=a+bx
!      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
c      SUBROUTINE fit(x,y,ndata,a,b)
c      INTEGER mwt,ndata
c      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
c      !USES gammq
c      INTEGER i
c      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
c      sx=0.
c      sy=0.
c      st2=0.
c      b=0.
c      if(mwt/=0) then 
c      ss=0.
c      do i=1,ndata
c       wt=1./(sig(i)**2)
c       ss=ss+wt
c       sx=sx+x(i)*wt
c       sy=sy+y(i)*wt
c      end do
c      else
c      do i=1,ndata
c       sx=sx+x(i)
c       sy=sy+y(i)
c      end do
c      ss=float(ndata)
c      endif
c      sxoss=sx/ss
c      if(mwt/=0) then
c      do i=1,ndata
c       t=(x(i)-sxoss)/sig(i)
c       st2=st2+t*t
c       b=b+t*y(i)/sig(i)
c      end do
c      else
c      do i=1,ndata
c      t=x(i)-sxoss
c      st2=st2+t*t
c      b=b+t*y(i)
c      end do
c      endif
c      b=b/st2
c      a=(sy-sx*b)/ss
c      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
c      sigb=sqrt(1./st2)
c      chi2=0.
c      if(mwt==0) then
c      do i=1,ndata
c      chi2=chi2+(y(i)-a-b*x(i))**2
c      end do
c      q=1.
c      sigdat=sqrt(chi2/(ndata-2))
c      siga=siga*sigdat
c      sigb=sigb*sigdat
c      else
c      do i=1,ndata
c      chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
c      end do
c      q=gammq(0.5*(ndata-2),0.5*chi2)
c      endif
c      END SUBROUTINE fit
C--------------------DDHRT-------------------------------------
C THIS SUBROUTINE IS USED TO SOLVE THE POLYNOMIAN EQUATION
C OF THE SCALE-DEPENDENT COEFFICIENT IN SGSVIS MODEL
C A AND B ARE THE LOWER AND UPPER LIMITS OF THE RANGE OF THE ROOT
C X ARRAY OF ROOT
C N ESTIMATION OF THE NUMBER OF THE ROOT

	SUBROUTINE DQRRT(A,N,BETA,EPS)
	DIMENSION B(N,N),XR(N),XI(N),A(N)
	DOUBLE PRECISION B,XR,XI,A
	DO 10 J=1,N
10	B(1,J)=-A(J)
	DO 20 I=2,N
	DO 20 J=1,N
20	B(I,J)=0.0
	DO 30 I=1,N-1
30	B(I+1,I)=1.0
	CALL CHHQR(B,N,XR,XI,EPS,L)
	
	
	BETA=0.
	DO J=1,N
!	print*,"XR(J)=",XR(J)
	IF(XR(J).GT.BETA.AND.XR(J).LT.2.0) BETA=XR(J)
	END DO
!	pause
	
	RETURN
	END
C--------------------CHHQR-------------------------------------
	SUBROUTINE CHHQR(A,N,U,V,EPS,JT)
	DIMENSION A(N,N),U(N),V(N)
	DOUBLE PRECISION A,U,V,W,B,C,X,Y,XY,P,Q,R,E,F,Z,G
	JT=1
	M=N
	IT=0
10	IF (M.EQ.0) RETURN
	L=M+1
40	L=L-1
	IF (ABS(A(L,L-1)).GT.EPS*(ABS(A(L-1,L-1))+ABS(A(L,L))).AND.
     *      L.GT.1) GOTO 40
	IF (L.EQ.M) THEN
	  U(M)=A(M,M)
	  V(M)=0.0
	  M=M-1
	  IT=0
	  GOTO 10
	END IF
	IF (L.EQ.M-1) THEN
	  B=-(A(M,M)+A(M-1,M-1))
	  C=A(M,M)*A(M-1,M-1)-A(M,M-1)*A(M-1,M)
	  W=B*B-4*C
	  Y=SQRT(ABS(W))
	  IF (W.GT.0.0) THEN
	    XY=1.0
	    IF (B.LT.0.0) XY=-1.0
	    U(M)=(-B-XY*Y)/2.0
	    U(M-1)=C/U(M)
	    V(M)=0.0
	    V(M-1)=0.0
	  ELSE
	    U(M)=-B/2.0
	    U(M-1)=U(M)
	    V(M)=Y/2.0
	    V(M-1)=-V(M)
	  END IF
	  M=M-2
	  IT=0
	  GOTO 10
	END IF
	IF (IT.GE.60) THEN
	  WRITE(*,50)
50	  FORMAT(1X,  'FAIL')
	  JT=0
	  RETURN
	END IF
	IT=IT+1
	DO 60 J=L+2,M
60	A(J,J-2)=0.0
	DO 70 J=L+3,M
70	A(J,J-3)=0.0
	DO 150 K=L,M-1
	  IF (K.NE.L) THEN
	    P=A(K,K-1)
	    Q=A(K+1,K-1)
	    R=0.0
	    IF (K.NE.M-1) R=A(K+2,K-1)
	  ELSE
	    X=A(M,M)+A(M-1,M-1)
	    Y=A(M-1,M-1)*A(M,M)-A(M-1,M)*A(M,M-1)
	    P=A(L,L)*(A(L,L)-X)+A(L,L+1)*A(L+1,L)+Y
	    Q=A(L+1,L)*(A(L,L)+A(L+1,L+1)-X)
	    R=A(L+1,L)*A(L+2,L+1)
	  END IF
	  IF (ABS(P)+ABS(Q)+ABS(R).NE.0.0) THEN
	    XY=1.0
	    IF (P.LT.0.0) XY=-1.0
	    S=XY*SQRT(P*P+Q*Q+R*R)
	    IF (K.NE.L) A(K,K-1)=-S
	    E=-Q/S
	    F=-R/S
	    X=-P/S
	    Y=-X-F*R/(P+S)
	    G=E*R/(P+S)
	    Z=-X-E*Q/(P+S)
	    DO 110 J=K,M
	      P=X*A(K,J)+E*A(K+1,J)
	      Q=E*A(K,J)+Y*A(K+1,J)
	      R=F*A(K,J)+G*A(K+1,J)
	      IF (K.NE.M-1) THEN
	        P=P+F*A(K+2,J)
	        Q=Q+G*A(K+2,J)
	        R=R+Z*A(K+2,J)
	        A(K+2,J)=R
	      END IF
	      A(K+1,J)=Q
	      A(K,J)=P
110	    CONTINUE
	    J=K+3
	    IF (J.GE.M) J=M
	    DO 120 I=L,J
	      P=X*A(I,K)+E*A(I,K+1)
	      Q=E*A(I,K)+Y*A(I,K+1)
	      R=F*A(I,K)+G*A(I,K+1)
	      IF (K.NE.M-1) THEN
	        P=P+F*A(I,K+2)
	        Q=Q+G*A(I,K+2)
	        R=R+Z*A(I,K+2)
	        A(I,K+2)=R
	      END IF
	      A(I,K+1)=Q
	      A(I,K)=P
120	    CONTINUE
	  END IF
150	CONTINUE
	GOTO 10
	END
	
!---------SOLVING THE POLYNOMIAL EQUATION FOR BETA----------
	SUBROUTINE DSRRT(A,BETA,N,M)

	DIMENSION A(N),XR(M),XI(M),B(N)
	IF (ABS(A(1))+1.0.EQ.1.0) THEN
	  L=0
	  WRITE(*,5)
	  RETURN
	END IF
5	FORMAT(1X,'  ERR')
	L=1
	K=M
	IS=0
	W=1.0
	DO 10 I=1,N
10	B(I)=A(I)/A(1)
20	PP=ABS(B(K+1))
	IF (PP.LT.1.0E-6) THEN
	  XR(K)=0.0
	  XI(K)=0.0
	  K=K-1
	  IF (K.EQ.1) THEN
	    XR(K)=-B(2)*W/B(1)
	    XI(K)=0.0
	    RETURN
	  END IF
	  GOTO 20
	END IF
	Q=PP**(1.0/K)
	P=Q
	W=W*P
	DO 30 I=1,K
	  B(I+1)=B(I+1)/Q
	  Q=Q*P
30	CONTINUE
	X=0.0001
	X1=X
	Y=0.2
	Y1=Y
	G=1.0E+37
	DX=1.0
40	U=B(1)
	V=0.0
	DO 50 I=1,K
	  P=U*X1
	  Q=V*Y1
	  PQ=(U+V)*(X1+Y1)
	  U=P-Q+B(I+1)
	  V=PQ-P-Q
50	CONTINUE
	G1=U*U+V*V
	IF (G1.LT.G) GOTO 105
	IF (IS.NE.0) GOTO 80
60	T=T/1.67
	X1=X-T*DX
	Y1=Y-T*DY
	IF (K.GE.50) THEN
	  P=SQRT(X1*X1+Y1*Y1)
	  Q=EXP(85.0/K)
	  IF (P.GE.Q) GOTO 60
	END IF
	IF (T.GE.1.0E-03) GOTO 40
	IF (G.LE.1.0E-5) GOTO 90
65	IS=1
	DD=SQRT(DX*DX+DY*DY)
	IF (DD.GT.1.0) DD=1.0
	DC=6.28/(K+4.5)
70	C=0.0
80	C=C+DC
	DX=DD*COS(C)
	DY=DD*SIN(C)
	X1=X+DX
	Y1=Y+DY
	IF (C.LE.6.29) GOTO 40
	DD=DD/1.67
	IF (DD.GT.1.0E-05) GOTO 70
90	IF (ABS(Y).LE.1.0E-06) THEN
	  P=-X
	  Y=0.0
	  Q=0.0
	ELSE
	  P=-2.0*X
	  Q=X*X+Y*Y
	  XR(K)=X*W
	  XI(K)=-Y*W
	  K=K-1
	END IF
	DO 100 I=1,K
	  B(I+1)=B(I+1)-B(I)*P
	  B(I+2)=B(I+2)-B(I)*Q
100	CONTINUE
	XR(K)=X*W
	XI(K)=Y*W
	K=K-1
	IF (K.EQ.1) THEN
	  XR(K)=-B(2)*W/B(1)
	  XI(K)=0.0
	  RETURN
	END IF
	GOTO 20
105	G=G1
	X=X1
	Y=Y1
	IS=0
	IF (G.LE.1.0E-5) GOTO 90
	U1=K*B(1)
	V1=0.0
	DO 110 I=2,K
	  P=U1*X
	  Q=V1*Y
	  PQ=(U1+V1)*(X+Y)
	  U1=P-Q+(K-I+1)*B(I)
	  V1=PQ-P-Q
110	CONTINUE
	P=U1*U1+V1*V1
	IF (P.LE.1.0E-5) GOTO 65
	DX=(U*U1+V*V1)/P
	DY=(U1*V-V1*U)/P
	T=1.0+4.0/K
	GOTO 60
	
	BETA=0.
	DO J=1,M
	print*,"XR(J)=",XR(J)
	IF(XR(J).GT.BETA) BETA=XR(J)
	END DO
	pause
	
	RETURN
	END
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!particle 
!------------------------------------------------------------
!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

      SUBROUTINE GET_PAR(U,NTIME)
C--------------------------------------------------------------
C  PURPOSE
C  
C  MAIN STEP OF CALCULTION
C
C--------------------------------------------------------------

	INCLUDE 'COMMON1.FI'
        REAL U(0:M1,0:M2,0:M3,3)
	CHARACTER(3) FILENAME
	REAL(8) TIME01,TIME02
	CALL ETIME(TIME01)


	IF(INIT_P.EQ.1) CALL INIT_PAR(NTIME)
	  CALL CALC_PAR(U,NTIME)
	  CALL OUTPUT_PAR(U,NTIME)

       IF(NTIME.LE.20000)FLAGPAR=0
       IF(NTIME.GT.20000)THEN
	    FLAGPAR=FLAGPAR+1
        CALL WRITEAVE_PAR(NTIME) !,FLAGPAR)
        ENDIF

	CALL ETIME(TIME02)
	TIME_PAR=TIME02-TIME01

!    WRITE(*,*) 'THE TIME COSTED IN GET_PAR :', TIME_PAR_MAX
	 WRITE(*,*) 'NUMBER OF PARTICLES        :',    NP_ALL
      RETURN
      END
C
C  ****************************** INIDEM **********************

      SUBROUTINE INIT_PAR(NTIME)
C--------------------------------------------------------------
C  PURPOSE
C  
C  INITIAL PHYSICL PARAMETERS OF PARTICLES
C
C--------------------------------------------------------------
	INCLUDE 'COMMON1.FI'
 
!-----GRID POINT FOR PARTICLES----------------------------------------
        
	INIT_P=0
	NTIMEP=NTIME
        DENS_P=2220
        DIAM_P=0.0002/HHH
        NSTEP_P=1
        PI=3.141592653
	DT_P=DT !DELTAT/REAL(NSTEP_P)
          write(*,*)'DT_P,DT',DT_P,DT
	MASS_P=DENS_P*PI*DIAM_P**3/6.
	TAU_P=DENS_P*DIAM_P**2/18.*RE
       N_PAR=0
       N_CAV=0
      WRITE(*,*)'INIT_PAR'
     
	RETURN
      END
C
C  ****************************** CALU_PAR **********************

      SUBROUTINE CALC_PAR(U,NTIME)
C--------------------------------------------------------------
C  PURPOSE
C  
C  1.ADD NEW PARTICLES
C  2.CALCULATE VELOCITY AND XYZ LOCATION OF PARTICLES
C  3.CALCULATE MAX CFL NUMBER OF PARTICLES TO MAKE SURE PRECISION
C
C--------------------------------------------------------------
	INCLUDE 'COMMON1.FI'
      REAL U(0:M1,0:M2,0:M3,3)

C  ADD NEW PARTICLES
C	IF(MOD(NTIME,NSTEP_P).EQ.0)CALL ADD_PAR(U)
	IF(NP_ALL.LT.10000)CALL ADD_PAR(U)
        NP_ALL=N_PAR 
C	IF(NTIME.LT.50)CALL ADD_PAR(U)
C  CALCULATE THE MAX CFL NUMBER OF ALL PARTICLES
	CALL CALC_CFL
     
C  CALCULATE THE VELOCITY AND LOCATION OF ALL PARTICLES 
C	CALL CALC_UX_EULER(U)
  	CALL CALC_UX_4RK  (U)
	CALL  BC_PAR(NTIME)
      RETURN
      END
C
C  ****************************** NEWDEM ********************
      SUBROUTINE ADD_PAR(U)
C------------------------------------------------------------
C  PURPOSE
C  
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

      REAL    U(0:M1,0:M2,0:M3,3)
	REAL    UF,VF,WF
	REAL    UP,VP,WP
	REAL    XP,YP,ZP
	INTEGER IP,JP,KP

      DOUBLE PRECISION R
	REAL NRND1


	REAL    UTAU_T
	INTEGER NUMP_NEW_FLOW
	REAL    U_ALL , C0

	CALL CPU_TIME(R)



	IF(1.EQ.1)THEN

	DO II=1,10000
C  INITIAL THE PARTICALES POSITION
	  XP=(ALX-1.0*DIAM_P)*NRND1(R)+0.5*DIAM_P
C	  ZP=ZG(2)+5*DIAM_P
	  YP=Y(1)+(0.1-1.0*DIAM_P)*NRND1(R)+0.5*DIAM_P
	  ZP=(ALZ-1.0*DIAM_P)*NRND1(R)+0.5*DIAM_P
	  N_PAR=N_PAR+1
	  IIP=N_PAR
        PAR(1,IIP)=XP
        PAR(2,IIP)=YP
        PAR(3,IIP)=ZP
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  INDEX_PAR(1,IIP)=IP
	  INDEX_PAR(2,IIP)=JP
	  INDEX_PAR(3,IIP)=KP
	
C  INITIAL THE PARTICALES VELOCITY
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)

!	  U_ALL=SQRT(2.*9.8*DIAM_P*2)/SIN(0.62832)
!	  UF=U_ALL*COS(0.62832)
!	  VF=-2.0*SQRT(2.*9.8*DIAM_P*2)
!	  WF=0.0
	  PAR(4,IIP)=UF
	  PAR(5,IIP)=VF
	  PAR(6,IIP)=WF
	ENDDO
	ENDIF

	      
	      
	OPEN(1,FILE=
	1       'RESULT_PAR/ADDPAR.plt')
	  WRITE(1,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "IF" "JF" "KF"'
	  WRITE(1,*) 'ZONE T=" ',N_PAR,'"'
	  DO II=1,N_PAR
	     IIP=II
	     XP=PAR(1,IIP)
	     YP=PAR(2,IIP)
	     ZP=PAR(3,IIP)
	     UP=PAR(4,IIP)
	     VP=PAR(5,IIP)
	     WP=PAR(6,IIP)
	     CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	     CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
	     WRITE(1,400) XP,YP,ZP,UP,VP,WP,IP,JP,KP
	  ENDDO
	  CLOSE(1) 
       NP_ALL=N_PAR 
      ! C0=N_PAR*3.1415926535*Diam_P**3/6./(ALX*ALZ*2)	     
100	FORMAT(A8,I6,A1)
200	FORMAT(A8,I6)
300	FORMAT(6F20.6)
400	FORMAT(6F20.6,3I6)	      
	      

 

      RETURN
      END
C
C  ****************************** CALC_CFL ******************
      SUBROUTINE CALC_CFL
C------------------------------------------------------------
C  PURPOSE
C  
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'
  

	REAL    UP,VP,WP
	REAL    XP,YP,ZP
	INTEGER IP,JP,KP

	REAL    CFLX ,CFLY ,CFLZ
	REAL    CFL_P

	CFLX=0.
	CFLY=0.
	CFLZ=0.
	DO II=1,N_PAR
	  IIP=II
	  XP=PAR(1,IIP)
	  YP=PAR(2,IIP)
	  ZP=PAR(3,IIP)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  UP=ABS(PAR(4,IIP))
	  VP=ABS(PAR(5,IIP))
	  WP=ABS(PAR(6,IIP))
	  CFLX=MAX(CFLX,UP*DT_P*DX1)
	  CFLY=MAX(CFLY,VP*DT_P/DY(JP))
	  CFLZ=MAX(CFLZ,WP*DT_P*DX3)
	ENDDO

	CFL_P=MAX(CFLX,CFLY,CFLZ)
	CFL_P_MAX=CFL_P
C	WRITE(65,*) 'MAX CFL NUMBER OF  PARTICLES: ',CFL_P_MAX,CFL_P
   
	RETURN
      END
C
C  ************************** CALC_UX ***********************
      SUBROUTINE CALC_UX_EULER(U)
C------------------------------------------------------------
C  EULER METHOD
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

      REAL U(0:M1,0:M2,0:M3,3)

	INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UP,VP,WP

	REAL    UF,VF,WF
	REAL    QU,QV,QW

	DO II=1,N_PAR
	  IIP=II
	  UP=PAR(4,IIP)
	  VP=PAR(5,IIP)
	  WP=PAR(6,IIP)
	  XP=PAR(1,IIP)
	  YP=PAR(2,IIP)
	  ZP=PAR(3,IIP)
	  CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
	  CALL FORCE_PAR(UP,VP,WP,UF,VF,WF,QU,QV,QW)
	  XP=PAR(1,IIP)+DT_P*UP
	  YP=PAR(2,IIP)+DT_P*VP
	  ZP=PAR(3,IIP)+DT_P*WP
	  UP=PAR(4,IIP)+DT_P*QU
	  VP=PAR(5,IIP)+DT_P*QV
	  WP=PAR(6,IIP)+DT_P*QW

	  PAR(1,IIP)=XP
	  PAR(2,IIP)=YP
	  PAR(3,IIP)=ZP
	  PAR(4,IIP)=UP
	  PAR(5,IIP)=VP
	  PAR(6,IIP)=WP
      ENDDO


	RETURN
      END
C
C  ************************** CALC_UX ***********************
      SUBROUTINE CALC_UX_4RK(U)
C------------------------------------------------------------
C  4 ORDER RUNGE-KUTTA METHOD
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

      REAL U(0:M1,0:M2,0:M3,3)

	INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UP,VP,WP

	REAL    UF,VF,WF
	REAL    QU,QV,QW
	REAL    QURK,QVRK,QWRK
	REAL    QXRK,QYRK,QZRK

	DO II=1,N_PAR
	  IIP=II

	  QXRK=0.
	  QYRK=0.
	  QZRK=0.
	  QURK=0.
	  QVRK=0.
	  QWRK=0.

	  XP=PAR(1,IIP)
	  YP=PAR(2,IIP)
	  ZP=PAR(3,IIP)
	  UP=PAR(4,IIP)
	  VP=PAR(5,IIP)
	  WP=PAR(6,IIP)
	  IF(YP.LT.0.5*DIAM_P) GOTO 100
	  IF(YP.GT.ALY-0.5*DIAM_P) GOTO 100
C       EAST
	  IF(XP.GT.(ALX-0.5*DIAM_P) )THEN
	   XP=1.0*DIAM_P+XP-ALX
	  ENDIF
C       WEST
	  IF(XP.LT.(0.5*DIAM_P   ))THEN
	   XP=ALX-DIAM_P+XP
        ENDIF
	
C       FRONT
	  IF( ZP.GT.(ALZ-0.5*DIAM_P) )THEN
	    ZP=1.0*DIAM_P+ZP-ALZ
	  ENDIF
C       BACK
	  IF( ZP.LT.(0.5*DIAM_P   ))THEN
	    ZP=ALZ-DIAM_P+ZP
          ENDIF 
	  
	  CALL INDEX_IJK (XP,YP,ZP,IP,JP,KP)
!	WRITE(*,*)'INDEX_IJK1',XP,YP,ZP,IP,JP,KP
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
!	WRITE(*,*)'INTERPOLAT1',IP,JP,KP,UF,VF,WF
 	  CALL FORCE_PAR (UP,VP,WP,UF,VF,WF,QU,QV,QW)
!	WRITE(*,*)'FORCE_PAR1',UP,VP,WP,UF,VF,WF,QU,QV,QW
 	  QXRK=QXRK+UP
	  QYRK=QYRK+VP
	  QZRK=QZRK+WP
	  QURK=QURK+QU
	  QVRK=QVRK+QV
	  QWRK=QWRK+QW





	  XP=PAR(1,IIP)+0.5*DT_P*UP
	  YP=PAR(2,IIP)+0.5*DT_P*VP
	  ZP=PAR(3,IIP)+0.5*DT_P*WP
	  UP=PAR(4,IIP)+0.5*DT_P*QU
	  VP=PAR(5,IIP)+0.5*DT_P*QV
	  WP=PAR(6,IIP)+0.5*DT_P*QW
	  IF(YP.LT.0.5*DIAM_P) GOTO 100
	  IF(YP.GT.ALY-0.5*DIAM_P) GOTO 100
C       EAST
	  IF(XP.GT.(ALX-0.5*DIAM_P) )THEN
	   XP=1.0*DIAM_P+XP-ALX
	  ENDIF
C       WEST
	  IF(XP.LT.(0.5*DIAM_P   ))THEN
	   XP=ALX-DIAM_P+XP
          ENDIF
C       FRONT
	  IF( ZP.GT.(ALZ-0.5*DIAM_P) )THEN
	    ZP=1.0*DIAM_P+ZP-ALZ
	  ENDIF
C       BACK
	  IF( ZP.LT.(0.5*DIAM_P   ))THEN
	    ZP=ALZ-DIAM_P+ZP
          ENDIF 
	  
	  CALL INDEX_IJK (XP,YP,ZP,IP,JP,KP)
!	WRITE(*,*)'INDEX_IJK2',XP,YP,ZP,IP,JP,KP
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
!	WRITE(*,*)'INTERPOLAT2',IP,JP,KP,UF,VF,WF
	  CALL FORCE_PAR (UP,VP,WP,UF,VF,WF,QU,QV,QW)
!	WRITE(*,*)'FORCE_PAR2',UP,VP,WP,UF,VF,WF,QU,QV,QW
	  QXRK=QXRK+2*UP
	  QYRK=QYRK+2*VP
	  QZRK=QZRK+2*WP
	  QURK=QURK+2*QU
	  QVRK=QVRK+2*QV
	  QWRK=QWRK+2*QW
	  XP=PAR(1,IIP)+0.5*DT_P*UP
	  YP=PAR(2,IIP)+0.5*DT_P*VP
	  ZP=PAR(3,IIP)+0.5*DT_P*WP
	  UP=PAR(4,IIP)+0.5*DT_P*QU
	  VP=PAR(5,IIP)+0.5*DT_P*QV
	  WP=PAR(6,IIP)+0.5*DT_P*QW






	  IF(YP.LT.0.5*DIAM_P) GOTO 100
	  IF(YP.GT.ALY-0.5*DIAM_P) GOTO 100
C       EAST
	  IF(XP.GT.(ALX-0.5*DIAM_P) )THEN
	   XP=1.0*DIAM_P+XP-ALX
	  ENDIF
C       WEST
	  IF(XP.LT.(0.5*DIAM_P   ))THEN
	   XP=ALX-DIAM_P+XP
          ENDIF
C       FRONT
	  IF( ZP.GT.(ALZ-0.5*DIAM_P) )THEN
	    ZP=1.0*DIAM_P+ZP-ALZ
	  ENDIF
C       BACK
	  IF( ZP.LT.(0.5*DIAM_P   ))THEN
	    ZP=ALZ-DIAM_P+ZP
          ENDIF 
	  CALL INDEX_IJK (XP,YP,ZP,IP,JP,KP)
!	WRITE(*,*)'INDEX_IJK3',XP,YP,ZP,IP,JP,KP
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
!	WRITE(*,*)'INTERPOLAT3',IP,JP,KP,UF,VF,WF
	  CALL FORCE_PAR (UP,VP,WP,UF,VF,WF,QU,QV,QW)
!	WRITE(*,*)'FORCE_PAR3',UP,VP,WP,UF,VF,WF,QU,QV,QW
	  QXRK=QXRK+2*UP
	  QYRK=QYRK+2*VP
	  QZRK=QZRK+2*WP
	  QURK=QURK+2*QU
	  QVRK=QVRK+2*QV
	  QWRK=QWRK+2*QW

	  XP=PAR(1,IIP)+DT_P*UP
	  YP=PAR(2,IIP)+DT_P*VP
	  ZP=PAR(3,IIP)+DT_P*WP
	  UP=PAR(4,IIP)+DT_P*QU
	  VP=PAR(5,IIP)+DT_P*QV
	  WP=PAR(6,IIP)+DT_P*QW






	  IF(YP.LT.0.5*DIAM_P) GOTO 100
  	  IF(YP.GT.ALY-0.5*DIAM_P) GOTO 100
C       EAST
	  IF(XP.GT.(ALX-0.5*DIAM_P) )THEN
	   XP=1.*DIAM_P+XP-ALX
	  ENDIF
C       WEST
	  IF(XP.LT.(0.5*DIAM_P   ))THEN
	   XP=ALX-DIAM_P+XP
          ENDIF
C       FRONT
	  IF( ZP.GT.(ALZ-0.5*DIAM_P) )THEN
	    ZP=1.*DIAM_P+ZP-ALZ
	  ENDIF
C       BACK
	  IF( ZP.LT.(0.5*DIAM_P   ))THEN
	    ZP=ALZ-DIAM_P+ZP
         ENDIF 
	  CALL INDEX_IJK (XP,YP,ZP,IP,JP,KP)
!	WRITE(*,*)'INDEX_IJK4',XP,YP,ZP,IP,JP,KP
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
!	WRITE(*,*)'INTERPOLAT4',IP,JP,KP,UF,VF,WF
	  CALL FORCE_PAR (UP,VP,WP,UF,VF,WF,QU,QV,QW)
!	WRITE(*,*)'FORCE_PAR4',UP,VP,WP,UF,VF,WF,QU,QV,QW
	  QXRK=QXRK+UP
	  QYRK=QYRK+VP
	  QZRK=QZRK+WP
	  QURK=QURK+QU
	  QVRK=QVRK+QV
	  QWRK=QWRK+QW

	  XP=PAR(1,IIP)+DT_P*QXRK/6.
	  YP=PAR(2,IIP)+DT_P*QYRK/6.
	  ZP=PAR(3,IIP)+DT_P*QZRK/6.
	  UP=PAR(4,IIP)+DT_P*QURK/6.
	  VP=PAR(5,IIP)+DT_P*QVRK/6.
	  WP=PAR(6,IIP)+DT_P*QWRK/6.
	  GOTO 200


100	  UP=PAR(4,IIP)
	  VP=PAR(5,IIP)
	  WP=PAR(6,IIP)
	  XP=PAR(1,IIP)
	  YP=PAR(2,IIP)
	  ZP=PAR(3,IIP)


	  CALL INDEX_IJK (XP,YP,ZP,IP,JP,KP)
	  CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
	  CALL FORCE_PAR (UP,VP,WP,UF,VF,WF,QU,QV,QW)
	  XP=PAR(1,IIP)+DT_P*UP
	  YP=PAR(2,IIP)+DT_P*VP
	  ZP=PAR(3,IIP)+DT_P*WP
	  UP=PAR(4,IIP)+DT_P*QU
	  VP=PAR(5,IIP)+DT_P*QV
	  WP=PAR(6,IIP)+DT_P*QW

200	  PAR(1,IIP)=XP
	  PAR(2,IIP)=YP
	  PAR(3,IIP)=ZP
	  PAR(4,IIP)=UP
	  PAR(5,IIP)=VP
	  PAR(6,IIP)=WP
      ENDDO


	RETURN
      END
C
C  ************************ FORCE_PAR ***********************
      SUBROUTINE FORCE_PAR(UP,VP,WP,UF,VF,WF,QU,QV,QW)
C------------------------------------------------------------
C  PURPOSE
C  
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

	REAL  UP,VP,WP
	REAL  UF,VF,WF
	REAL  QU,QV,QW
	REAL  VISCOS
		VISCOS=0.000015 !1.5e-5
		
	REP =SQRT((UF-UP)**2+(VF-VP)**2+(WF-WP)**2)*DIAM_P*RE
C-----STOKES 1851 DRAG FORCE
C	FCOR=1.0
C-----SCHILLER AND NAUMANN,1933 DRAG FORCE
C	FCOR=1.+0.15*REP**0.687
C-----CLIFT AND GAUVIN,1970 DRAG FORCE
	FCOR=1.+0.15*REP**0.687+0.0175*(1.+42500*REP**(-1.16))**(-1)

	QU=FCOR*(UF-UP)/TAU_P
	! QV=FCOR*(VF-VP)/TAU_P !-9.8
	QV=-9.8/(VISCOS*RE/HHH)**2 !-9.8*LY**3/(VISCOS*RE)**2/8.0
	QV=FCOR*(VF-VP)/TAU_P +QV
	QW=FCOR*(WF-WP)/TAU_P
	! WRITE(*,*)QU,QV,QW
	RETURN
      END
C
C  ************************** BC_PAR ************************
      SUBROUTINE BC_PAR(NTIME)
C------------------------------------------------------------
C  PURPOSE
C  
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

	REAL    UP,VP,WP
	REAL    XP,YP,ZP
	INTEGER IIP

	REAL    SEED
	INTEGER ISEED

	CHARACTER(3) FILENAME01


	CALL CPU_TIME(SEED)
	ISEED=-NINT(RAN2(SEED)*1.E08)

	NP_IMPACT =0
	NP_EJECT  =0
	N_CAV     =0

	DO II=1,N_PAR
	  IIP=II
	  
	  

C       EAST
	  IF( PAR(1,IIP).GT.(ALX-0.5*DIAM_P) )THEN
	    PAR(1,IIP)=1.0*DIAM_P+PAR(1,IIP)-ALX
	  ENDIF
C       WEST
	  IF( PAR(1,IIP).LT.(0.5*DIAM_P   ))THEN
	    PAR(1,IIP)=ALX-DIAM_P+PAR(1,IIP)
          ENDIF
	
C       FRONT
	  IF( PAR(3,IIP).GT.(ALZ-0.5*DIAM_P) )THEN
	    PAR(3,IIP)=1.0*DIAM_P+PAR(3,IIP)-ALZ
	  ENDIF
C       BACK
	  IF( PAR(3,IIP).LT.(0.5*DIAM_P   ))THEN
	    PAR(3,IIP)=ALZ-DIAM_P+PAR(3,IIP)
          ENDIF 

C       TOP
	  IF( PAR(2,IIP).GT.(ALY-0.5*DIAM_P) .AND. PAR(5,IIP).GE.0.0 )THEN
	    PAR(5,IIP)=-PAR(5,IIP)
	    PAR(2,IIP)=2.*ALY-DIAM_P-PAR(2,IIP)
	  ENDIF

C       BOTTOM
!	  IF( PAR(2,IIP).LT.(0.5*DIAM_P) .AND. PAR(5,IIP).LE.0.0 )THEN
!!!!!!!!!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	  IF( PAR(2,IIP).LT.(0.5*DIAM_P).AND. PAR(5,IIP).LE.0.0 )THEN
	    IF(1.EQ.0)THEN
	    !PAR(5,IIP)=-PAR(5,IIP)
	    PAR(2,IIP)=DIAM_P-PAR(2,IIP)
	      UP=PAR(4,IIP)
	      VP=PAR(5,IIP)
 	      U_ALL=SQRT(UP**2+VP**2)
	      PAR(4,IIP)=U_ALL*0.8
	      PAR(5,IIP)=U_ALL*0.6
	    !set the up velocity angel is 37 degree     
             ENDIF

	    NP_IMPACT=NP_IMPACT+1

!	    IF(TIMET.GT.300.0 .AND. TIMET.LE.310.0)THEN
!		WRITE(FILENAME01,'(I3)') 100
!	    OPEN(1,FILE=
!     >       '../RESULT_PAR/IMPACT_'//FILENAME01 ,ACCESS='APPEND')
!	    WRITE(1,'(6F12.6)') (PAR(L,IIP),L=1,6)
!	    CLOSE(1)
 !             ENDIF


       	 	CALL SPLASH(IIP,ISEED)
	     IF(1.EQ.0)THEN
	      YP=PAR(2,IIP)
	      UP=PAR(4,IIP)
	      VP=PAR(5,IIP)
	      WP=PAR(6,IIP)
	      U_ALL=(1.0-0.50)*SQRT(UP**2+VP**2+WP**2)
	      PAR(4,IIP)=U_ALL*COS(0.62832)*UP/SQRT(UP**2+WP**2)
	      PAR(5,IIP)=U_ALL*SIN(0.62832)
	      PAR(6,IIP)=U_ALL*COS(0.62832)*WP/SQRT(UP**2+WP**2)
	      PAR(2,IIP)=DIAM_P-YP
		ENDIF
	  ENDIF

	ENDDO
	N_PAR=N_PAR+NP_EJECT

	IF(1.EQ.1)THEN
	  OPEN (1,FILE='RESULT_PAR/SPLASH.plt',ACCESS='APPEND')
	  WRITE(1,*) NTIME,NP_IMPACT,NP_EJECT,N_CAV
	  CLOSE(1)
	ENDIF

C  UPDATE LINKED LIST
      IF(1.EQ.1)THEN
	DO II=N_CAV,1,-1
	  ID_CAV=INDEX_CAV(II)
	  ID_PAR=N_PAR
	  PAR(1,ID_CAV)=PAR(1,ID_PAR)
	  PAR(2,ID_CAV)=PAR(2,ID_PAR)
	  PAR(3,ID_CAV)=PAR(3,ID_PAR)
	  PAR(4,ID_CAV)=PAR(4,ID_PAR)
	  PAR(5,ID_CAV)=PAR(5,ID_PAR)
	  PAR(6,ID_CAV)=PAR(6,ID_PAR)
	  N_PAR=N_PAR-1
	ENDDO
      ENDIF
      

	RETURN
      END
C
C  ************************** BC_PAR ************************
      SUBROUTINE SPLASH(IIP,ISEED)
C------------------------------------------------------------
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

	INTEGER IIP,ISEED

	REAL    UP,VP,WP
	REAL    XP,YP,ZP

	REAL    V_IMP,A_IMP,MEAN
	REAL    P_REB,V_REB,A_REB
	REAL    UP_REB,VP_REB,WP_REB
	INTEGER N_EJE
	REAL    V_EJE,A_EJE
	REAL    UP_EJE,VP_EJE,WP_EJE
	REAL    U_ALL

	REAL    RAND_UNI,RAND_EXP,RAND_NOR
	REAL    PI
	XP=PAR(1,IIP)
	YP=PAR(2,IIP)
	ZP=PAR(3,IIP)
	UP=PAR(4,IIP)
	VP=PAR(5,IIP)
	WP=PAR(6,IIP)
	PI=3.14159265357
	
	V_IMP=SQRT(UP**2+VP**2+WP**2)
	A_IMP=ATAN(ABS(VP)/SQRT(UP**2+WP**2))

C  REBOUD PARTICLE
C  PROBABILITY
	P_REB=0.95*(1.-EXP(-2.*V_IMP*RE/100000.0*1.5/HHH))
	RAND_UNI=RAN1(ISEED)
	IF(P_REB.GT.RAND_UNI)THEN
C  ANGLE
	  RAND_NOR=GASDEV(ISEED)
	  VARIANCE=15
	  MEAN    =30
	  RAND_NOR=RAND_NOR*(VARIANCE)+MEAN
!	  RAND_NOR=RAND_NOR*SQRT(VARIANCE)+MEAN
	  A_REB   =RAND_NOR
      
C  VELOCITY
100	  RAND_NOR=GASDEV(ISEED)
	  MEAN    =0.6 *V_IMP
	  VARIANCE=0.25*V_IMP !  VARIANCE=0.25*V_IMP
	  RAND_NOR=RAND_NOR*VARIANCE+MEAN
	!  RAND_NOR=RAND_NOR*SQRT(VARIANCE)+MEAN
	  V_REB   =RAND_NOR
    !           WRITE(*,*)V_REB
	  IF(V_REB.LE.0.0) GOTO 100 
	!  WRITE(*,*)V_IMP,V_REB
	  A_REB   =PI/180*A_REB
	  U_ALL   =SQRT(UP**2+WP**2)
	  UP_REB  =V_REB*COS(A_REB)*(UP/U_ALL)
	  VP_REB  =V_REB*SIN(A_REB)
	  WP_REB  =V_REB*COS(A_REB)*(WP/U_ALL)
  
	  PAR(4,IIP)=UP_REB
	  PAR(5,IIP)=VP_REB
	  PAR(6,IIP)=WP_REB
C  LOCATION
	!   PAR(2,IIP)=0.001-YP
	 PAR(2,IIP)=DIAM_P-YP
	ELSE
	  N_CAV=N_CAV+1
	  INDEX_CAV(N_CAV)=IIP
	ENDIF



C  EJECTION PARTICLE
C  NUMBER
	!N_EJE=NINT(0.01/SQRT(9.8*DIAM_P)*V_IMP)
	
	!N_EJE=NINT(0.01/SQRT(9.8*DIAM_P)*V_IMP*RE/100000.0*1.5)  !WULIANGGANG
	N_EJE=NINT(0.01/SQRT(9.8*DIAM_P*HHH)*V_IMP*RE/100000.0*1.5/HHH)  !WULIANGGANG	
	IF(N_EJE.GT.0)THEN
	  DO II=N_PAR+NP_EJECT+1,N_PAR+NP_EJECT+N_EJE

C  ANGLE
300	     RAND_NOR=GASDEV(ISEED)
	     VARIANCE=15
	     MEAN    =60
	     RAND_NOR=RAND_NOR*VARIANCE+MEAN
!	     RAND_NOR=RAND_NOR*SQRT(VARIANCE)+MEAN
	     A_EJE   =RAND_NOR
	     IF(A_EJE.LE.0.) GOTO 300
C  VELOCITY
200	     RAND_EXP=EXPDEV(ISEED)  
	     MEAN    =0.08*V_IMP   ! MEAN=0.08*V_IMP
	    ! RAND_EXP=RAND_EXP+MEAN
	     V_EJE   =RAND_EXP*MEAN
            
	     IF(V_EJE.LE.0.) GOTO 200
          !      write(*,*)'V_eje',A_EJE,A_IMP,RAND_EXP



	     A_EJE   =PI/180*A_EJE
	     U_ALL   =SQRT(UP**2+WP**2)
	     UP_EJE  =V_EJE*COS(A_EJE)*(UP/U_ALL)
	     VP_EJE  =V_EJE*SIN(A_EJE)
	     WP_EJE  =V_EJE*COS(A_EJE)*(WP/U_ALL)

	     PAR(4,II)=UP_EJE
	     PAR(5,II)=VP_EJE
	     PAR(6,II)=WP_EJE
C  LOCATION
	     PAR(1,II)=XP
	 !     PAR(2,II)=0.001-YP
	    PAR(2,II)=DIAM_P-YP
	     PAR(3,II)=ZP
	  ENDDO
	  NP_EJECT=NP_EJECT+N_EJE
	ENDIF

	RETURN
      END
C
C  ****************************** PASS_PAR ******************    
C  ****************************** PASS_PAR ******************
C
C  ****************************** INTPL **********************
      SUBROUTINE INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
C------------------------------------------------------------
C  PURPOSE
C  
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'
      COMMON/INDX/IPA(M1),IMA(M1),KPA(M3),KMA(M3)
      COMMON/INDX2/JPA(M2),JMU(M2),JMV(M2)

	INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UF,VF,WF
        REAL    U(0:M1,0:M2,0:M3,3)

	INTEGER IUM,IUP
	INTEGER JUM,JUP
	INTEGER KUM,KUP

	REAL XUM,XUP
	REAL YUM,YUP
	REAL ZUM,ZUP

	REAL PHA,PHB,PHC,PHD,PHE,PHF,PHG,PHH
	REAL PHAB,PHCD,PHEF,PHGH
	REAL PHABCD,PHEFGH
      REAL I3
      
      DO I3=1,3
      IF(I3.EQ.1)THEN
        IUM=IP
	  IUP=IPA(IP)
	  XUM=XG(IP)
	  XUP=XG(IUP)
	  IF(IP.EQ.N1M) XUP=XG(N1)
      ELSE
	IF( XP.GE.XC(IP) )THEN
	  IUM=IP
	  IUP=IPA(IP)
	  XUM=XC(IP)
	  XUP=XC(IUP)
	   IF(IP.EQ.N1M) XUP=XC(N1M)+1./DX1
	ELSE
	  IUM=IMA(IP)
	  IUP=IP
	  XUM=XC(IUM)
	  XUP=XC(IP)	
	  IF(IP.EQ.1) XUM=XC(1)-1./DX1  
	ENDIF
      ENDIF
      
      IF(I3.EQ.2)THEN
      
        JUM=JP
	  JUP=JPA(JP)
	  YUM=YG(JP)
	  YUP=YG(JUP)
		  IF(JP.EQ.N2M)YUP=YC(N2)  
      ELSE
	IF( YP.GE.YC(JP) )THEN
	  JUM=JP
	  JUP=JPA(JP)
	  YUM=YC(JP)
	  YUP=YC(JUP)
	  IF(JP.EQ.N2M)YUP=YC(N2)
	ELSE
	  JUM=JMU(JP)
	  JUP=JP
	  YUM=YC(JUM)
	  YUP=YC(JP)
	 IF(JP.EQ.1)YUM=YG(1)
	ENDIF
	ENDIF
	
	IF(I3.EQ.3)THEN
	  KUM=KP
	  KUP=KPA(KP)
	  ZUM=ZG(KP)
	  ZUP=ZG(KUP)
	   IF(KP.EQ.N3M) ZUP=ZG(N3)  
	ELSE
	IF( ZP.GE.ZC(KP) )THEN
	  KUM=KP
	  KUP=KPA(KP)
	  ZUM=ZC(KP)
	  ZUP=ZC(KUP)
	   IF(KP.EQ.N3M) ZUP=ZC(N3M)+1./DX3	 	  
	ELSE
	  KUM=KMA(KP)
	  KUP=KP
	  ZUM=ZC(KUM)
	  ZUP=ZC(KP)
	  	   IF(KP.EQ.1) ZUM=ZC(1)-1./DX3	 
	ENDIF
	ENDIF

      IF(I3.EQ.1)THEN
	PHA=U(IUM,JUM,KUM,1)
	PHB=U(IUP,JUM,KUM,1)
	PHC=U(IUP,JUP,KUM,1)
	PHD=U(IUM,JUP,KUM,1)
	PHE=U(IUM,JUM,KUP,1)
	PHF=U(IUP,JUM,KUP,1)
	PHG=U(IUP,JUP,KUP,1)
	PHH=U(IUM,JUP,KUP,1)
	
	PHAB=(PHB-PHA)/(XUP-XUM)*(XP-XUM)+PHA
	PHCD=(PHC-PHD)/(XUP-XUM)*(XP-XUM)+PHD
	PHEF=(PHF-PHE)/(XUP-XUM)*(XP-XUM)+PHE
	PHGH=(PHG-PHH)/(XUP-XUM)*(XP-XUM)+PHH
	PHABCD=(PHCD-PHAB)/(YUP-YUM)*(YP-YUM)+PHAB
	PHEFGH=(PHGH-PHEF)/(YUP-YUM)*(YP-YUM)+PHEF
	
	UF=(PHEFGH-PHABCD)/(ZUP-ZUM)*(ZP-ZUM)+PHABCD
	ENDIF
	
	
	IF(I3.EQ.2)THEN

	PHA=U(IUM,JUM,KUM,2)
	PHB=U(IUP,JUM,KUM,2)
	PHC=U(IUP,JUP,KUM,2)
	PHD=U(IUM,JUP,KUM,2)
	PHE=U(IUM,JUM,KUP,2)
	PHF=U(IUP,JUM,KUP,2)
	PHG=U(IUP,JUP,KUP,2)
	PHH=U(IUM,JUP,KUP,2)
	
	PHAB=(PHB-PHA)/(XUP-XUM)*(XP-XUM)+PHA
	PHCD=(PHC-PHD)/(XUP-XUM)*(XP-XUM)+PHD
	PHEF=(PHF-PHE)/(XUP-XUM)*(XP-XUM)+PHE
	PHGH=(PHG-PHH)/(XUP-XUM)*(XP-XUM)+PHH
	PHABCD=(PHCD-PHAB)/(YUP-YUM)*(YP-YUM)+PHAB
	PHEFGH=(PHGH-PHEF)/(YUP-YUM)*(YP-YUM)+PHEF
	VF=(PHEFGH-PHABCD)/(ZUP-ZUM)*(ZP-ZUM)+PHABCD
      ENDIF
      
      IF(I3.EQ.3)THEN
	PHA=U(IUM,JUM,KUM,3)
	PHB=U(IUP,JUM,KUM,3)
	PHC=U(IUP,JUP,KUM,3)
	PHD=U(IUM,JUP,KUM,3)
	PHE=U(IUM,JUM,KUP,3)
	PHF=U(IUP,JUM,KUP,3)
	PHG=U(IUP,JUP,KUP,3)
	PHH=U(IUM,JUP,KUP,3)
	PHAB=(PHB-PHA)/(XUP-XUM)*(XP-XUM)+PHA
	PHCD=(PHC-PHD)/(XUP-XUM)*(XP-XUM)+PHD
	PHEF=(PHF-PHE)/(XUP-XUM)*(XP-XUM)+PHE
	PHGH=(PHG-PHH)/(XUP-XUM)*(XP-XUM)+PHH
	PHABCD=(PHCD-PHAB)/(YUP-YUM)*(YP-YUM)+PHAB
	PHEFGH=(PHGH-PHEF)/(YUP-YUM)*(YP-YUM)+PHEF
	WF=(PHEFGH-PHABCD)/(ZUP-ZUM)*(ZP-ZUM)+PHABCD
	ENDIF
	ENDDO

  
	RETURN
      END
C
C  ----------------------------- INDEX_IJK ------------------
      SUBROUTINE INDEX_IJK (XP,YP,ZP,IP,JP,KP)
C------------------------------------------------------------
C  PURPOSE
C  
C
C------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

	REAL    XP,YP,ZP
	INTEGER IP,JP,KP


	IF(XP.LT.XC(1  )) XP=ANINT(XP*1E06+0.5)/1E06
	IF(XP.GT.XC(N1M)) XP=ANINT(XP*1E06-0.5)/1E06
	IF(YP.LT.YC(1  )) YP=ANINT(YP*1E06+0.5)/1E06
	IF(YP.GT.YC(N2M)) YP=ANINT(YP*1E06-0.5)/1E06
	IF(ZP.LT.ZC(1  )) ZP=ANINT(ZP*1E06+0.5)/1E06
	IF(ZP.GT.ZC(N3M)) ZP=ANINT(ZP*1E06-0.5)/1E06
       !  WRITE(*,*)M1,M2,M3
        ! WRITE(*,*)N1M,N2M,N3M,XC(N1M),YC(N2M),ZC(N3M)
	IP=-1
	DO I=1,N1M
	  IF( XP.GE.XG(I) .AND. XP.LE.XG(I+1) )THEN
	    IP=I
	  ENDIF
	  IF( IP.NE.-1)EXIT
      ENDDO

	JP=-1
	DO J=1,N2M
	  IF( YP.GE.Y(J) .AND. YP.LE.Y(J+1))THEN
	    JP=J
	  ENDIF
	  IF( JP.NE.-1)EXIT
      ENDDO

	KP=-1
	DO K=1,N3M
	  IF( ZP.GE.ZG(K) .AND. ZP.LE.ZG(K+1))THEN
	    KP=K
	  ENDIF
	  IF( KP.NE.-1)EXIT
      ENDDO

	RETURN
      END
C

C
C  *************************** POST_PAR ***********************
    
C  *************************** POST_PAR ***********************
   
C
C  ************************** OUTPUT_PAR **********************
      SUBROUTINE OUTPUT_PAR(U,NTIME)
C--------------------------------------------------------------
C  PURPOSE
C  
C  
C
C--------------------------------------------------------------
	INCLUDE 'COMMON1.FI'

      REAL    U(0:M1,0:M2,0:M3,3)

	! CHARACTER(3) FILENAME
	CHARACTER(6) FILENAME

	INTEGER IP,JP,KP
	REAL    XP,YP,ZP
	REAL    UP,VP,WP
	REAL    TERMP
	 NPRINT_P=4000
	 
 	IF(1.EQ.1)THEN
	 !  UTAU_TEMP=UTAU_WALL(1)
	  UTAU_TEMP=SQRT(ABS(WSM(1)))
	  ST=TAU_P*UTAU_TEMP**2*RE
	  OPEN(1,FILE='RESULT_PAR/PAR_T.plt',ACCESS='APPEND')
	!  WRITE(1,'(3F15.6,3I10)') NTIMET,CFL_P_MAX,ST,
   !  >                           NP_ALL,N_MAX,NP_MAX_BC
	  WRITE(1,*) NTIME, NP_ALL
	  CLOSE(1)
	ENDIF




	IF( NTIME.GE.100)THEN
	IF(MOD(NTIME,NPRINT_P).EQ.0)THEN

!	  WRITE(FILENAME,'(I6)') 100
	  WRITE(FILENAME,'(I6)') NTIME+100000

C	  OPEN(1,FILE=
C	1       '../RESULT_PAR/PAR'//FILENAME01//'_'//FILENAME//'.plt')
	  OPEN(115,FILE=
     >       'RESULT_PAR/PAR_'//FILENAME//'.plt')
	  WRITE(115,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "UF" "VF" "WF"'
	  WRITE(115,100) 'ZONE T="',NTIME,'"'
!	  WRITE(115,200) 'I= ',N_PAR
        TERMP=RE*1.5/100000.0/HHH
	  DO II=1,N_PAR
	     IIP=II
	     XP=PAR(1,IIP)
	     YP=PAR(2,IIP)
	     ZP=PAR(3,IIP)
	     UP=PAR(4,IIP)*TERMP
	     VP=PAR(5,IIP)*TERMP
	     WP=PAR(6,IIP)*TERMP
	     CALL INDEX_IJK(XP,YP,ZP,IP,JP,KP)
	     CALL INTERPOLAT(XP,YP,ZP,IP,JP,KP,UF,VF,WF,U)
	     UF=UF*TERMP
	     VF=VF*TERMP
	     WF=WF*TERMP
	     
	     WRITE(115,400) XP,YP,ZP,UP,VP,WP,UF,VF,WF
	  ENDDO
	  CLOSE(115)

	ENDIF
	ENDIF


100	FORMAT(A8,I6,A1)
200	FORMAT(A3,I6)
300	FORMAT(6F20.6)
400	FORMAT(9F20.6)

      RETURN
      END


C***********************************************************************
C     PORTABLE RADOM NUMBER GENERATORS
C     FUNCTION RAN1   --- UNIFORM DEVIATE BETWEEN 0.0 AND 1.0
C     FUNCTION EXPDEV --- EXPONENTIALLY DISTRIBUTED POSITIVE RANDOM 
C                         DEVIATE OF UNIT MEAN
C     FUNCTION GASDEV --- NORMALLY DISTRIBUTED DEVIATE WITH ZERO MEAN 
C                         AND UNIT VARIANCE
C                         RANDOM*SQRT(VARIANCE)+MEAN
C     FUNCTION RAN2   --- UNIFORM DEVIATE BETWEEN 0.0 AND 1.0
C***********************************************************************
      FUNCTION RAN1(IDUM)
C
C***********************************************************************
      INTEGER IDUM,IA,IM,IQ,IR,NTAB,NDIV
      REAL RAN1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER J,K,IV(NTAB),IY
      SAVE IV,IY
      DATA IV /NTAB*0/, IY /0/
      IF (IDUM.LE.0.OR.IY.EQ.0) THEN
        IDUM=MAX(-IDUM,1)
        DO 11 J=NTAB+8,1,-1
          K=IDUM/IQ
          IDUM=IA*(IDUM-K*IQ)-IR*K
          IF (IDUM.LT.0) IDUM=IDUM+IM
          IF (J.LE.NTAB) IV(J)=IDUM
11      CONTINUE
        IY=IV(1)
      ENDIF
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF (IDUM.LT.0) IDUM=IDUM+IM
      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=MIN(AM*IY,RNMX)
      RETURN
      END
C
C***********************************************************************
      FUNCTION EXPDEV(IDUM)
C
C***********************************************************************
      INTEGER IDUM
      REAL EXPDEV
CU    USES RAN1
      REAL DUM,RAN1
1     DUM=RAN1(IDUM)
      IF(DUM.EQ.1.0)GOTO 1
      EXPDEV=-LOG(1.-DUM)
      RETURN
      END
C
C***********************************************************************
      FUNCTION GASDEV(IDUM)
C
C***********************************************************************
      INTEGER IDUM
      REAL GASDEV
CU    USES RAN1
      INTEGER ISET
      REAL FAC,GSET,RSQ,V1,V2,RAN1
      SAVE ISET,GSET
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN1(IDUM)-1.
        V2=2.*RAN1(IDUM)-1.
        RSQ=V1**2+V2**2
        IF(RSQ.GE.1..OR.RSQ.EQ.0.)GOTO 1
        FAC=SQRT(-2.*LOG(RSQ)/RSQ)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END
C
C***********************************************************************
	FUNCTION RAN2(R) 
C***********************************************************************
	REAL S, U, V, R
	INTEGER M
      S=65536.0
      U=2053.0
      V=13849.0
      M=R/S
      R=R-M*S
      R=U*R+V
      M=R/S
      R=R-M*S
      RAN2=R/S
      RETURN
      END
C
C***********************************************************************
      REAL FUNCTION NRND1(R) 
C***********************************************************************
	DOUBLE PRECISION S, U, V, R
	INTEGER M
      S=65536.0
      U=2053.0
      V=13849.0
      M=R/S
      R=R-M*S
      R=U*R+V
      M=R/S
      R=R-M*S
      NRND1=R/S
      RETURN
      END
C
      SUBROUTINE WRITEAVE_PAR(NTIME) !,FLAGPAR)
       INCLUDE 'COMMON1.FI'
       ! INTEGER FLAGPAR
	Integer N_par,NP_all
       
	Parameter (Max_zone=500000)
	Integer   Num_par_zone(Max_zone)

	REAL      Num_par_K(M2)
	REAL      Num_par(M1,M2,M3)
	REAL      Q (M1,M2,M3)
	REAL      WX(M1,M2,M3)
	REAL      C (M1,M2,M3)


	Real      Cm (M2), Crms (M2)
	Real      Um (M2), Urms (M2)
	Real      Vm (M2), Vrms (M2)
	Real      Wm (M2), Wrms (M2)
	Real      Qm (M2), Qrms (M2)
	Real      Wxm(M2), Wxrms(M2)
	Real      UPm(M2), UPrms(M2)
	Real      VPm(M2), VPrms(M2)
	Real      WPm(M2), WPrms(M2)
	
	Real     UV(M2), VW(M2), UW(M2)



	
	REAL   UTAU,RE1,RE2

	CHARACTER(3) FILENAME
	CHARACTER(6) FILENAME1
	CHARACTER(1) FILENAME2

       
       QQ=600.0

C	WRITE(*,*) '--------- RE ------------'                
	RE1=SQRT(WSM(1))*RE
	RE2=RE*1.5/100000.0/HHH   !/RE1
	NPRINT_P=2000
	
C-----Calculate Means----------------------------------------------------
C-----INIT Means VABIS----------------------------------------------------
        IF(FLAGPAR.EQ.1)THEN
	 ! Index_ZONE=1
	Do I=1,Max_zone
	   Num_par_zone(I)=0
	Enddo
	Do K=1,M2
	UPrms1(K)=0.0
        VPrms1(K)=0.0
        WPrms1(K)=0.0
	UPVP1(K)=0.0
        VPWP1(K)=0.0
        UPWP1(K)=0.0
	Num_par_K1(K)=0
	Enddo
	ENDIF

	  Do K=1,M3
	  Do J=1,M2
	  Do I=1,M1
	     Num_par(I,J,K)=0
	           C(I,J,K)=0.
	  Enddo
	  Enddo
	  Enddo


	  Do II=1,N_PAR
             IIP=II
	     Num_par_zone(FLAGPAR)=Num_par_zone(FLAGPAR)+1
	     XP=Par(1,IIP)
	     YP=Par(2,IIP)
	     ZP=Par(3,IIP)
	     UP=Par(4,IIP)
	     VP=Par(5,IIP)
	     WP=Par(6,IIP)
	     Call Index_IJK (XP,YP,ZP,IP,JP,KP)

	     Num_par_K1(JP)=Num_par_K1(JP)+1
	     Num_par(IP,JP,KP)=Num_par(IP,JP,KP)+1
	     UPrms1(JP) = UPrms1(JP)+UP*UP
	     VPrms1(JP) = VPrms1(JP)+VP*VP
	     WPrms1(JP) = WPrms1(JP)+WP*WP
	     UPm1(JP) = UPm1(JP)+UP
	     VPm1(JP) = VPm1(JP)+VP
	     WPm1(JP) = WPm1(JP)+WP
	     UPVP1(JP)  =UPVP1(JP) + UP*VP
	     VPWP1(JP)  =VPWP1(JP) + VP*WP
	     UPWP1(JP)  =UPWP1(JP) + UP*WP    
	  End do 
	  Do K=1,N3M
	  Do J=1,N2M
	  Do I=1,N1M
	     C(I,J,K)=3.1415926535*Diam_P**3/6.*Num_par(I,J,K)*QQ

	1                            /(XG(I+1)-XG(I))
	1                            /(YG(J+1)-YG(J))
	1                            /(ZG(K+1)-ZG(K))
          Cm1(J)=Cm1(J)+ C(I,J,K)
	  Enddo
	  Enddo
	  Enddo

	IF(MOD(NTIME,NPRINT_P).EQ.O)THEN
	Do K=1,N2M

	
	   If(Num_par_K1(K).EQ.0)Then
	      UPm(K)  = 0.
	      VPm(K)  = 0.
	      WPm(K)  = 0.
	      UPrms(K)  = 0.
	      VPrms(K)  = 0.
	      WPrms(K)  = 0.	  
	      UV(K)=0.0
             VW(K)=0.0
             UW(K)=0.0	      
	   Else   
	      Temp=1.0/Real(Num_par_K1(K))
	      UPm(K)  = UPm1(K) *Temp
	      VPm(K)  = VPm1(K) *Temp
	      WPm(K)  = WPm1(K) *Temp
	      
	      UPrms(K)  = UPrms1(K) *Temp
	      VPrms(K)  = VPrms1(K) *Temp
	      WPrms(K)  = WPrms1(K) *Temp
	 
	      UV(K)  = UPVP1(K) *Temp-UPm(K)*VPm(K)
	      VW(K)  = VPWP1(K) *Temp-VPm(K)*WPm(K)
	      UW(K)  = UPWP1(K) *Temp-UPm(K)*WPm(K)	 	      
	      
	      UPrms(K)  = SQRT(UPrms(K)-UPm(K)**2)
	      VPrms(K)  = SQRT(VPrms(K)-VPm(K)**2)
	      WPrms(K)  = SQRT(WPrms(K)-WPm(K)**2)
	           
	      
	   Endif
	    Num_par_K(K)= Num_par_K1(K)/REAL(FLAGPAR)
	    Num_par_K(K)= Num_par_K(K)/(ALX*ALZ)/DY(K)
	   
	   Cm(K)= Cm1(K)/REAL(FLAGPAR)/real(N1M*N3M)  
	   Enddo
	
	   Do K=1,N2M

	      UPm(K)  = UPm(K) *RE2
	      VPm(K)  = VPm(K) *RE2
	      WPm(K)  = WPm(K) *RE2
	      UV(K)  = UV(K) *RE2 *RE2
	      VW(K)  = VW(K) *RE2 *RE2
	      UW(K)  = UW(K) *RE2 *RE2		      
	      
	      UPrms(K)  = UPrms(K) *RE2
	      VPrms(K)  = VPrms(K) *RE2
	      WPrms(K)  = WPrms(K) *RE2	      
  
	       Cm(K)= Cm(K)   
	Enddo
	If(0.EQ.0)Then
	  Open (1,file='RESULT1/0_4Means_par.plt')
	  WRITE(1,*)
	1'VARIABLES="Y","UP","VP","WP","UPrms","VPrms","WPrms" ,"UV"'
	Write(1,*)'ZONE T="',IZONE,'" K=',M2-1
	DO K=1,M2-1
	Write(1,'(10E20.12)')YC(K)*HHH,UPm(K),VPm(K),
     > WPm(K),UPrms(K),VPrms(K),WPrms(K),UV(K)
	 
	  Enddo
	  Close(1)
	Endif
	If(0.EQ.0)Then
	  Open (1,file='RESULT1/0_4Means_par2.plt')
	  WRITE(1,*)
	1  'VARIABLES="Y","N","C","UP","VP","WP","UPrms","VPrms","WPrms" ,"UV"'
	  Write(1,*)'ZONE T="',IZONE,'" K=',(M2-1)/2
	  DO K=1,(M2-1)/2
	Write(1,'(10E20.12)')YC(K)*HHH,Num_par_K(M2-K)*QQ,Cm(M2-K),
     > UPm(M2-K),VPm(M2-K),WPm(M2-K),UPrms(M2-K),VPrms(M2-K),WPrms(M2-K)
     > ,UV(M2-K)
	  Enddo
	  Close(1)
	Endif
	ENDIF
100	FORMAT(A8,I6,A1)
200	FORMAT(A3,I6)
300	FORMAT(9F20.6)
      End
        
        
