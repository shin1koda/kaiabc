module model

  !$ use omp_lib
  implicit none


  real(8) :: Lam=0.0d0
  real(8),parameter :: DeltaT=0.01d0,TZERO=-273.15d0,AmpActEne=1000.0d0
  integer,dimension(1:8) :: EvalFlag=0
  integer,dimension(60)  :: RSTFlag=0
  real(8) :: TSt=0.0d0,TSp=90.0d0,Rat=1.0d0
  integer :: NShift=0
  integer :: Fuck=0,FEvalOut=0
  integer,parameter :: NList=25
  real(8),dimension(NList):: EvalDetails=0.0d0,WDetails=0.0d0
  real(8) :: BetaB=1.0d0

  integer,dimension(1:8,-1:6) :: HTable=0
  integer,dimension(0:7,0:1715)  :: IndTable=0

  real(8),dimension(1716,1716) :: UncoM2h=0.0d0,UncoM2ex=0.0d0,UncoM2Sp=0.0d0,&
                                  UncoM2Sm=0.0d0,UncoM2Tp=0.0d0,UncoM2Tm=0.0d0
  real(8),dimension(4800) :: M2h=0.0d0,M2ex=0.0d0
  integer,dimension(4800) :: Col2h=0,Col2ex=0
  real(8),dimension(2838) :: M2Sp=0.0d0,M2Sm=0.0d0,M2Tp=0.0d0,M2Tm=0.0d0
  integer,dimension(2838) :: Col2Sp=0,Col2Sm=0,Col2Tp=0,Col2Tm=0
  integer,dimension(1717) :: Row2h=0,Row2ex=0,Row2Sp=0,Row2Sm=0,Row2Tp=0,Row2Tm=0

  real(8),dimension(0:1715,0:1,0:1) :: PhosMask=0.0d0
  real(8),dimension(0:1715,0:1,0:1,0:1) :: NucPhosMask=0.0d0
  real(8),dimension(0:1715) :: C2ATPMask=0.0d0
  real(8),dimension(0:6) :: C1ATPMask=0.0d0


  real(8) :: lamhex=6.0d0

  type mx
    real(8) :: time=0.0d0
    real(8) :: amp=0.0d0
  end type mx

  type prm
    character(len=10) :: n=""
    integer :: f=0
    integer :: fen=0
    real(8) :: v=0.0d0
    real(8) :: iv=0.0d0
    real(8) :: v100=0.0d0
    real(8) :: en=0.0d0
    real(8) :: ien=0.0d0
    type(prm),pointer :: next=>null()
  end type prm

  type ext
    integer :: f=0
    integer :: i=0
    real(8) :: v=0.0d0
  end type ext

  type(prm),pointer :: First,Head,Prev,A0,B0,C0,&
                       k1h,k1eIR,k1eIF,k1eAR,k1eAF,k1IA0,k1IA6,k1AIB0,k1AIB6,&
                       k2h,k2eI,k2eA,&
                       k2SpR,k2SmR,k2TpR,k2TmR,k2SpTm,k2SmTp,&
                       k2SpF,k2SmF,k2TpF,k2TmF,&
                       K1A,K1gsB,K1fsB,&
                       K2RIA,K2FIA,K2A,&
                       K2RFTUU,K2RFTPU,K2RFTUP,K2RFTPP,&
                       K2RFDUU,K2RFDPU,K2RFDUP,K2RFDPP,&
                       kbpF,kbmF,kbpC,kbmC
  !$OMP threadprivate( First,Head,Prev,A0,B0,C0,&
  !$OMP                k1h,k1eIR,k1eIF,k1eAR,k1eAF,k1IA0,k1IA6,k1AIB0,k1AIB6,&
  !$OMP                k2h,k2eI,k2eA,&
  !$OMP                k2SpR,k2SmR,k2TpR,k2TmR,k2SpTm,k2SmTp,&
  !$OMP                k2SpF,k2SmF,k2TpF,k2TmF,&
  !$OMP                K1A,K1gsB,K1fsB,&
  !$OMP                K2RIA,K2FIA,K2A,&
  !$OMP                K2RFTUU,K2RFTPU,K2RFTUP,K2RFTPP,&
  !$OMP                K2RFDUU,K2RFDPU,K2RFDUP,K2RFDPP,&
  !$OMP                kbpF,kbmF,kbpC,kbmC)
  real(8),dimension(0:1) :: K1B=1.0d0
  !$OMP threadprivate(K1B)
  real(8),dimension(0:6) :: k1IA=1.0d0
  !$OMP threadprivate(k1IA)
  real(8),dimension(0:1715) :: K2RF=1.0d0,K2Atil=1.0d0,Gam2A0=0.0d0
  !$OMP threadprivate(K2RF,K2Atil,Gam2A0)
  real(8) :: ratio=1.0d0
  !$OMP threadprivate(ratio)
  integer,parameter :: NVM=60
  character(len=100),parameter :: fmt1="(A10,ES12.5e2,I5)"
  character(len=100),parameter :: fmt1t="(A10,ES12.5e2,I5,ES15.5e2,I5)"
  character(len=100),parameter :: fmt2="(A5,f7.2,f10.5)"

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine initStates(X,XMat,F)

    real(8),dimension(-1:1),intent(out) :: X
    real(8),dimension(1716,14),intent(out) :: XMat
    integer,intent(in),optional :: F
    integer :: I, I1
    integer,dimension(0:7) :: I2

    X=0.0d0
    X(0)=B0%v*(kbmF%v/(kbpF%v+kbmF%v))
    X(1)=B0%v-X(0)

    XMat=0.0d0
    if (present(F)) then
      I=XInd(0,(/0,0,0,0,2,0,2,2/))
    else
      I=XInd(0,(/0,0,0,0,6,0,0,0/))
    end if
    !I=XInd(0,(/0,0,0,0,0,6,0,0/))
    XMat(I+1,7)=k1AIB0%v/(k1AIB0%v+k1IA(6))
    XMat(I+1,14)=k1IA(6)/(k1AIB0%v+k1IA(6))

  end subroutine initStates
!!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine initPrmset()

    character(len=10) :: prmname

    !$OMP parallel private(prmname)
    allocate(First)
    Head => First

    prmname="A0"
    call appendPrm(A0    ,prmname,0,0.6d0)

    prmname="B0"
    call appendPrm(B0    ,prmname,0,3.5d0)

    prmname="C0"
    call appendPrm(C0    ,prmname,0,3.5d0)

    !C1
    prmname="k1hyd"
    call appendPrm(k1h   ,prmname,0,0.15d0)

    prmname="k1exInRig"
    call appendPrm(k1eIR ,prmname,0,0.000000d0)

    prmname="k1exInFle"
    call appendPrm(k1eIF ,prmname,1,10.0d0)

    prmname="k1IA0"
    call appendPrm(k1IA0 ,prmname,1,0.01d0)

    prmname="k1IA6"
    call appendPrm(k1IA6 ,prmname,1,10.0d0)

    prmname="k1AIB0"
    call appendPrm(k1AIB0 ,prmname,1,100.0d0)

    prmname="k1AIB6"
    call appendPrm(k1AIB6 ,prmname,1,100.0d0)

    prmname="k1exAcRig"
    call appendPrm(k1eAR ,prmname,0,0.000000d0)

    prmname="k1exAcFle"
    call appendPrm(k1eAF ,prmname,1,10.0d0)

    prmname="K1gsB"
    call appendPrm(K1gsB   ,prmname,1,0.1d0)

    prmname="K1fsB"
    call appendPrm(K1fsB   ,prmname,0,0.1d0)

    prmname="kb+Free"
    call appendPrm(kbpF ,prmname,1,0.01d0)

    prmname="kb-Free"
    call appendPrm(kbmF ,prmname,1,10.0d0)

    prmname="kb+Comp"
    call appendPrm(kbpC ,prmname,1,0.01d0)

    prmname="kb-Comp"
    call appendPrm(kbmC ,prmname,1,10.0d0)

    prmname="K1A"
    call appendPrm(K1A   ,prmname,1,0.00001d0)

    !C2
    prmname="k2hyd"
    call appendPrm(k2h   ,prmname,0,0.6d0)

    prmname="k2exIna"
    call appendPrm(k2eI ,prmname,0,0.000000d0)

    prmname="k2exAct"
    call appendPrm(k2eA ,prmname,1,40d0)

    prmname="k2S+R"
    call appendPrm(k2SpR  ,prmname,1,0.08d0)

    prmname="k2S-R"
    call appendPrm(k2SmR  ,prmname,1,0.08d0)

    prmname="k2S+F"
    call appendPrm(k2SpF  ,prmname,1,0.08d0)

    prmname="k2S-F"
    call appendPrm(k2SmF  ,prmname,1,0.08d0)

    prmname="k2T+R"
    call appendPrm(k2TpR  ,prmname,1,0.8d0)

    prmname="k2T-R"
    call appendPrm(k2TmR  ,prmname,1,0.8d0)

    prmname="k2T+F"
    call appendPrm(k2TpF  ,prmname,1,0.8d0)

    prmname="k2T-F"
    call appendPrm(k2TmF  ,prmname,1,0.8d0)

    prmname="k2S+T-"
    call appendPrm(k2SpTm,prmname,1,0.000012d0)

    prmname="k2S-T+"
    call appendPrm(k2SmTp,prmname,0,0.000012d0)

!Rigid or Flexible
    prmname="K2RFTUU"
    call appendPrm(K2RFTUU ,prmname,1,0.1d0)

    prmname="K2RFTUP"
    call appendPrm(K2RFTUP ,prmname,1,0.1d0)

    prmname="K2RFTPP"
    call appendPrm(K2RFTPP ,prmname,0,10d0)

    prmname="K2RFTPU"
    call appendPrm(K2RFTPU ,prmname,1,10d0)

    prmname="K2RFDUU"
    call appendPrm(K2RFDUU ,prmname,1,0.1d0)

    prmname="K2RFDUP"
    call appendPrm(K2RFDUP ,prmname,0,0.1d0)

    prmname="K2RFDPP"
    call appendPrm(K2RFDPP ,prmname,0,10d0)

    prmname="K2RFDPU"
    call appendPrm(K2RFDPU ,prmname,0,10d0)

!Buried or Exposed
    prmname="K2RIA"
    call appendPrm(K2RIA ,prmname,1,100000d0)

    prmname="K2FIA"
    call appendPrm(K2FIA ,prmname,1,100d0)

    prmname="K2A"
    call appendPrm(K2A   ,prmname,1,0.005d0)

    !$OMP end parallel
  contains

    subroutine appendPrm(Alias,Nam,Flag,Val,IniVal,FlagEne,Energy,IniEne)

      type(prm),pointer :: Alias,Next
      character(len=10) :: Nam
      integer,optional :: Flag,FlagEne
      real(8),optional :: Val,IniVal,Energy,IniEne


      Alias => Head
      allocate(Next)
      Head => Next
      Alias%n  = Nam
      if (present(Flag)) Alias%f  = Flag
      if (present(Val))    Alias%v  = Val
      if (present(IniVal)) then
        Alias%iv = IniVal
      else
        Alias%iv = Alias%v
      end if
      if (present(FlagEne)) Alias%fen  = FlagEne
      if (present(Energy))    Alias%en  = Energy
      if (present(IniEne)) then
        Alias%ien = IniEne
      else
        Alias%ien = Alias%en
      end if
      Alias%next => Head

    end subroutine appendPrm

  end subroutine initPrmset
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine reloadPrmset(E,TempFlag,Te)

    real(8),dimension(1:),intent(in) :: E
    type(prm),pointer :: P
    integer,optional,intent(in) :: TempFlag
    integer :: Flag
    real(8),optional,intent(in) :: Te
    real(8) :: G
    integer :: I,NVar

    if (present(TempFlag)) then
      Flag=TempFlag
    else
      Flag=0
    endif

    I=0
    P=>First
    do while(associated(P%next))

      G=0.0d0

      if (P%f/=0) then
        I=I+1
        G=G+E(I)
      end if

      if ((Flag/=0).and.(P%fen/=0)) then
        I=I+1
        P%en=P%ien-AmpActEne*E(I)
      end if

      if (Flag/=0) then
        G=G+P%en*(1.0d0/(Te-TZERO)-1.0d0/(30.0d0-TZERO))
      end if

      P%v=P%iv*exp(-G)
      P=>P%next
    end do

    P=>First
    do while(associated(P%next))
      if ((P%f==0).and.((Flag==0).or.(P%fen==0))) then
        P%v=getVal(P,Flag,Te)
      end if
      P=>P%next
    end do

    contains

      function getVal(Pnt,TFlag,Te)

        real(8) :: getVal
        type(prm),pointer :: Pnt
        integer,intent(in) :: TFlag
        real(8),optional,intent(in) :: Te

        select case (trim(Pnt%n))
          case("k1exAcFle")
            getVal = k1eIF%v*Lam
          case("k1AIB6")
            getVal = k1AIB0%v
          !case("K1fsB")
          !  getVal = kbpF%v / kbmF%v * kbmC%v /kbpC%v * K1gsB%v
          case("k2S-T+")
            getVal = k2SpTm%v * k2SmR%v * k2TpR%v /k2SpR%v /k2TmR%v
          !case("K2FIA")
          !  getVal = K2RIA%v*1.0d0
          case("K2RFTUP")
            getVal = K2RFTUU%v/K2RFDUU%v * K2RFDUP%v
          case("K2RFDPP")
            getVal = K2RFDUU%v * K2RFDUP%v/K2RFDUU%v * K2RFDPU%v/K2RFDUU%v
          case("K2RFTPP")
            getVal = K2RFTUU%v * K2RFDUP%v/K2RFDUU%v * K2RFDPU%v/K2RFDUU%v
          case("K2RFTPU")
            getVal = K2RFTUU%v/K2RFDUU%v * K2RFDPU%v
          case("k2S+R")
            getVal = (K2RFDPU%v/K2RFTUU%v)**(1.0d0/6.0d0) *k2SpF%v/k2SmF%v * k2SmR%v
          case("k2T+R")
            getVal = (K2RFDUP%v/K2RFTUU%v)**(1.0d0/6.0d0) *k2TpF%v/k2TmF%v * k2TmR%v
          case("k2S-R")
            getVal = (K2RFTUU%v/K2RFDPU%v)**(1.0d0/6.0d0) *k2SmF%v/k2SpF%v * k2SpR%v
          case("k2T-R")
            getVal = (K2RFTUU%v/K2RFDUP%v)**(1.0d0/6.0d0) *k2TmF%v/k2TpF%v * k2TpR%v
          case default
            getVal=Pnt%iv
            if (TFlag>0) getVal=getVal*exp(-Pnt%ien*(1.0d0/(Te-TZERO)-1.0d0/(30.0d0-TZERO)))
        end select

      end function getVal

  end subroutine reloadPrmset

  subroutine reloadPrmsetTemp(E,Te)
    real(8),dimension(1:),intent(in) :: E
    type(prm),pointer :: P
    real(8),optional,intent(in) :: Te

    call reloadPrmset(E,1,100.0d0)

    P=>First
    do while(associated(P%next))
      P%v100=P%v
      P=>P%next
    end do

    call reloadPrmset(E,1,Te)

    P=>First
    do while(associated(P%next))
      if ((P%f==0).and.(P%fen==0)) then
        if (P%v/=0.0d0) then
          P%en = log(P%v100/P%v)/(1.0d0/(Te-TZERO)-1.0d0/(100.0d0-TZERO))
        else
          P%en = 0.0d0
        end if
      end if
      P=>P%next
    end do
  end subroutine reloadPrmsetTemp
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine printPrmset(Filename,TempFlag)

    character(len=100),intent(in) :: Filename
    integer,optional,intent(in) :: TempFlag
    integer :: Flag
    type(prm),pointer :: P

    if (present(TempFlag)) then
      Flag=TempFlag
    else
      Flag=0
    endif

    open(10,file=trim(Filename),status="unknown")

    P=>First
    do while(associated(P%next))
      if (Flag>0) then
        write(10,fmt1t) P%n,P%v,P%f,P%en,P%fen
        write( *,fmt1t) P%n,P%v,P%f,P%en,P%fen
      else
        write(10,fmt1) P%n,P%v,P%f
        write( *,fmt1) P%n,P%v,P%f
      end if
      P=>P%next
    end do

    close(10)


  end subroutine printPrmset
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine readPrmset(Filename,TempFlag)

    character(len=100),intent(in) :: Filename
    integer,optional,intent(in) :: TempFlag
    integer :: Flag,I
    type(prm),pointer :: P

    if (present(TempFlag)) then
      Flag=TempFlag
    else
      Flag=0
    endif

    !$OMP parallel private(I,P)
    I=10
    !$ I=I+omp_get_thread_num()
    open(I,file=trim(Filename),status="old")

    P=>First
    do while(associated(P%next))
      if (Flag>0) then
        read(I,fmt1t) P%n,P%v,P%f,P%en,P%fen
        P%iv=P%v
        P%ien=P%en
        !$OMP master
        write(*,fmt1t) P%n,P%v,P%f,P%en,P%fen
        !$OMP end master
      else
        read(I,fmt1) P%n,P%v,P%f
        P%iv=P%v
        !$OMP master
        write(*,fmt1) P%n,P%v,P%f
        !$OMP end master
      end if
      P=>P%next
    end do

    close(I)
    !$OMP end parallel

  end subroutine readPrmset
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine mergePrmsets(Filename1,Te1,Filename2,Te2,Filename)

    character(len=100),intent(in) :: Filename1,Filename2,Filename
    real(8),optional,intent(in) :: Te1,Te2
    type(prm),pointer :: P
    real(8) :: Val1,Val2
    integer :: F1,F2


    open(10,file=trim(Filename),status="unknown")
    open(11,file=trim(Filename1),status="old")
    open(12,file=trim(Filename2),status="old")

    P=>First
    do while(associated(P%next))
      read(11,fmt1) P%n,Val1,F1
      read(12,fmt1) P%n,Val2,F2
      P%f=max(F1,F2)
      if ((Val1==0.0d0).and.(Val2==0.0d0)) then
        P%en=0.0d0
      else
        P%en = -log(Val1/Val2)/(1.0d0/(Te1-TZERO)-1.0d0/(Te2-TZERO))
      end if
      P%v = Val2 * exp(-P%en* (1.0d0/(30.0d0-TZERO)-1.0d0/(Te2-TZERO)) ) 
      write(10,fmt1t) P%n,P%v,P%f,P%en,P%f
      write( *,fmt1t) P%n,P%v,P%f,P%en,P%f
      P=>P%next
    end do

    close(10)
    close(11)
    close(12)
    
  end subroutine mergePrmsets
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine getRSTFlag(RFlag,TempFlag)

    integer,dimension(60),intent(inout) :: RFlag
    integer,optional,intent(in) :: TempFlag
    integer :: Flag,I
    type(prm),pointer :: P
    character(len=10) :: CName

    RFlag=0

    if (present(TempFlag)) then
      Flag=TempFlag
    else
      Flag=0
    endif

    I=0
    P=>First
    do while(associated(P%next))
      CName = P%n
      if (P%f/=0) then
        I=I+1
        if ((Flag==0).and.(CName(1:1)=="k")) RFlag(I)=1
      end if
      if ((Flag>0).and.(P%fen/=0)) then
        I=I+1
        if (CName(1:1)=="k") RFlag(I)=1
      end if
      P=>P%next
    end do



  end subroutine getRSTFlag
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  function getNVar(TempFlag)

    integer :: getNVar
    integer,optional,intent(in) :: TempFlag
    integer :: Flag
    type(prm),pointer :: P

    if (present(TempFlag)) then
      Flag=TempFlag
    else
      Flag=0
    endif

    getNVar=0
    P=>First
    do while(associated(P%next))
      if (P%f/=0) getNVar=getNVar+1
      if ((Flag>0).and.(P%fen/=0)) getNVar=getNVar+1
      P=>P%next
    end do

  end function getNVar
!-------------------------------------------------------------------------------

  subroutine getDataPar(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs)

    integer,intent(in) :: NTime,Span
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos
    real(8),dimension(0:1,0:1,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMatSta=0.0d0
    real(8) :: Af
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0

    call reloadPrmsetTemp(E,30.0d0)
    call initEqConst()
    !call getStationary(XMatSta)


    call setOMPMKL(4,4)
    !$OMP parallel private(FA,FB,ITime,NT,X,XMat,Af,Bf)
    !$ call reloadPrmsetTemp(E,30.0d0)
    !$ call initEqConst()
    !$ FB=mod(omp_get_thread_num(),2)
    !$ FA=mod(omp_get_thread_num()/2,2)
    !$ A0%v = dble(FA)*A0%iv
    !$ B0%v = dble(FB)*B0%iv
    !$ call initStates(X,XMat)
    !$ !XMat=XMatSta
    !$OMP barrier
    !$ NT=NTime
    !$ if (omp_get_thread_num()==3) NT=NTime
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT*0.5d0)
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT*0.5d0)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     Phos(FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
    !$     ATP(FA,FB,1:2,ITime/Span)=getATP(XMat)
    !$     call getATPase(X,XMat,ATPase(FA,FB,1:6,ITime/Span))
    !$     call getAfBf(X,XMat,Afs(FA,FB,1,ITime/Span),Bfs(FA,FB,0:1,ITime/Span))
    !$     !if (omp_get_thread_num()==3) then
    !$     !  call getAfBf(X,XMat,Af,Bf)
    !$     !  write(*,*) dble(ITime)*DeltaT,sum(Bf),sum(Bf(0:1)/K1B(0:1))
    !$     !end if
    !$OMP end critical
    !$   end if
    !$ end do
    !$OMP end parallel


  end subroutine getDataPar

  subroutine getDataPar_TOut(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs,Conf)

    integer,intent(in) :: NTime,Span
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos,Conf
    real(8),dimension(0:7,0:1,0:1,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:7,0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Conf=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(8)
    do FA=1,1
    do FB=1,1
      if (FA*2+FB==0) then
        write(*,*) "getting data C"
      else if (FA*2+FB==1) then
        write(*,*) "getting data CB"
      else if (FA*2+FB==2) then
        write(*,*) "getting data CA"
      else if (FA*2+FB==3) then
        write(*,*) "getting data CAB"
      end if
      if (FA*2+FB==3) then
        NT=NTime
      else
        NT=4800
      end if
      !$OMP parallel private(ITime,X,XMat,Af,Bf,ITe,Te)
      !$ ITe=omp_get_thread_num()
      !$ Te=dble(ITe)*2.0d0+26.0d0
      !$ call reloadPrmsetTemp(E,Te)
      !$ call initEqConst()
      !$ A0%v = dble(FA)*A0%iv
      !$ B0%v = dble(FB)*B0%iv
      !$ call initStates(X,XMat)
      !$ if (FA*2+FB==3) XMat=XMat0C
      !$ !XMat=XMat0C
      !$ !call resetNuc(XMat)
      !$OMP barrier
      !$ do ITime=0,NT
      !$   if (Fuck==1) exit
      !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
      !$   if (mod(ITime,Span)==0) then
      !$OMP critical
      !$     Phos(ITe,FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
      !$     ATP(ITe,FA,FB,1:2,ITime/Span)=getATP(XMat)
      !$     call getATPase(X,XMat,ATPase(ITe,FA,FB,1:6,ITime/Span))
      !$     call getAfBf(X,XMat,Afs(ITe,FA,FB,1,ITime/Span),&
      !$       Bfs(Ite,FA,FB,0:1,ITime/Span))
      !$     Conf(ITe,FA,FB,0:1,0:1,ITime/Span)=getConf(X,XMat)
      !$OMP end critical
      !$   end if
      !$ end do
      !$OMP end parallel
    end do
    end do


  end subroutine getDataPar_TOut


  subroutine getDataPar_Out30c(NTime,Span,E,Phos,ATP,ATPase,ADPRel,Afs,Bfs,Conf,C1ATP,ITemp,IState)

    integer,intent(in) :: NTime,Span
    integer,intent(in) :: ITemp
    integer,optional,intent(in) :: IState
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos,Conf
    real(8),dimension(0:7,0:1,0:1,2,0:1200),intent(out) :: ATP,ADPRel
    real(8),dimension(0:7,0:1,0:1,0:6,0:1200),intent(out) :: C1ATP
    real(8),dimension(0:7,0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    ADPRel=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Conf=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(4)
    if (FA*2+FB==0) then
      write(*,*) "getting data C"
    else if (FA*2+FB==1) then
      write(*,*) "getting data CB"
    else if (FA*2+FB==2) then
      write(*,*) "getting data CA"
    else if (FA*2+FB==3) then
      write(*,*) "getting data CAB"
    end if
    ITe=2
    !$OMP parallel private(ITime,X,XMat,Af,Bf,FA,FB,NT)
    !$ FB=mod(omp_get_thread_num(),2)
    !$ FA=omp_get_thread_num()/2
    !$ if (FA*2+FB==3) then
    !$   NT=NTime
    !$   if (present(IState)) NT=4800
    !$ else
    !$   NT=4800
    !$ end if
    !$ call reloadPrmsetTemp(E,dble(ITemp))
    !$ call initEqConst()
    !$ A0%v = dble(FA)*A0%iv
    !$ B0%v = dble(FB)*B0%iv
    !$ call initStates(X,XMat)
    !$ if (present(IState)) XMat=XMat0C
    !$OMP barrier
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     Phos(ITe,FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
    !$     ATP(ITe,FA,FB,1:2,ITime/Span)=getATP(XMat)
    !$     C1ATP(ITe,FA,FB,0:6,ITime/Span)=getC1ATP(XMat)
    !$     call getATPase(X,XMat,ATPase(ITe,FA,FB,1:6,ITime/Span))
    !$     call getADPRelease(X,XMat,ADPRel(ITe,FA,FB,1:2,ITime/Span))
    !$     call getAfBf(X,XMat,Afs(ITe,FA,FB,1,ITime/Span),&
    !$       Bfs(Ite,FA,FB,0:1,ITime/Span))
    !$     Conf(ITe,FA,FB,0:1,0:1,ITime/Span)=getConf(X,XMat)
    !$OMP end critical
    !$   end if
    !$ end do
    !$OMP end parallel


  end subroutine getDataPar_Out30c


  subroutine getDataPar_Out30c_Complex(NTime,Span,E,Comp,ITemp)

    integer,intent(in) :: NTime,Span,ITemp
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:1,0:1,0:6,0:6,0:1,0:1200),intent(out) :: Comp

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Comp=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(4)
    if (FA*2+FB==0) then
      write(*,*) "getting data C"
    else if (FA*2+FB==1) then
      write(*,*) "getting data CB"
    else if (FA*2+FB==2) then
      write(*,*) "getting data CA"
    else if (FA*2+FB==3) then
      write(*,*) "getting data CAB"
    end if
    ITe=2
    !$OMP parallel private(ITime,X,XMat,Af,Bf,FA,FB,NT)
    !$ FB=mod(omp_get_thread_num(),2)
    !$ FA=omp_get_thread_num()/2
    !$ if (FA*2+FB==3) then
    !$   NT=NTime
    !$ else
    !$   NT=4800
    !$ end if
    !$ call reloadPrmsetTemp(E,dble(ITemp))
    !$ call initEqConst()
    !$ A0%v = dble(FA)*A0%iv
    !$ B0%v = dble(FB)*B0%iv
    !$ call initStates(X,XMat)
    !$ !if (FA*2+FB==3) XMat=XMat0C
    !$ !XMat=XMat0C
    !$ !call resetNuc(XMat)
    !$OMP barrier
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     call getComplex(X,XMat,Comp(FA,FB,0:6,0:6,0:1,ITime/Span))
    !$OMP end critical
    !$   end if
    !$ end do
    !$OMP end parallel


  end subroutine getDataPar_Out30c_Complex



  subroutine getDataPar_TOutReset(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs,Conf)

    integer,intent(in) :: NTime,Span
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos,Conf
    real(8),dimension(0:7,0:1,0:1,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:7,0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Conf=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(4)
    do FA=0,0
    do FB=0,0
      if (FA*2+FB==0) then
        write(*,*) "getting data C"
      else if (FA*2+FB==1) then
        write(*,*) "getting data CB"
      else if (FA*2+FB==2) then
        write(*,*) "getting data CA"
      else if (FA*2+FB==3) then
        write(*,*) "getting data CAB"
      end if
      if (FA*2+FB==3) then
        NT=NTime
      else
        NT=4800
      end if
      !$OMP parallel private(ITime,X,XMat,Af,Bf,ITe,Te)
      !$ ITe=omp_get_thread_num()
      !$ Te=dble(ITe)*5.0d0+25.0d0
      !$ call reloadPrmsetTemp(E,Te)
      !$ call initEqConst()
      !$ A0%v = dble(FA)*A0%iv
      !$ B0%v = dble(FB)*B0%iv
      !$ call initStates(X,XMat)
      !$ XMat=XMat0C
      !$ call resetNuc(XMat)
      !$OMP barrier
      !$ do ITime=0,NT
      !$   if (Fuck==1) exit
      !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
      !$   if (mod(ITime,Span)==0) then
      !$OMP critical
      !$     Phos(ITe,FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
      !$     ATP(ITe,FA,FB,1:2,ITime/Span)=getATP(XMat)
      !$     call getATPase(X,XMat,ATPase(ITe,FA,FB,1:6,ITime/Span))
      !$     call getAfBf(X,XMat,Afs(ITe,FA,FB,1,ITime/Span),&
      !$       Bfs(Ite,FA,FB,0:1,ITime/Span))
      !$     Conf(ITe,FA,FB,0:1,0:1,ITime/Span)=getConf(X,XMat)
      !$OMP end critical
      !$   end if
      !$ end do
      !$OMP end parallel
    end do
    end do


  end subroutine getDataPar_TOutReset

  subroutine getDataPar_Outk1h(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs,Conf,C1ATP,Ik1h,Flag)

    integer,intent(in) :: NTime,Span,Ik1h,Flag
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos,Conf
    real(8),dimension(0:7,0:1,0:1,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:7,0:1,0:1,0:6,0:1200),intent(out) :: C1ATP
    real(8),dimension(0:7,0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Conf=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(8)
    FA=0
    FB=0
    !$OMP parallel private(ITime,X,XMat,Af,Bf,NT,ITe)
    !$ NT=NTime
    !$ call reloadPrmsetTemp(E,30.0d0)
    !$ call initEqConst()
    !$ ITe=omp_get_thread_num()
    !$ if (Flag==1) then
    !$   k1h%v = dble(ITe+Ik1h)*0.025d0
    !$ else if (Flag==2) then
    !$   k2h%v = dble(ITe+Ik1h)*0.025d0
    !$ end if
    !$ call initStates(X,XMat)
    !$ !if (FA*2+FB==3) XMat=XMat0C
    !$ !XMat=XMat0C
    !$ !call resetNuc(XMat)
    !$OMP barrier
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     Phos(ITe,FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
    !$     ATP(ITe,FA,FB,1:2,ITime/Span)=getATP(XMat)
    !$     C1ATP(ITe,FA,FB,0:6,ITime/Span)=getC1ATP(XMat)
    !$     call getATPase(X,XMat,ATPase(ITe,FA,FB,1:6,ITime/Span))
    !$     call getAfBf(X,XMat,Afs(ITe,FA,FB,1,ITime/Span),&
    !$       Bfs(Ite,FA,FB,0:1,ITime/Span))
    !$     Conf(ITe,FA,FB,0:1,0:1,ITime/Span)=getConf(X,XMat)
    !$OMP end critical
    !$   end if
    !$ end do
    !$OMP end parallel


  end subroutine getDataPar_Outk1h

  subroutine getDataPar_OutATot(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs,Conf,C1ATP,IATot)

    integer,intent(in) :: NTime,Span,IATot
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos,Conf
    real(8),dimension(0:7,0:1,0:1,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:7,0:1,0:1,0:6,0:1200),intent(out) :: C1ATP
    real(8),dimension(0:7,0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Conf=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(8)
    FA=0
    FB=0
    !$OMP parallel private(ITime,X,XMat,Af,Bf,NT,ITe)
    !$ NT=NTime
    !$ call reloadPrmsetTemp(E,30.0d0)
    !$ call initEqConst()
    !$ ITe=omp_get_thread_num()
    !$ A0%v = dble(ITe+IATot)*0.025d0
    !$ call initStates(X,XMat)
    !$ !if (FA*2+FB==3) XMat=XMat0C
    !$ !XMat=XMat0C
    !$ !call resetNuc(XMat)
    !$OMP barrier
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     Phos(ITe,FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
    !$     ATP(ITe,FA,FB,1:2,ITime/Span)=getATP(XMat)
    !$     C1ATP(ITe,FA,FB,0:6,ITime/Span)=getC1ATP(XMat)
    !$     call getATPase(X,XMat,ATPase(ITe,FA,FB,1:6,ITime/Span))
    !$     call getAfBf(X,XMat,Afs(ITe,FA,FB,1,ITime/Span),&
    !$       Bfs(Ite,FA,FB,0:1,ITime/Span))
    !$     Conf(ITe,FA,FB,0:1,0:1,ITime/Span)=getConf(X,XMat)
    !$OMP end critical
    !$   end if
    !$ end do
    !$OMP end parallel


  end subroutine getDataPar_OutATot

  subroutine getDataPar_OutBtot(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs,Conf,C1ATP,IBtot)

    integer,intent(in) :: NTime,Span,IBtot
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200),intent(out) :: Phos,Conf
    real(8),dimension(0:7,0:1,0:1,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:7,0:1,0:1,0:6,0:1200),intent(out) :: C1ATP
    real(8),dimension(0:7,0:1,0:1,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,0:1,0:1,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1,0:1,0:1200),intent(out) :: Bfs

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af,Te
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,ITe

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Conf=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call setOMPMKL(8)
    FA=0
    FB=0
    !$OMP parallel private(ITime,X,XMat,Af,Bf,NT,ITe)
    !$ NT=NTime
    !$ call reloadPrmsetTemp(E,30.0d0)
    !$ call initEqConst()
    !$ ITe=omp_get_thread_num()
    !$ B0%v = dble(ITe+IBtot)*0.1d0
    !$ call initStates(X,XMat)
    !$ !if (FA*2+FB==3) XMat=XMat0C
    !$ !XMat=XMat0C
    !$ !call resetNuc(XMat)
    !$OMP barrier
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     Phos(ITe,FA,FB,0:1,0:1,ITime/Span)=getPhos(XMat)
    !$     ATP(ITe,FA,FB,1:2,ITime/Span)=getATP(XMat)
    !$     C1ATP(ITe,FA,FB,0:6,ITime/Span)=getC1ATP(XMat)
    !$     call getATPase(X,XMat,ATPase(ITe,FA,FB,1:6,ITime/Span))
    !$     call getAfBf(X,XMat,Afs(ITe,FA,FB,1,ITime/Span),&
    !$       Bfs(Ite,FA,FB,0:1,ITime/Span))
    !$     Conf(ITe,FA,FB,0:1,0:1,ITime/Span)=getConf(X,XMat)
    !$OMP end critical
    !$   end if
    !$ end do
    !$OMP end parallel


  end subroutine getDataPar_OutBtot

  function getPeriod_k1h(k1,E)

    real(8) :: getPeriod_k1h
    real(8),intent(in) :: k1
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0
    real(8),dimension(0:1,0:1) :: Phos=0.0d0
    real(8),dimension(0:2) :: U=0.0d0
    integer ::I,Cnt
    integer,dimension(0:1)::IMax

    call reloadPrmsetTemp(E,30.0d0)
    call initEqConst()
    k1h%v=k1
    call initStates(X,XMat)
    Phos=getPhos(XMat)
    U(0)=Phos(0,0)
    Cnt=0
    do I=1,10000
      call R_Runge4thCustomized(X,XMat,DeltaT)
      Phos=getPhos(XMat)
      U(mod(I,3))=Phos(0,0)
      if ((I>1000).and.(U(mod(I,3))<U(mod(I-1,3))).and.(U(mod(I-1,3))>U(mod(I-2,3)))) then
        IMax(cnt)=I-1
        cnt = cnt+1
      end if
      if (cnt==2) then
        getPeriod_k1h = dble(IMax(1)-IMax(0))*DeltaT
      end if
    end do

  end function getPeriod_k1h


  subroutine getDataPar_T(NTime,Span,E,Phos,ATP,ATPase,Afs,Bfs,Q10,Phos0,Peri04)

    integer,intent(in) :: NTime,Span
    real(8),dimension(1:),intent(in) :: E
    real(8),dimension(0:7,0:1,0:1,0:1200),intent(out) :: Phos
    real(8),dimension(0:7,2,0:1200),intent(out) :: ATP
    real(8),dimension(0:7,6,0:1200),intent(out) :: ATPase
    real(8),dimension(0:7,1,0:1200),intent(out) :: Afs
    real(8),dimension(0:7,0:1,0:1200),intent(out) :: Bfs
    real(8),dimension(2,0:1),intent(out) :: Q10
    real(8),dimension(0:1,0:1),intent(out) :: Phos0
    real(8),intent(out) :: Peri04
    real(8),dimension(6) :: TempATPase

    real(8),dimension(-1:1) :: X=0.0d0
    real(8),dimension(1716,14) :: XMat=0.0d0,XMat0c=0.0d0
    real(8) :: Af
    real(8),dimension(0:1) :: Bf
    integer :: FA,FB,ITime,NT,I

    Phos=0.0d0
    ATP=0.0d0
    ATPase=0.0d0
    Afs=0.0d0
    Bfs=0.0d0
    Q10=0.0d0
    Peri04=0.0d0

    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)

    call reloadPrmsetTemp(E,40.0d0)
    call initEqConst()
    A0%v=0.0d0
    B0%v=0.0d0
    call initStates(X,XMat)
    XMat=XMat0c
    call resetNuc(XMat)
    do I=1,10
      call R_Runge4thCustomized(X,XMat,DeltaT)
    end do
    call getATPase(X,XMat,TempATPase)
    Q10(1,0)=sum(TempATPase)
    Q10(2,0)=Q10(1,0)
    call getStationary(XMat)
    call getATPase(X,XMat,TempATPase)
    Q10(1,1)=sum(TempATPase)
    Q10(2,1)=Q10(1,1)

    call reloadPrmsetTemp(E,35.0d0)
    call initEqConst()
    A0%v=0.0d0
    B0%v=0.0d0
    call initStates(X,XMat)
    XMat=XMat0c
    call resetNuc(XMat)
    do I=1,10
      call R_Runge4thCustomized(X,XMat,DeltaT)
    end do
    call getATPase(X,XMat,TempATPase)
    Q10(2,0)=Q10(2,0)/sum(TempATPase)
    call getStationary(XMat)
    call getATPase(X,XMat,TempATPase)
    Q10(2,1)=Q10(2,1)/sum(TempATPase)

    call reloadPrmsetTemp(E,25.0d0)
    call initEqConst()
    A0%v=0.0d0
    B0%v=0.0d0
    call initStates(X,XMat)
    XMat=XMat0c
    call resetNuc(XMat)
    do I=1,10
      call R_Runge4thCustomized(X,XMat,DeltaT)
    end do
    call getATPase(X,XMat,TempATPase)
    Q10(1,0)=Q10(1,0)/sum(TempATPase)
    call getStationary(XMat)
    call getATPase(X,XMat,TempATPase)
    Q10(1,1)=Q10(1,1)/sum(TempATPase)


    call reloadPrmsetTemp(E,0.0d0)
    call initEqConst()
    call getStationary(XMat0c)
    Phos0 = getPhos(XMat0c)

    call setOMPMKL(8)
    !$OMP parallel private(FA,FB,ITime,NT,X,XMat,Af,Bf)
    !$ if (omp_get_thread_num()<4) then
    !$   call reloadPrmsetTemp(E,30.0d0)
    !$   call initEqConst()
    !$   FA=mod(omp_get_thread_num()/2,2)
    !$   FB=mod(omp_get_thread_num(),2)
    !$   A0%v = dble(FA)*A0%iv
    !$   B0%v = dble(FB)*B0%iv
    !$   if ((FA==0).and.(FB==1)) A0%v = 3.0d0*A0%iv
    !$ else if (omp_get_thread_num()==4) then
    !$   call reloadPrmsetTemp(E,38.0d0)
    !$   call initEqConst()
    !$ else if (omp_get_thread_num()==5) then
    !$   call reloadPrmsetTemp(E,26.0d0)
    !$   call initEqConst()
    !$ else if (omp_get_thread_num()==6) then
    !$   call reloadPrmsetTemp(E,28.0d0)
    !$   call initEqConst()
    !$ else if (omp_get_thread_num()==7) then
    !$   call reloadPrmsetTemp(E,40.0d0)
    !$   call initEqConst()
    !$ end if
    !$ call initStates(X,XMat)
    !$ if (omp_get_thread_num()>2) XMat=XMat0C
    !$ !if (omp_get_thread_num()>2) call initStates(X,XMat,1)
    !$OMP barrier
    !$ if (omp_get_thread_num()<8) then
    !$ NT=NTime
    !$ do ITime=0,NT
    !$   if (Fuck==1) exit
    !$   if ((omp_get_thread_num()==0).and.(ITime>4800)) exit
    !$   if (ITime/=0) call R_Runge4thCustomized(X,XMat,DeltaT)
    !$   if (mod(ITime,Span)==0) then
    !$OMP critical
    !$     Phos(omp_get_thread_num(),0:1,0:1,ITime/Span)=getPhos(XMat)
    !$     ATP(omp_get_thread_num(),1:2,ITime/Span)=getATP(XMat)
    !$     call getATPase(X,XMat,ATPase(omp_get_thread_num(),1:6,ITime/Span))
    !$     call getAfBf(X,XMat,Afs(omp_get_thread_num(),1,ITime/Span),&
    !$       Bfs(omp_get_thread_num(),0:1,ITime/Span))
    !$     if (omp_get_thread_num()==3) then
    !$       call getAfBf(X,XMat,Af,Bf)
    !$       !write(*,*) dble(ITime)*DeltaT,sum(Bf),sum(Bf(0:1)/K1B(0:1))
    !$       !write(*,*) dble(ITime)*DeltaT,Bf(0),Bf(1)
    !$       !write(*,*) dble(ITime)*DeltaT, ATPase(3,1,Itime/Span)/ATP(3,1,ITime/Span)/24.0d0
    !$     end if
    !$OMP end critical
    !$   end if
    !$ end do
    !$   if (omp_get_thread_num()==0) Peri04=getPeriod_k1h(0.4d0,E)
    !$ end if
    !$OMP end parallel


  end subroutine getDataPar_T

  function EvalFuncPar_T(E)

    real(8) :: EvalFuncPar_T
    real(8),dimension(NList) :: EvalList,WList
    real(8),dimension(1:),intent(in) :: E
    type dataexp
      real(8),dimension(0:100) :: t=0.0d0
      real(8),dimension(0:100) :: d=0.0d0
      real(8),dimension(0:100) :: std=0.0d0
      real(8),dimension(0:100) :: w=0.0d0
    end type
    type(dataexp),dimension(0:1,0:1,1:8),save :: PhosExp
    real(8),dimension(0:7,0:1,0:1,0:1200) :: Phos
    real(8),dimension(0:7,2,0:1200) :: ATP
    real(8),dimension(0:7,6,0:1200) :: ATPase
    real(8),dimension(0:7,1,0:1200) :: Afs
    real(8),dimension(0:7,0:1,0:1200) :: Bfs
    real(8),dimension(2,0:1) :: Q10
    real(8),dimension(0:1,0:1) :: Phos0
    real(8),dimension(0:7,0:1200) :: TotalATPase
    real(8),dimension(0:1200) :: TempData
    type(ext),dimension(10) :: Extremum
    real(8) :: TotalW,Temp,DataMax,DataMin,Al,Period,Peri04
    real(8),dimension(5) :: Periods
    integer :: I,J,ITe,ITime,NT
    character(len=2) :: CTe
    character(len=10) :: CName
    type(prm),pointer :: P
    integer,save :: Flag=0


    if (Flag==0) then
      write(*,*) "Reading experimental data"
      do ITe=1,8
        write(CTe,'(I2.2)') ITe*2+24
        open(10,file="data/exp/tempcomp/phos"//CTe//"c.dat",status="old")
        do I=0,40
          read(10,"(F5.0,8(F9.6))") PhosExp(0,0,Ite)%t(I), &
            PhosExp(1,1,Ite)%d(I),   PhosExp(0,1,Ite)%d(I),&
            PhosExp(1,0,Ite)%d(I),   PhosExp(0,0,Ite)%d(I), &
            PhosExp(1,1,Ite)%std(I), PhosExp(0,1,Ite)%std(I),&
            PhosExp(1,0,Ite)%std(I), PhosExp(0,0,Ite)%std(I)
          PhosExp(0,1,ITe)%t(I)=PhosExp(0,0,Ite)%t(I)
          PhosExp(1,0,ITe)%t(I)=PhosExp(0,0,Ite)%t(I)
          PhosExp(1,1,ITe)%t(I)=PhosExp(0,0,Ite)%t(I)
        end do
        close(10)

        !ST
        PhosExp(0,0,ITe)%w(0:40)=5.0d0
        !SpT
        PhosExp(0,1,ITe)%w(0:40)=3.0d0
        !pSpT
        PhosExp(1,1,ITe)%w(0:40)=1.0d0
        !pST
        PhosExp(1,0,ITe)%w(0:40)=3.0d0

        do I=0,1
        do J=0,1
          !PhosExp(I,J,ITe)%w( 0: 0)=PhosExp(I,J,ITe)%w( 0: 0)*25.0d0
          PhosExp(I,J,ITe)%w( 0: 5)=PhosExp(I,J,ITe)%w( 0: 5)*3.0d0
          !PhosExp(I,J,ITe)%w( 0: 1)=PhosExp(I,J,ITe)%w( 0: 1)*0.0d0
          !PhosExp(I,J,ITe)%w(26:30)=PhosExp(I,J,ITe)%w(26:30)*10.0d0
          !PhosExp(I,J,ITe)%w(21:40)=PhosExp(I,J,ITe)%w(21:40)*0.0d0
          !if ((ITe==1).or.(Ite==3)) PhosExp(I,J,ITe)%w(0:30)=PhosExp(I,J,ITe)%w(0:30)*2.0d0
        end do
        end do
      end do

      Flag=1
    end if

    EvalFuncPar_T=0.0d0
    EvalList=0.0d0
    WList=0.0d0
    write(*,*) "Hey"

    call reloadPrmsetTemp(E,30.0d0)
    P=>First
    do while(associated(P%next))
      CName = P%n
      if ((CName(1:1)=="k").and.(P%en<0.0d0)) then
        EvalFuncPar_T=1000.0d0
        return
      end if
      P=>P%next
    end do

    NT=40
    call getDataPar_T(NT*300,10,E,Phos,ATP,ATPase,Afs,Bfs,Q10,Phos0,Peri04)

    if (Fuck==1) then
      EvalFuncPar_T=100000.0d0
      Fuck=0
      return
    end if

! EvalList
! 1:PhosCAB30c,14:26c, 15:28c,12:PhosCAB38c,16:40c,
! 2:PhosCA, 3:PhosC(U,T),
! 4:NucC, 5:NucCA,11:NucCB, 6:NucCAB,
! 7:ATPaseC, 8:ATPaseCA, 9:ATPaseCB, 10:ATPaseCAB,
! 13:Bfree,17:Q10 t=0, 18:Q10 t=inf
! 20:Period@Atot=2.4
    WList(1)  = 0.5d0
    WList(14) = 1.0d0
    WList(15) = 0.5d0
    WList(12) = 0.5d0
    WList(16) = 0.5d0
    WList(2)  = 0.1d0
    WList(3)  = 0.25d0
    WList(4)  = 0.0d0
    WList(5)  = 0.0d0
    WList(11) = 0.0d0
    WList(6)  = 0.0d0
    WList(7)  = 0.1d0
    WList(8)  = 0.0d0
    WList(9)  = 0.0d0
    WList(10) = 0.0d0
    WList(24) = 0.1d0
    WList(13) = 0.1d0
    WList(17) = 0.2d0
    WList(18) = 0.05d0
    WList(19) = 1.0d0
    WList(20) = 0.3d0
    WList(21) = 0.05d0
    WList(22) = 0.02d0
    WList(23) = 0.01d0
    WList(25) = 0.30d0

    EvalFuncPar_T=0.0d0
    TotalW=0.0d0
    do ITime=0,NT
      do I=0,1
      do J=0,1
        EvalFuncPar_T = EvalFuncPar_T &
          + PhosExp(I,J,3)%w(ITime) * (PhosExp(I,J,3)%d(ITime)-Phos(3,I,J,ITime*30))**2
        TotalW = TotalW + PhosExp(I,J,3)%w(ITime)
      end do
      end do
    end do
    EvalList(1) = EvalFuncPar_T/TotalW

    EvalFuncPar_T=0.0d0
    TotalW=0.0d0
    do ITime=0,NT
      do I=0,1
      do J=0,1
        EvalFuncPar_T = EvalFuncPar_T &
          + PhosExp(I,J,7)%w(ITime) * (PhosExp(I,J,7)%d(ITime)-Phos(4,I,J,ITime*30))**2
        TotalW = TotalW + PhosExp(I,J,7)%w(ITime)
      end do
      end do
    end do
    EvalList(12) = EvalFuncPar_T/TotalW

    EvalFuncPar_T=0.0d0
    TotalW=0.0d0
    do ITime=0,NT
      do I=0,1
      do J=0,1
        EvalFuncPar_T = EvalFuncPar_T &
          + PhosExp(I,J,1)%w(ITime) * (PhosExp(I,J,1)%d(ITime)-Phos(5,I,J,ITime*30))**2
        TotalW = TotalW + PhosExp(I,J,1)%w(ITime)
      end do
      end do
    end do
    EvalList(14) = EvalFuncPar_T/TotalW

    EvalFuncPar_T=0.0d0
    TotalW=0.0d0
    do ITime=0,NT
      do I=0,1
      do J=0,1
        EvalFuncPar_T = EvalFuncPar_T &
          + PhosExp(I,J,2)%w(ITime) * (PhosExp(I,J,2)%d(ITime)-Phos(6,I,J,ITime*30))**2
        TotalW = TotalW + PhosExp(I,J,2)%w(ITime)
      end do
      end do
    end do
    EvalList(15) = EvalFuncPar_T/TotalW

    EvalFuncPar_T=0.0d0
    TotalW=0.0d0
    do ITime=0,NT
      do I=0,1
      do J=0,1
        EvalFuncPar_T = EvalFuncPar_T &
          + PhosExp(I,J,8)%w(ITime) * (PhosExp(I,J,8)%d(ITime)-Phos(7,I,J,ITime*30))**2
        TotalW = TotalW + PhosExp(I,J,8)%w(ITime)
      end do
      end do
    end do
    EvalList(16) = EvalFuncPar_T/TotalW

    EvalFuncPar_T = (Phos(2,0,0,200)-0.22d0)**2 &
                + (Phos(2,0,1,200)-0.22d0)**2 &
                + (Phos(2,1,1,200)-0.38d0)**2 &
                + (Phos(2,1,0,200)-0.18d0)**2
    EvalList(2) = EvalFuncPar_T*0.25d0

    EvalList(3) = (Phos(0,1,1,120)**2+8.0d0*Phos(0,1,0,120)**2+Phos(0,0,1,480)**2)/10.0d0

    EvalList(4) = (0.5d0*(ATP(0,1,120)+ATP(0,2,120))-0.4d0)**2
    !EvalList(4) = ((ATP(0,1,120)-0.8d0)**2+ ATP(0,2,120)**2)*0.5d0
    EvalList(11) = (ATP(1,1,480)-0.5d0)**2
    EvalList(5) = (0.5d0*(ATP(2,1,120)+ATP(2,2,120))-0.70d0)**2

    TempData = 0.5d0*(ATP(3,1,0:1200)+ATP(3,2,0:1200))
    call getExtremum(TempData(300:600),Extremum)

    DataMax=TempData(600)
    DataMin=TempData(600)
    do I=1,10
      if (Extremum(I)%f==0) exit
      if (Extremum(I)%f==1) then
        DataMax=Extremum(I)%v
      else if (Extremum(I)%f==-1) then
        DataMin=Extremum(I)%v
      end if
    end do

    EvalList(6) = 0.5d0*((DataMax-0.8d0)**2+(DataMin-0.25d0)**2)

    TotalATPase = sum(ATPase,dim=2)/30.0d0

    !EvalList(7) = (TotalATPase(0,480)-0.48d0)**2
    EvalList(7) = (TotalATPase(0,480)-0.40d0)**2
    EvalList(8) = (TotalATPase(2,480)-0.60d0)**2
    EvalList(9) = (TotalATPase(1,480)-0.30d0)**2

    EvalList(24) &
      =(sum(TotalATPase(3,0:480))/sum(TotalATPase(2,0:480))-0.9d0)**2 &
      +(sum(TotalATPase(0,0:480))/sum(TotalATPase(2,0:480))-0.8d0)**2
    EvalList(24) = 0.5d0*EvalList(24)

    TempData = TotalATPase(3,0:1200)
    call getExtremum(TempData(300:600),Extremum)

    DataMax=TempData(600)
    DataMin=TempData(600)
    do I=1,10
      if (Extremum(I)%f==0) exit
      if (Extremum(I)%f==1) then
        DataMax=Extremum(I)%v
      else if (Extremum(I)%f==-1) then
        DataMin=Extremum(I)%v
      end if
    end do

    EvalList(10) = 0.5d0*((DataMax-0.93d0)**2+(DataMin-0.18d0)**2)


    TempData = (Bfs(3,0,0:1200)+Bfs(3,1,0:1200))/3.5d0
    call getExtremum(TempData(300:600),Extremum)

    DataMax=TempData(600)
    DataMin=TempData(600)
    do I=1,10
      if (Extremum(I)%f==0) exit
      if (Extremum(I)%f==1) then
        DataMax=Extremum(I)%v
      else if (Extremum(I)%f==-1) then
        DataMin=Extremum(I)%v
      end if
    end do

    EvalList(13) = (DataMin-0.38d0)**2


    EvalList(17) = min(Q10(1,0)-1.4d0,0.0d0)**2+min(Q10(2,0)-1.2d0,0.0d0)**2
    EvalList(17) = 0.5d0*EvalList(17)
    EvalList(18) = (Q10(1,1)-1.05d0)**2

    EvalList(19) = 2.0d0*(Phos0(1,0)-0.06d0)**2+(Phos0(0,1)-0.3d0)**2+(Phos0(1,1)-0.3d0)**2
    EvalList(19) = EvalList(19)/4.0d0

    TempData = Phos(1,0,0,0:1200)
    call getExtremum(TempData(300:1200),Extremum)

    Period=0.1d0*dble(Extremum(3)%i-Extremum(1)%i)
    if ( abs(Extremum(3)%v-Extremum(2)%v)<0.05d0) Period=0.0d0
    EvalList(20) = (Period/24.0d0-1.0d0)**2

    TempData = Phos(3,0,0,0:1200)
    call getExtremum(TempData(300:1200),Extremum)

    Period=0.1d0*dble(min(Extremum(2)%i-Extremum(1)%i,Extremum(3)%i-Extremum(2)%i))
    if ( abs(Extremum(3)%v-Extremum(2)%v)<0.05d0) Period=0.0d0
    EvalList(21) = (0.46d0-Period/25.0d0)**2

    EvalList(22) = (max(log(K2RFDUP%v/K2RFDUU%v)/log(10.0d0)+0.7d0,0.0d0)/4.0d0)**2

    Periods=0.0d0
    do I=1,5
      TempData = Phos(I+2,0,0,0:1200)
      call getExtremum(TempData(300:1200),Extremum)

      Period=0.1d0*dble(Extremum(3)%i-Extremum(1)%i)
      if ( abs(Extremum(3)%v-Extremum(2)%v)<0.05d0) Period=0.0d0
      Periods(I)=Period
    end do
    EvalList(23)=((maxval(Periods)-minval(Periods))/24.0d0)**2

    EvalList(25) = (min(0.0d0,Peri04/24.0d0-1.35d0))**2

    EvalDetails=EvalList**0.5d0
    WDetails=WList
    if (FEvalOut==1) then
      call writeEvalDetails(6)
    end if

    EvalFuncPar_T = (sum(EvalList*WList)/sum(WList))**0.5d0

  end function EvalFuncPar_T

  subroutine writeEvalDetails(Uni)
    integer,intent(in) :: Uni
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosCAB30",EvalDetails(1) ,WDetails(1)
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosCAB26",EvalDetails(14),WDetails(14)
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosCAB28",EvalDetails(15),WDetails(15)
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosCAB38",EvalDetails(12),WDetails(12)
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosCAB40",EvalDetails(16),WDetails(16)
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosCA   ",EvalDetails(2) ,WDetails(2)
    write(Uni,'(A10,ES13.5e2,F6.2)') "PhosC    ",EvalDetails(3) ,WDetails(3)
    write(Uni,'(A10,ES13.5e2,F6.2)') "NucC     ",EvalDetails(4) ,WDetails(4)
    write(Uni,'(A10,ES13.5e2,F6.2)') "NucCA    ",EvalDetails(5) ,WDetails(5)
    write(Uni,'(A10,ES13.5e2,F6.2)') "NucCB    ",EvalDetails(11),WDetails(11)
    write(Uni,'(A10,ES13.5e2,F6.2)') "NucCAB   ",EvalDetails(6) ,WDetails(6)
    write(Uni,'(A10,ES13.5e2,F6.2)') "ATPaseC  ",EvalDetails(7) ,WDetails(7)
    write(Uni,'(A10,ES13.5e2,F6.2)') "ATPaseCA ",EvalDetails(8) ,WDetails(8)
    write(Uni,'(A10,ES13.5e2,F6.2)') "ATPaseCB ",EvalDetails(9) ,WDetails(9)
    write(Uni,'(A10,ES13.5e2,F6.2)') "ATPaseCAB",EvalDetails(10),WDetails(10)
    write(Uni,'(A10,ES13.5e2,F6.2)') "ADP Prod ",EvalDetails(24),WDetails(24)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Bf       ",EvalDetails(13),WDetails(13)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Q10 0    ",EvalDetails(17),WDetails(17)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Q10 inf  ",EvalDetails(18),WDetails(18)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Phos0c   ",EvalDetails(19),WDetails(19)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Period@3A",EvalDetails(20),WDetails(20)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Tau p    ",EvalDetails(21),WDetails(21)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Mu pT/T  ",EvalDetails(22),WDetails(22)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Diff Peri",EvalDetails(23),WDetails(23)
    write(Uni,'(A10,ES13.5e2,F6.2)') "Perik1h04",EvalDetails(25),WDetails(25)
  end subroutine writeEvalDetails

  function EvalFuncPar(E)

    real(8) :: EvalFuncPar
    real(8),dimension(NList) :: EvalList,WList
    real(8),dimension(1:),intent(in) :: E
    type dataexp
      real(8),dimension(0:100) :: t=0.0d0
      real(8),dimension(0:100) :: d=0.0d0
      real(8),dimension(0:100) :: w=0.0d0
    end type
    type(dataexp),dimension(0:1,0:1),save :: PhosExp
    real(8),dimension(0:1,0:1,0:1,0:1,0:1200) :: Phos
    real(8),dimension(0:1,0:1,2,0:1200) :: ATP
    real(8),dimension(0:1,0:1,6,0:1200) :: ATPase
    real(8),dimension(0:1,0:1,1,0:1200) :: Afs
    real(8),dimension(0:1,0:1,0:1,0:1200) :: Bfs
    real(8),dimension(0:1,0:1,0:1200) :: TotalATPase
    real(8),dimension(0:1200) :: TempData
    type(ext),dimension(10) :: Extremum
    real(8) :: TotalW,Temp,DataMax,DataMin,Al
    integer :: I,J,ITime
    integer,save :: Flag=0


    if (Flag==0) then
      open(10,file="data/exp/murayama/phos_exp.dat",status="old")
      read(10,"()")
      do I=0,45
        read(10,"(F4.0,3X,F9.7,3(4X,F9.7))") &
          PhosExp(0,0)%t(I),PhosExp(1,1)%d(I),PhosExp(0,1)%d(I),PhosExp(1,0)%d(I),PhosExp(0,0)%d(I)
        PhosExp(0,1)%t(I)=PhosExp(0,0)%t(I)
        PhosExp(1,0)%t(I)=PhosExp(0,0)%t(I)
        PhosExp(1,1)%t(I)=PhosExp(0,0)%t(I)
      end do
      close(10)

      !ST
      PhosExp(0,0)%w(1:45)=5.0d0
      !SpT
      PhosExp(0,1)%w(1:45)=3.0d0
      !pSpT
      PhosExp(1,1)%w(1:45)=1.0d0
      !pST
      PhosExp(1,0)%w(1:45)=3.0d0

      do I=0,1
      do J=0,1
        PhosExp(I,J)%w(1:10)=PhosExp(I,J)%w(1:10)*5.0d0
        PhosExp(I,J)%w(21:30)=PhosExp(I,J)%w(21:30)*5.0d0
      end do
      end do

      Flag=1
    end if


    EvalFuncPar=0.0d0
    EvalList=0.0d0
    WList=0.0d0
    write(*,*) "Hey"

    call getDataPar(6000,10,E,Phos,ATP,ATPase,Afs,Bfs)

    if (Fuck==1) then
      EvalFuncPar=100000.0d0
      Fuck=0
      return
    end if

! EvalList
! 1:PhosCAB, 2:PhosCA, 3:PhosC(U,T), 4:NucC, 5:NucCA,
! 6:NucCAB, 7:ATPaseC, 8:ATPaseCA, 9:ATPaseCB, 10:ATPaseCAB
    WList(1)  = 1.0d0
    WList(12) = 0.0d0
    WList(2)  = 0.1d0
    WList(3)  = 0.1d0
    WList(4)  = 0.05d0
    WList(5)  = 0.05d0
    WList(11) = 0.01d0
    WList(6)  = 0.1d0
    WList(7)  = 0.05d0
    WList(8)  = 0.025d0
    WList(9)  = 0.0d0
    WList(10) = 0.025d0
    WList(13) = 0.1d0

    EvalFuncPar=0.0d0
    TotalW=0.0d0
    do ITime=0,30
      do I=0,1
      do J=0,1
        EvalFuncPar = EvalFuncPar + PhosExp(I,J)%w(ITime) * (PhosExp(I,J)%d(ITime)-Phos(1,1,I,J,ITime*20))**2
        TotalW = TotalW + PhosExp(I,J)%w(ITime)
      end do
      end do
    end do

    EvalList(1) = EvalFuncPar/TotalW


    EvalFuncPar = (Phos(1,0,0,0,200)-0.22d0)**2 &
                + (Phos(1,0,0,1,200)-0.22d0)**2 &
                + (Phos(1,0,1,1,200)-0.38d0)**2 &
                + (Phos(1,0,1,0,200)-0.18d0)**2
    EvalList(2) = EvalFuncPar*0.25d0

    EvalList(3) = (Phos(0,0,1,1,120)**2+Phos(0,0,1,0,120)**2+Phos(0,0,0,1,480)**2)/3.0d0

    !EvalList(4) = (0.5d0*(ATP(0,0,1,120)+ATP(0,0,2,120))-0.38d0)**2
    EvalList(4) = ((ATP(0,0,1,120)-0.8d0)**2+ ATP(0,0,2,120)**2)*0.5d0
    EvalList(11) = (ATP(0,1,1,480)-0.5d0)**2
    EvalList(5) = (0.5d0*(ATP(1,0,1,120)+ATP(1,0,2,120))-0.70d0)**2

    TempData = 0.5d0*(ATP(1,1,1,0:1200)+ATP(1,1,2,0:1200))
    call getExtremum(TempData(300:600),Extremum)

    DataMax=TempData(600)
    DataMin=TempData(600)
    do I=1,10
      if (Extremum(I)%f==0) exit
      if (Extremum(I)%f==1) then
        DataMax=Extremum(I)%v
      else if (Extremum(I)%f==-1) then
        DataMin=Extremum(I)%v
      end if
    end do

    EvalList(6) = 0.5d0*((DataMax-0.8d0)**2+(DataMin-0.25d0)**2)

    TotalATPase = sum(ATPase,dim=3)/30.0d0

    EvalList(7) = (TotalATPase(0,0,480)-0.48d0)**2
    EvalList(8) = (TotalATPase(1,0,480)-0.60d0)**2
    EvalList(9) = (TotalATPase(0,1,480)-0.30d0)**2

    TempData = TotalATPase(1,1,0:1200)
    call getExtremum(TempData(300:600),Extremum)

    DataMax=TempData(600)
    DataMin=TempData(600)
    do I=1,10
      if (Extremum(I)%f==0) exit
      if (Extremum(I)%f==1) then
        DataMax=Extremum(I)%v
      else if (Extremum(I)%f==-1) then
        DataMin=Extremum(I)%v
      end if
    end do

    EvalList(10) = 0.5d0*((DataMax-0.93d0)**2+(DataMin-0.18d0)**2)


    TempData = (Bfs(1,1,0,0:1200)+Bfs(1,1,1,0:1200))/3.5d0
    call getExtremum(TempData(300:600),Extremum)

    DataMax=TempData(600)
    DataMin=TempData(600)
    do I=1,10
      if (Extremum(I)%f==0) exit
      if (Extremum(I)%f==1) then
        DataMax=Extremum(I)%v
      else if (Extremum(I)%f==-1) then
        DataMin=Extremum(I)%v
      end if
    end do

    EvalList(13) = (DataMin-0.38d0)**2

    if (FEvalOut==1) then
      EvalDetails=EvalList**0.5d0
      write(*,'(A10,ES13.5e2)') "PhosCAB  ",EvalList(1)**0.5
      !write(*,'(A10,ES13.5e2)') "PhosCAB- ",EvalList(12)**0.5
      write(*,'(A10,ES13.5e2)') "PhosCA   ",EvalList(2)**0.5
      write(*,'(A10,ES13.5e2)') "PhosC    ",EvalList(3)**0.5
      write(*,'(A10,ES13.5e2)') "NucC     ",EvalList(4)**0.5
      write(*,'(A10,ES13.5e2)') "NucCA    ",EvalList(5)**0.5
      write(*,'(A10,ES13.5e2)') "NucCB    ",EvalList(11)**0.5
      write(*,'(A10,ES13.5e2)') "NucCAB   ",EvalList(6)**0.5
      write(*,'(A10,ES13.5e2)') "ATPaseC  ",EvalList(7)**0.5
      write(*,'(A10,ES13.5e2)') "ATPaseCA ",EvalList(8)**0.5
      write(*,'(A10,ES13.5e2)') "ATPaseCB ",EvalList(9)**0.5
      write(*,'(A10,ES13.5e2)') "ATPaseCAB",EvalList(10)**0.5
      write(*,'(A10,ES13.5e2)') "Bf       ",EvalList(13)**0.5
    end if

    EvalFuncPar = (sum(EvalList*WList)/sum(WList))**0.5d0

  end function EvalFuncPar

  subroutine initEqConst()

    integer :: I,J,K,L,I1
    integer,dimension(0:7) :: I2
    integer,dimension(3) :: IC2
    real(8) :: Temp

    K1B(0)=K1gsB%v
    K1B(1)=k1fsB%v

    do I=0,6
      k1IA(I) = k1IA0%v**(dble(6-I)/6.0d0) * k1IA6%v**(dble(I)/6.0d0)
    end do


    do I=0,1715
      call decompInd(I,I1,I2)
      call I2toIC2(I2,IC2)
      K2RF(I) =   K2RFDUU%v & 
                 * (K2RFTUU%v/K2RFDUU%v)**(dble(IC2(1))/6.0d0) &
                 * (K2RFDPU%v/K2RFDUU%v)**(dble(IC2(2))/6.0d0) &
                 * (K2RFDUP%v/K2RFDUU%v)**(dble(IC2(3))/6.0d0) 
    end do

    Gam2A0= (1.0d0/K2FIA%v + K2RF/K2RIA%v)/(1.0d0+1.0d0/K2FIA%v + K2RF*(1.0d0+1.0d0/K2RIA%v))
    K2Atil=K2A%v/Gam2A0

  end subroutine initEqConst

  function GamR(Af)

    real(8),dimension(0:1715) :: GamR
    real(8),intent(in) :: Af
    real(8) :: TempR,TempF

    TempR=1.0d0+1.0d0/K2RIA%v*(1.0d0+Af/K2A%v)
    TempF=1.0d0+1.0d0/K2FIA%v*(1.0d0+Af/K2A%v)

    GamR=0.0d0
    GamR=K2RF*TempR/(K2RF*TempR+TempF)

  end function

  function Gam2A(Af)

    real(8),dimension(0:1715) :: Gam2A
    real(8),intent(in) :: Af
    real(8) :: TempR,TempF

    TempR=1.0d0+1.0d0/K2RIA%v*(1.0d0+Af/K2A%v)
    TempF=1.0d0+1.0d0/K2FIA%v*(1.0d0+Af/K2A%v)

    Gam2A=0.0d0
    Gam2A=1.0d0-(K2RF+1.0d0)/(K2RF*TempR+TempF)

  end function

  subroutine getAfBf(X,XMat,Af,Bf)

    real(8),dimension(-1:1),intent(in) :: X
    real(8),dimension(1716,14),intent(in) :: XMat
    real(8),intent(out) :: Af
    real(8),dimension(0:1),intent(out) :: Bf
    real(8),dimension(0:13) :: C1
    real(8),dimension(0:1715) :: C2
    real(8) :: AfPrev,Tf, TfPrev,Eps,DAf,DTf,Temp,Tf0,Tf1
    integer :: I,J,K,L,N,Nmax,M

    C1(0:13) = sum(XMat,Dim=1)
    C2(0:1715) = sum(XMat,Dim=2)

    Eps=0.00001d0

    Af=0.0d0
    Bf=0.0d0
    Tf=0.0d0
    AfPrev=Af
    TfPrev=Tf

    NMax=2000


    !!! Tf = gsBf/Kgs + fsBf/Kfs
    do N=1,NMax
      if (B0%v<0.001d0) exit

      if (N==1) then
        Tf0=0.0d0
        Tf1=B0%v*sum(1.0d0/K1B)
        do M=1,10
          if (FuncT0(0.5d0*(Tf0+Tf1))<0.0d0) then
            Tf0=0.5d0*(Tf0+Tf1)
          else
            Tf1=0.5d0*(Tf0+Tf1)
          end if
        end do
        Tf=0.5d0*(Tf0+Tf1)
      end if

      TfPrev=Tf
      Temp= FuncT0(Tf)
      DTf = -Temp/FuncT1(Tf)

      Tf = TfPrev + DTf
!      write(*,*) Tf
      if ((abs(1.0d0-TfPrev/Tf)<Eps)) then
!         write(*,*) N
         exit
      end if
      if (N==NMax) then
        write(*,*) "FuckT"
        Fuck=1
      end if
    end do

    Bf = X(0:1)/(1.0d0+f0(Tf)/K1B)


    do N=1,NMax

      AfPrev=Af
      Temp=FuncA(0)
      DAf = -Temp/FuncA(1)

      Af = AfPrev + DAf

      if ((abs(Temp)<Eps)) then
!         write(*,*) N
         exit
      end if
      if (N==NMax) then
        Fuck=1
        write(*,*) "FuckA"
      end if
    end do
  contains

    function funcA(Ind)
      real(8) :: funcA
      integer,intent(in) :: Ind
      real(8) :: Temp

      Temp = 6.0d0*C0%v*sum(C1(7:13))
      Temp = Temp * (Bf(1)/K1B(1)/(1.0d0+sum(Bf/K1B)))**6

      if (Ind==0) funcA = - A0%v + Af + Temp*Af/(K1A%v+Af) + C0%v*sum(C2*Af/(K2Atil+Af))
      if (Ind==1) funcA = 1.0d0 + Temp*K1A%v/(K1A%v+Af)**2 + C0%v*sum(C2*K2Atil/(K2Atil+Af)**2)
      !if (Ind==0) funcA = - A0%v + Af + X(1)*Af/(K1A%v+Af) + C0%v*sum(C2*Af/(K2Atil+Af))
      !if (Ind==1) funcA = 1.0d0 + X(1)*K1A%v/(K1A%v+Af)**2 + C0%v*sum(C2*K2Atil/(K2Atil+Af)**2)

    end function funcA

    function f0(t)
      real(8) :: f0
      real(8),intent(in) :: t
      real(8) :: Temp

      f0=0.0d0
      Temp=1.0d0+t
      f0=6.0d0*C0%v*sum(C1(7:13))/Temp

    end function f0

    function f1(t)
      real(8) :: f1
      real(8),intent(in) :: t
      real(8) :: Temp

      f1=0.0d0
      Temp=1.0d0+t
      f1= - 6.0d0*C0%v*sum(C1(7:13)) / (Temp)**2

    end function f1

    function funcT0(t)
      real(8) :: funcT0
      real(8),intent(in) :: t

      funcT0=0.0d0
      funcT0=t-sum(X(0:1)/(K1B+f0(t)))

    end function funcT0

    function funcT1(t)
      real(8) :: funcT1
      real(8),intent(in) :: t

      funcT1=0.0d0
      funcT1= 1.0d0+f1(t)*sum(X(0:1)/(K1B+f0(t))**2)

    end function funcT1

  end subroutine getAfBf

  subroutine getComplex(X,XMat,Comp)

    real(8),dimension(-1:1),intent(in) :: X
    real(8),dimension(1716,14),intent(in) :: XMat
    real(8),dimension(0:6,0:6,0:1),intent(out) :: Comp
    real(8) :: Af
    real(8),dimension(0:1) :: Bf,Ba
    real(8),dimension(0:1715,0:13) :: CMat
    real(8),dimension(0:1715) :: C2Bu,C2Ba
    real(8),dimension(0:6) :: Alpha1,Beta1
    real(8) :: Temp
    real(8),dimension(0:1,0:1715) :: Alpha2
    integer :: I,J,K,L,N,Nmax,M

    CMat(0:1715,0:13) = XMat
    C2Bu = sum(CMat(0:1715,0:6) ,Dim=2)
    C2Ba = sum(CMat(0:1715,7:13),Dim=2)
    call getAfBf(X,XMat,Af,Bf)

    Temp=Af/(K1A%v+Af)
    Alpha1(0)= 1.0d0 * (1.0d0-Temp)**6 * Temp**0
    Alpha1(1)= 6.0d0 * (1.0d0-Temp)**5 * Temp**1
    Alpha1(2)=15.0d0 * (1.0d0-Temp)**4 * Temp**2
    Alpha1(3)=20.0d0 * (1.0d0-Temp)**3 * Temp**3
    Alpha1(4)=15.0d0 * (1.0d0-Temp)**2 * Temp**4
    Alpha1(5)= 6.0d0 * (1.0d0-Temp)**1 * Temp**5
    Alpha1(6)= 1.0d0 * (1.0d0-Temp)**0 * Temp**6

    Temp=Bf(1)/K1B(1)/(1.0d0+sum(Bf/K1B))
    Beta1(0)= 1.0d0 * (1.0d0-Temp)**6 * Temp**0
    Beta1(1)= 6.0d0 * (1.0d0-Temp)**5 * Temp**1
    Beta1(2)=15.0d0 * (1.0d0-Temp)**4 * Temp**2
    Beta1(3)=20.0d0 * (1.0d0-Temp)**3 * Temp**3
    Beta1(4)=15.0d0 * (1.0d0-Temp)**2 * Temp**4
    Beta1(5)= 6.0d0 * (1.0d0-Temp)**1 * Temp**5
    Beta1(6)= 1.0d0 * (1.0d0-Temp)**0 * Temp**6

    Alpha2(0,0:1715) = K2Atil/(K2Atil+Af)
    Alpha2(1,0:1715) = Af/(K2Atil+Af)

    Comp=0.0d0
    Comp(0,0,0)=Comp(0,0,0)+sum(C2Bu*Alpha2(0,0:1715))
    Comp(0,0,1)=Comp(0,0,1)+sum(C2Bu*Alpha2(1,0:1715))
    Ba(0)=sum(C2Ba*Alpha2(0,0:1715))
    Ba(1)=sum(C2Ba*Alpha2(1,0:1715))
    do I=0,6
    do J=0,6
      Temp=Beta1(J)
      if (J==6) then
        Temp = Temp*Alpha1(I)
      else if (I/=0) then
        Temp = 0.0d0
      end if
      Comp(I,J,0)=Comp(I,J,0)+Temp*Ba(0)
      Comp(I,J,1)=Comp(I,J,1)+Temp*Ba(1)
    end do 
    end do 


  end subroutine getComplex

  function ExNucRatio(Time,Ind)
    real(8) :: ExNucRatio
    real(8),intent(in) :: Time
    integer,intent(in) :: Ind

    ExNucRatio=1.0d0

    if ((Time>=12.0d0).and.(mod(int(Time)-4*Ind+36,24)<12)) ExNucRatio=0.5d0

  end function ExNucRatio

!!-------------------------------------------------------------------------------
!! OED of the model. X(0):gsKaiB, X(1):fsKaiB,X(-1):time.
!!-------------------------------------------------------------------------------
  subroutine RHS_Model(X,XMat,RHS,RHSMat)

    real(8),dimension(-1:1),intent(in) :: X
    real(8),dimension(1716,14),intent(in) :: XMat
    real(8),dimension(-1:1),intent(out) :: RHS
    real(8),dimension(1716,14),intent(out) :: RHSMat
    real(8),dimension(1716,7) :: TempRHSMat

    real(8) :: Af,Gam,kbp,kbm,GamB
    real(8),dimension(0:1) :: Bf
    real(8),dimension(0:1715)::GammaR,Gamma2A,k1exA,k1exI,k2ex,k2Sp,k2Sm,k2Tp,k2Tm
    real(8),dimension(1716,7) :: Mk1exI,Mk1exA

    real(8),dimension(4800) :: Mk2h,Mk2ex
    real(8),dimension(2838) :: Mk2Sp,Mk2Sm,Mk2Tp,Mk2Tm
    character(len=1),dimension(6) :: matdescra=(/"G"," "," ","F"," "," "/)

    integer :: I

    call getAfBf(X,XMat,Af,Bf)

    RHS=0.0d0
    RHSMat=0.0d0
    RHS(-1) = 1.0d0

    Gam = Bf(0)
    kbp = (X(0)-Gam)*kbpC%v+Gam*kbpF%v
    Gam = Bf(1)
    kbm = (X(1)-Gam)*kbmC%v+Gam*kbmF%v

    RHS(0) = - kbp  + kbm
    RHS(1) =   kbp  - kbm



    GammaR=GamR(Af)
    Gamma2A=Gam2A(Af)

    k1exA = k1eAF%v + GammaR*(k1eAR%v-k1eAF%v)
    k1exI = k1eIF%v + GammaR*(k1eIR%v-k1eIF%v)
    k1exA = k1exA*ratio
    k1exI = k1exI*ratio
    Mk1exA=0.0d0
    Mk1exI=0.0d0
    do I=1,6
      Mk1exI(1:1716,I) = dble(7-I)*k1exI(0:1715)
      Mk1exA(1:1716,I) = dble(7-I)*k1exA(0:1715)
    end do



    k2ex = k2eI%v + Gamma2A*(k2eA%v-k2eI%v)
    k2ex = k2ex*ratio
    k2Sp = k2SpF%v + GammaR*(k2SpR%v-k2SpF%v)
    k2Sm = k2SmF%v + GammaR*(k2SmR%v-k2SmF%v)
    k2Tp = k2TpF%v + GammaR*(k2TpR%v-k2TpF%v)
    k2Tm = k2TmF%v + GammaR*(k2TmR%v-k2TmF%v)

    Mk2h = k2h%v * M2h
    do I=1,4800
      Mk2ex(I) = k2ex(Col2ex(I)-1) * M2ex(I)
    end do
    do I=1,2838
      Mk2Sp(I) = k2Sp(Col2Sp(I)-1) * M2Sp(I)
    end do
    do I=1,2838
      Mk2Sm(I) = k2Sm(Col2Sm(I)-1) * M2Sm(I)
    end do
    do I=1,2838
      Mk2Tp(I) = k2Tp(Col2Tp(I)-1) * M2Tp(I)
    end do
    do I=1,2838
      Mk2Tm(I) = k2Tm(Col2Tm(I)-1) * M2Tm(I)
    end do

    do I=2,7
      RHSMat(1:1716,I)   = RHSMat(1:1716,I)   -dble(I-1)*k1h%v*XMat(1:1716,I)
      RHSMat(1:1716,I+7) = RHSMat(1:1716,I+7) -dble(I-1)*k1h%v*XMat(1:1716,I+7)
      RHSMat(1:1716,I-1) = RHSMat(1:1716,I-1) +dble(I-1)*k1h%v*XMat(1:1716,I)
      RHSMat(1:1716,I+6) = RHSMat(1:1716,I+6) +dble(I-1)*k1h%v*XMat(1:1716,I+7)
    end do

    TempRHSMat = Mk1exI*XMat(1:1716,1:7)
    RHSMat(1:1716,1:6) = RHSMat(1:1716,1:6)-TempRHSMat(1:1716,1:6)
    RHSMat(1:1716,2:7) = RHSMat(1:1716,2:7)+TempRHSMat(1:1716,1:6)
    TempRHSMat = Mk1exA*XMat(1:1716,8:14)
    RHSMat(1:1716,8:13) = RHSMat(1:1716,8:13)-TempRHSMat(1:1716,1:6)
    RHSMat(1:1716,9:14) = RHSMat(1:1716,9:14)+TempRHSMat(1:1716,1:6)

    GamB=(1.0d0+(k1AIB6%v/k1AIB0%v)**(1.0d0/6.0d0)*sum(X(0:1)/K1B(0:1)))/(1.0d0+sum(X(0:1)/K1B(0:1)))
    GamB=GamB**6
    do I=1,7
      RHSMat(1:1716,I)   = RHSMat(1:1716,I)   - k1IA(I-1) * XMat(1:1716,I)
      RHSMat(1:1716,I+7) = RHSMat(1:1716,I+7) + k1IA(I-1) * XMat(1:1716,I)
      RHSMat(1:1716,I+7) = RHSMat(1:1716,I+7) - GamB*k1AIB0%v * XMat(1:1716,I+7)
      RHSMat(1:1716,I  ) = RHSMat(1:1716,I  ) + GamB*k1AIB0%v * XMat(1:1716,I+7)
    end do

    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2h,Col2h,Row2h(1:1716),Row2h(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2ex,Col2ex,Row2ex(1:1716),Row2ex(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Sp,Col2Sp,Row2Sp(1:1716),Row2Sp(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Sm,Col2Sm,Row2Sm(1:1716),Row2Sm(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Tp,Col2Tp,Row2Tp(1:1716),Row2Tp(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Tm,Col2Tm,Row2Tm(1:1716),Row2Tm(2:1717),XMat,1716,1.0d0,RHSMat,1716)

  end subroutine RHS_Model

  subroutine getStationary(XMat)

    real(8),dimension(-1:1) :: X
    real(8),dimension(1716,14),intent(out) :: XMat

    real(8) :: Af,Gam,kbp,kbm,GamB
    real(8),dimension(0:1) :: Bf
    real(8),dimension(0:1715)::GammaR,Gamma2A,k1exA,k1exI,k2ex,k2Sp,k2Sm,k2Tp,k2Tm
    real(8),dimension(1716,7) :: Mk1exI,Mk1exA

    real(8),dimension(4800) :: Mk2h,Mk2ex
    real(8),dimension(2838) :: Mk2Sp,Mk2Sm,Mk2Tp,Mk2Tm
    real(8),dimension(14388) :: MC2Rate1=0.0d0,MC2Rate2=0.0d0
    integer,dimension(14388) :: ColC2Rate1=0,ColC2Rate2=0
    integer,dimension(1717) :: RowC2Rate1=0,RowC2Rate2=0
    real(8),dimension(201432) :: MFullC2Rate=0.0d0
    integer,dimension(201432) :: ColFullC2Rate=0
    integer,dimension(24025) :: RowFullC2Rate=0
    real(8),dimension(266640) :: MatRate=0.0d0,MatRate2=0.0d0
    integer,dimension(266640) :: ColRate=0,ColRate2=0
    integer,dimension(24025) ::  RowRate=0,RowRate2=0
    real(8),dimension(65208) :: MatC1ND=0.0d0
    integer,dimension(65208) :: ColC1ND=0
    integer,dimension(24025) :: RowC1ND=0
    real(8),dimension(24024) :: MatC1D=0.0d0
    integer,dimension(24024) :: ColC1D=0
    integer,dimension(24025) :: RowC1D=0
    real(8),dimension(145332) :: MatL=0.0d0
    integer,dimension(145332) :: ColL=0
    integer,dimension(24025)  :: RowL=0
    real(8),dimension(121308) :: MatU=0.0d0
    integer,dimension(121308) :: ColU=0
    integer,dimension(24025)  :: RowU=0
    real(8),dimension(24024) :: XVec=0.0d0,XVec2=0.0d0

    integer :: I,J,Info,NZL,NZU,NZ


    X=0.0d0
    Af=0.0d0
    Bf=0.0d0

    GammaR=GamR(Af)
    Gamma2A=Gam2A(Af)

    k1exA = k1eAF%v + GammaR*(k1eAR%v-k1eAF%v)
    k1exI = k1eIF%v + GammaR*(k1eIR%v-k1eIF%v)

    k2ex = k2eI%v + Gamma2A*(k2eA%v-k2eI%v)
    k2Sp = k2SpF%v + GammaR*(k2SpR%v-k2SpF%v)
    k2Sm = k2SmF%v + GammaR*(k2SmR%v-k2SmF%v)
    k2Tp = k2TpF%v + GammaR*(k2TpR%v-k2TpF%v)
    k2Tm = k2TmF%v + GammaR*(k2TmR%v-k2TmF%v)

    Mk2h = k2h%v * M2h
    do I=1,4800
      Mk2ex(I) = k2ex(Col2ex(I)-1) * M2ex(I)
    end do
    do I=1,2838
      Mk2Sp(I) = k2Sp(Col2Sp(I)-1) * M2Sp(I)
    end do
    do I=1,2838
      Mk2Sm(I) = k2Sm(Col2Sm(I)-1) * M2Sm(I)
    end do
    do I=1,2838
      Mk2Tp(I) = k2Tp(Col2Tp(I)-1) * M2Tp(I)
    end do
    do I=1,2838
      Mk2Tm(I) = k2Tm(Col2Tm(I)-1) * M2Tm(I)
    end do

    !!!C2
    call mkl_dcsradd('N',0,0,1716,1716, Mk2h,Col2h,Row2h,1.0d0,&
      Mk2ex,Col2ex,Row2ex,MC2Rate1,ColC2Rate1,RowC2Rate1,14388,Info)
    call mkl_dcsradd('N',0,0,1716,1716, Mk2Sp,Col2Sp,Row2Sp,1.0d0,&
      MC2Rate1,ColC2Rate1,RowC2Rate1,MC2Rate2,ColC2Rate2,RowC2Rate2,14388,Info)
    call mkl_dcsradd('N',0,0,1716,1716, Mk2Sm,Col2Sm,Row2Sm,1.0d0,&
      MC2Rate2,ColC2Rate2,RowC2Rate2,MC2Rate1,ColC2Rate1,RowC2Rate1,14388,Info)
    call mkl_dcsradd('N',0,0,1716,1716, Mk2Tp,Col2Tp,Row2Tp,1.0d0,&
      MC2Rate1,ColC2Rate1,RowC2Rate1,MC2Rate2,ColC2Rate2,RowC2Rate2,14388,Info)
    call mkl_dcsradd('N',0,0,1716,1716, Mk2Tm,Col2Tm,Row2Tm,1.0d0,&
      MC2Rate2,ColC2Rate2,RowC2Rate2,MC2Rate1,ColC2Rate1,RowC2Rate1,14388,Info)
    do I=1,14
      MFullC2Rate((I-1)*14388+1:I*14388)=MC2Rate1
      ColFullC2Rate((I-1)*14388+1:I*14388)=ColC2Rate1+(I-1)*1716
      RowFullC2Rate((I-1)*1716+1:I*1716)=RowC2Rate1(1:1716)+(I-1)*14388
    end do
    RowFullC2Rate(24025)=RowC2Rate1(1717)+13*14388


    !!!C1
    GamB=(1.0d0+(k1AIB6%v/k1AIB0%v)**(1.0d0/6.0d0)*sum(X(0:1)/K1B(0:1)))/(1.0d0+sum(X(0:1)/K1B(0:1)))
    GamB=GamB**6


    NZ=1
    do I=1,7
      do J=1,1716
        RowC1ND((I-1)*1716+J) = NZ

        if (I>1) then
          MatC1ND(NZ)=dble(8-I)*k1exI(J-1)
          ColC1ND(NZ)=(I-2)*1716+J
          NZ=NZ+1
        end if
        if (I<7) then
          MatC1ND(NZ)=dble(I)*k1h%v
          ColC1ND(NZ)=(I)*1716+J
          NZ=NZ+1
        end if
        MatC1ND(NZ)=GamB*k1AIB0%v
        ColC1ND(NZ)=(I+6)*1716+J
        NZ=NZ+1
      end do
    end do

    do I=1,7
      do J=1,1716
        RowC1ND((I+6)*1716+J) = NZ

        MatC1ND(NZ)=k1IA(I-1)
        ColC1ND(NZ)=(I-1)*1716+J
        NZ=NZ+1
        if (I>1) then
        if (k1exA(J-1)/=0.0d0) then
          MatC1ND(NZ)=dble(8-I)*k1exA(J-1)
          ColC1ND(NZ)=(I+5)*1716+J
          NZ=NZ+1
        end if
        end if
        if (I<7) then
          MatC1ND(NZ)=dble(I)*k1h%v
          ColC1ND(NZ)=(I+7)*1716+J
          NZ=NZ+1
        end if
      end do
    end do
    RowC1ND(24025) = NZ


    do I=1,24025
      if (I<24025) ColC1D(I)=I
      RowC1D(I)=I
    end do

    MatC1D=0.0d0
    do I=1,NZ-1
      MatC1D(ColC1ND(I))=MatC1D(ColC1ND(I))-MatC1ND(I)
    end do
    call mkl_dcsradd('N',0,0,24024,24024, MatC1D,ColC1D,RowC1D,1.0d0,&
      MatC1ND(1:NZ-1),ColC1ND(1:NZ-1),RowC1ND,MatRate2,ColRate2,RowRate2,266640,Info)
    call mkl_dcsradd('N',0,0,24024,24024,MFullC2Rate,ColFullC2Rate,RowFullC2Rate,1.0d0,&
      MatRate2,ColRate2,RowRate2,MatRate,ColRate,RowRate,266640,Info)

    !!!Decomp to L and U triangle matrixes
    NZU=1
    NZL=1
    do I=1,24024
      RowL(I)=NZL
      RowU(I)=NZU
      do NZ=RowRate(I),RowRate(I+1)-1
        J=ColRate(NZ)
        if (I<J) then
          MatU(NZU)=MatRate(NZ)
          ColU(NZU)=J
          NZU=NZU+1
        else
          MatL(NZL)= - MatRate(NZ)
          ColL(NZL)=J
          NZL=NZL+1
        end if
      end do
      if (I==24024) then
        RowL(I+1)=NZL
        RowU(I+1)=NZU
      end if
    end do
    !write(*,*) maxval(MatL),maxval(MatU)

    !!! Gauss-Seidel
    XVec(24024)=1.0d0
    do I=1,1000
      call mkl_dcsrgemv('N',24024,MatU,RowU,ColU,XVec,XVec2)
      call mkl_dcsrtrsv('L','N','N',24024,MatL,RowL,ColL,XVec2,XVec)
      XVec=XVec/sum(XVec)
      !write(*,*) I,maxval(XVec),minval(XVec)
    end do



    XMat=0.0d0
    do I=1,14
      XMat(1:1716,I)=XVec((I-1)*1716+1:I*1716)
    end do


  end subroutine getStationary


  subroutine getATPase(X,XMat,ATPase)

    real(8),dimension(-1:1),intent(in) :: X
    real(8),dimension(1716,14),intent(in) :: XMat
    real(8),dimension(6),intent(out) :: ATPase
    real(8),dimension(1716,14) :: RHSMat

    real(8) :: Af
    real(8),dimension(0:1) :: Bf
    real(8),dimension(0:1715)::GammaR,Gamma2A,k2Sp,k2Sm,k2Tp,k2Tm

    real(8),dimension(4800) :: Mk2h
    real(8),dimension(2838) :: Mk2Sp,Mk2Sm,Mk2Tp,Mk2Tm
    character(len=1),dimension(6) :: matdescra=(/"G"," "," ","F"," "," "/)

    integer :: I

    call getAfBf(X,XMat,Af,Bf)


    GammaR=GamR(Af)
    Gamma2A=Gam2A(Af)


    k2Sp = k2SpF%v + GammaR*(k2SpR%v-k2SpF%v)
    k2Sm = k2SmF%v + GammaR*(k2SmR%v-k2SmF%v)
    k2Tp = k2TpF%v + GammaR*(k2TpR%v-k2TpF%v)
    k2Tm = k2TmF%v + GammaR*(k2TmR%v-k2TmF%v)

    Mk2h = k2h%v * M2h
    do I=1,2838
      Mk2Sp(I) = k2Sp(Col2Sp(I)-1) * M2Sp(I)
    end do
    do I=1,2838
      Mk2Sm(I) = k2Sm(Col2Sm(I)-1) * M2Sm(I)
    end do
    do I=1,2838
      Mk2Tp(I) = k2Tp(Col2Tp(I)-1) * M2Tp(I)
    end do
    do I=1,2838
      Mk2Tm(I) = k2Tm(Col2Tm(I)-1) * M2Tm(I)
    end do

    Mk2h=abs(Mk2h)
    Mk2Sp=abs(Mk2Sp)
    Mk2Sm=abs(Mk2Sm)
    Mk2Tp=abs(Mk2Tp)
    Mk2Tm=abs(Mk2Tm)

    RHSMat=0.0d0
    do I=2,7
      RHSMat(1:1716,I)   = RHSMat(1:1716,I)   +dble(I-1)*k1h%v*XMat(1:1716,I)
      RHSMat(1:1716,I+7) = RHSMat(1:1716,I+7) +dble(I-1)*k1h%v*XMat(1:1716,I+7)
      RHSMat(1:1716,I-1) = RHSMat(1:1716,I-1) +dble(I-1)*k1h%v*XMat(1:1716,I)
      RHSMat(1:1716,I+6) = RHSMat(1:1716,I+6) +dble(I-1)*k1h%v*XMat(1:1716,I+7)
    end do
    ATPase(1)=sum(RHSMat)/12.0d0


    RHSMat=0.0d0
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2h,Col2h,Row2h(1:1716),Row2h(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    ATPase(2)=sum(RHSMat)/12.0d0

    RHSMat=0.0d0
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Sp,Col2Sp,Row2Sp(1:1716),Row2Sp(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    ATPase(3)=sum(RHSMat)/12.0d0

    RHSMat=0.0d0
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Sm,Col2Sm,Row2Sm(1:1716),Row2Sm(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    ATPase(4)=-sum(RHSMat)/12.0d0

    RHSMat=0.0d0
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Tp,Col2Tp,Row2Tp(1:1716),Row2Tp(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    ATPase(5)=sum(RHSMat)/12.0d0

    RHSMat=0.0d0
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2Tm,Col2Tm,Row2Tm(1:1716),Row2Tm(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    ATPase(6)=-sum(RHSMat)/12.0d0

    ATPase = ATPase * 24.0d0

  end subroutine getATPase

  subroutine getADPRelease(X,XMat,ADPRelease)

    real(8),dimension(-1:1),intent(in) :: X
    real(8),dimension(1716,14),intent(in) :: XMat
    real(8),dimension(2),intent(out) :: ADPRelease
    real(8),dimension(1716,14) :: RHSMat
    real(8),dimension(1716,7) :: TempRHSMat

    real(8) :: Af,Gam,kbp,kbm,GamB
    real(8),dimension(0:1) :: Bf
    real(8),dimension(0:1715)::GammaR,Gamma2A,k1exA,k1exI,k2ex,k2Sp,k2Sm,k2Tp,k2Tm
    real(8),dimension(1716,7) :: Mk1exI,Mk1exA

    real(8),dimension(4800) :: Mk2h,Mk2ex
    real(8),dimension(2838) :: Mk2Sp,Mk2Sm,Mk2Tp,Mk2Tm
    character(len=1),dimension(6) :: matdescra=(/"G"," "," ","F"," "," "/)

    integer :: I

    call getAfBf(X,XMat,Af,Bf)

    GammaR=GamR(Af)
    Gamma2A=Gam2A(Af)

    k1exA = k1eAF%v + GammaR*(k1eAR%v-k1eAF%v)
    k1exI = k1eIF%v + GammaR*(k1eIR%v-k1eIF%v)
    k1exA = k1exA*ratio
    k1exI = k1exI*ratio
    Mk1exA=0.0d0
    Mk1exI=0.0d0
    do I=1,6
      Mk1exI(1:1716,I) = dble(7-I)*k1exI(0:1715)
      Mk1exA(1:1716,I) = dble(7-I)*k1exA(0:1715)
    end do

    k2ex = k2eI%v + Gamma2A*(k2eA%v-k2eI%v)
    k2ex = k2ex*ratio
    do I=1,4800
      Mk2ex(I) = k2ex(Col2ex(I)-1) * M2ex(I)
    end do
    Mk2ex=abs(Mk2ex)

    ADPRelease(1) = sum(Mk1exI*XMat(1:1716,1:7 ))&
                  + sum(Mk1exA*XMat(1:1716,8:14))

    RHSMat=0.0d0
    call mkl_dcsrmm("n",1716,14,1716,1.0d0,matdescra,&
                    Mk2ex,Col2ex,Row2ex(1:1716),Row2ex(2:1717),XMat,1716,1.0d0,RHSMat,1716)
    ADPRelease(2)=sum(RHSMat)*0.5d0

    ADPRelease = ADPRelease/6.0d0*24.0d0

  end subroutine getADPRelease


  subroutine resetNuc(XMat)
    real(8),dimension(1716,14),intent(inout) :: XMat
    real(8),dimension(0:1715) :: C2
    integer :: I,I1,J
    integer,dimension(0:7) :: I2,J2

    C2(0:1715)=sum(XMat(1:1716,1:14),dim=2)

    do I=0,1715
      call decompInd(I,I1,I2)
      J2(4)=I2(4)+I2(0)
      J2(5)=I2(5)+I2(1)
      J2(6)=I2(6)+I2(2)
      J2(7)=I2(7)+I2(3)
      J2(0)=0
      J2(1)=0
      J2(2)=0
      J2(3)=0
      J=XInd(0,J2)

      if (I/=J) then
        C2(J)=C2(J)+C2(I)
        C2(I)=0.0d0
      end if
    end do

    XMat=0.0d0
    XMat(1:1716,7)=k1AIB0%v/(k1AIB0%v+k1IA(6))*C2(0:1715)
    XMat(1:1716,14)=k1IA(6)/(k1AIB0%v+k1IA(6))*C2(0:1715)

  end subroutine resetNuc


  subroutine initMasks()

    integer :: I,I1
    integer,dimension(0:7) :: I2

    do I=0,1715
      call decompInd(I,I1,I2)

      PhosMask(I,0,0) = dble(I2(0)+I2(4))/6.0d0
      PhosMask(I,0,1) = dble(I2(1)+I2(5))/6.0d0
      PhosMask(I,1,1) = dble(I2(3)+I2(7))/6.0d0
      PhosMask(I,1,0) = dble(I2(2)+I2(6))/6.0d0

      NucPhosMask(I,0,0,0) = dble(I2(0))/6.0d0
      NucPhosMask(I,0,0,1) = dble(I2(1))/6.0d0
      NucPhosMask(I,0,1,1) = dble(I2(3))/6.0d0
      NucPhosMask(I,0,1,0) = dble(I2(2))/6.0d0
      NucPhosMask(I,1,0,0) = dble(I2(4))/6.0d0
      NucPhosMask(I,1,0,1) = dble(I2(5))/6.0d0
      NucPhosMask(I,1,1,1) = dble(I2(7))/6.0d0
      NucPhosMask(I,1,1,0) = dble(I2(6))/6.0d0

      C2ATPMask(I) = dble(I2(4)+I2(5)+I2(6)+I2(7))/6.0d0
    end do

    do I=0,6
      C1ATPMask(I) = dble(I)/6.0d0
    end do

  end subroutine initMasks

  subroutine initRateMats()

    integer :: I,J,J1,K,Info
    integer,dimension(0:7) :: I2,J2
    integer,dimension(3) :: JC2

    do J=0,1715
      call decompInd(J,J1,J2)
      call I2toIC2(J2,JC2)

      if (JC2(1)/=0) then
        UncoM2h(J+1,J+1) = UncoM2h(J+1,J+1) - dble(JC2(1))
        do K=4,7
          if (J2(K)/=0) then
            I2=J2
            I2(K)=I2(K)-1
            I2(K-4)=I2(K-4)+1
            I = XInd(0,I2)
            UncoM2h(I+1,J+1) = UncoM2h(I+1,J+1) + dble(J2(K))
          end if
        end do
      end if

      if (JC2(1)/=6) then
        UncoM2ex(J+1,J+1) = UncoM2ex(J+1,J+1) - dble(6-JC2(1))
        do K=0,3
          if (J2(K)/=0) then
            I2=J2
            I2(K)=I2(K)-1
            I2(K+4)=I2(K+4)+1
            I = XInd(0,I2)
            UncoM2ex(I+1,J+1) = UncoM2ex(I+1,J+1) + dble(J2(K))
          end if
        end do
      end if

      if ((J2(4)+J2(5))/=0) then
        UncoM2Sp(J+1,J+1) = UncoM2Sp(J+1,J+1) -dble((J2(4)+J2(5)))
        do K=4,5
          if (J2(K)/=0) then
            I2=J2
            I2(K)=I2(K)-1
            I2(K-2)=I2(K-2)+1
            I = XInd(0,I2)
            UncoM2Sp(I+1,J+1) = UncoM2Sp(I+1,J+1) + dble(J2(K))
          end if
        end do
      end if

      if ((J2(2)+J2(3))/=0) then
        UncoM2Sm(J+1,J+1) = UncoM2Sm(J+1,J+1) - dble((J2(2)+J2(3)))
        do K=2,3
          if (J2(K)/=0) then
            I2=J2
            I2(K)=I2(K)-1
            I2(K+2)=I2(K+2)+1
            I = XInd(0,I2)
            UncoM2Sm(I+1,J+1) = UncoM2Sm(I+1,J+1) + dble(J2(K))
          end if
        end do
      end if

      if ((J2(4)+J2(6))/=0) then
        UncoM2Tp(J+1,J+1) = UncoM2Tp(J+1,J+1) - dble((J2(4)+J2(6)))
        do K=4,6,2
          if (J2(K)/=0) then
            I2=J2
            I2(K)=I2(K)-1
            I2(K-3)=I2(K-3)+1
            I = XInd(0,I2)
            UncoM2Tp(I+1,J+1) = UncoM2Tp(I+1,J+1) + dble(J2(K))
          end if
        end do
      end if

      if ((J2(1)+J2(3))/=0) then
        UncoM2Tm(J+1,J+1) = UncoM2Tm(J+1,J+1) - dble((J2(1)+J2(3)))
        do K=1,3,2
          if (J2(K)/=0) then
            I2=J2
            I2(K)=I2(K)-1
            I2(K+3)=I2(K+3)+1
            I = XInd(0,I2)
            UncoM2Tm(I+1,J+1) = UncoM2Tm(I+1,J+1) + dble(J2(K))
          end if
        end do
      end if
    end do


    call mkl_ddnscsr((/0,1,1,2,4800,1/),1716,1716,UncoM2h, 1716,M2h, Col2h, Row2h, Info)
    call mkl_ddnscsr((/0,1,1,2,4800,1/),1716,1716,UncoM2ex,1716,M2ex,Col2ex,Row2ex,Info)
    call mkl_ddnscsr((/0,1,1,2,2838,1/),1716,1716,UncoM2Sp,1716,M2Sp,Col2Sp,Row2Sp,Info)
    call mkl_ddnscsr((/0,1,1,2,2838,1/),1716,1716,UncoM2Sm,1716,M2Sm,Col2Sm,Row2Sm,Info)
    call mkl_ddnscsr((/0,1,1,2,2838,1/),1716,1716,UncoM2Tp,1716,M2Tp,Col2Tp,Row2Tp,Info)
    call mkl_ddnscsr((/0,1,1,2,2838,1/),1716,1716,UncoM2Tm,1716,M2Tm,Col2Tm,Row2Tm,Info)

!    write(*,*) Row2h
!
!
!    write(*,*) minval(Col2h),maxval(Col2h)
!    write(*,*) minval(Col2ex),maxval(Col2ex)
!    write(*,*) minval(Col2Sp),maxval(Col2Sp)
!    write(*,*) minval(Col2Sm),maxval(Col2Sm)
!    write(*,*) minval(Col2Tp),maxval(Col2Tp)
!    write(*,*) minval(Col2Tm),maxval(Col2Tm)
  end subroutine initRateMats




  subroutine getExtremum(Dat,Extremum)

    real(8),dimension(:),intent(in) :: Dat
    type(ext),dimension(10),intent(out) :: Extremum
    integer :: I,Cnt
    real(8) :: RFlag

    Cnt=0
    do I=lbound(Dat,1)+1,ubound(Dat,1)-1
      RFlag = (Dat(I+1)-Dat(I))*(Dat(I)-Dat(I-1))

      if (RFlag<0.0d0) then
        Cnt=Cnt+1
        if (Dat(I)>Dat(I-1)) then
          Extremum(Cnt)%f=1
        else
          Extremum(Cnt)%f=-1
        end if
        Extremum(Cnt)%i=I
        Extremum(Cnt)%v=Dat(I)
      end if

      if (Cnt==10) exit
    end do

  end subroutine getExtremum

function getConf(X,XMat)

  real(8),dimension(0:1,0:1) :: getConf
  real(8),dimension(-1:1),intent(in) :: X
  real(8),dimension(1716,14),intent(in) :: XMat
  real(8),dimension(0:1715) :: C2_1I,C2_1A,GammaR
  real(8) :: Af
  real(8),dimension(0:1) :: Bf
  integer :: I,J

  getConf=0.0d0
  call getAfBf(X,XMat,Af,Bf)
  GammaR=GamR(Af)

  C2_1I(0:1715)=sum(XMat(1:1716,1:7),dim=2)
  C2_1A(0:1715)=sum(XMat(1:1716,8:14),dim=2)

  getConf(0,1)=sum(C2_1I*GammaR)
  getConf(1,1)=sum(C2_1A*GammaR)
  getConf(0,0)=sum(C2_1I)-getConf(0,1)
  getConf(1,0)=sum(C2_1A)-getConf(1,1)

end function getConf

function getPhos(XMat)

  real(8),dimension(0:1,0:1) :: getPhos
  real(8),dimension(1716,14),intent(in) :: XMat
  real(8),dimension(0:1715) :: C2
  integer :: I,J

  getPhos=0.0d0

  C2(0:1715)=sum(XMat(1:1716,1:14),dim=2)
  do I=0,1
  do J=0,1
    getPhos(I,J) = sum(PhosMask(0:1715,I,J)*C2(0:1715))
  end do
  end do

end function getPhos

function getNucPhos(XMat)

  real(8),dimension(0:1,0:1,0:1) :: getNucPhos
  real(8),dimension(1716,14),intent(in) :: XMat
  real(8),dimension(0:1715) :: C2
  integer :: I,J,K

  getNucPhos=0.0d0

  C2(0:1715)=sum(XMat(1:1716,1:14),dim=2)
  do K=0,1
  do I=0,1
  do J=0,1
    getNucPhos(K,I,J) = sum(NucPhosMask(0:1715,K,I,J)*C2(0:1715))
  end do
  end do
  end do

end function getNucPhos


function getATP(XMat)

  real(8),dimension(1:2) :: getATP
  real(8),dimension(1716,14),intent(in) :: XMat
  real(8),dimension(0:6) :: C1
  real(8),dimension(0:13) :: C1Full
  real(8),dimension(0:1715) :: C2

  getATP=0.0d0

  C1Full(0:13)=sum(XMat(1:1716,1:14),dim=1)
  C1=C1Full(0:6)+C1Full(7:13)
  getATP(1) = sum(C1ATPMask(0:6)*C1(0:6))
  C2(0:1715)=sum(XMat(1:1716,1:14),dim=2)
  getATP(2) = sum(C2ATPMask(0:1715)*C2(0:1715))

end function getATP

function getC1ATP(XMat)

  real(8),dimension(0:6) :: getC1ATP
  real(8),dimension(1716,14),intent(in) :: XMat
  real(8),dimension(0:13) :: C1Full

  getC1ATP=0.0d0

  C1Full(0:13)=sum(XMat(1:1716,1:14),dim=1)
  getC1ATP=C1Full(0:6)+C1Full(7:13)

end function getC1ATP


subroutine I2toIC2(I2,IC2)

  integer,dimension(0:7),intent(in) :: I2
  integer,dimension(3),intent(out) :: IC2

  IC2=0

  IC2(1)=I2(4)+I2(5)+I2(6)+I2(7)
  IC2(2)=I2(2)+I2(3)+I2(6)+I2(7)
  IC2(3)=I2(1)+I2(3)+I2(5)+I2(7)

end subroutine I2toIC2

subroutine initIndTable()

  integer :: N,R,XI,I,I1,J,K,Temp
  integer,dimension(0:7) :: I2


  do XI=0,1715

    N  = mod(XI,1716)

    R=6

    do I=0,6
      do J=0,R
        if (N<HTable(7-I,R-J)) then
          IndTable(I,XI)=J
          R=R-J
          exit
        else
          N=N-HTable(7-I,R-J)
        end if
      end do
    end do

    IndTable(7,XI)=6-sum(IndTable(0:6,XI))

  end do

end subroutine


subroutine decompInd(XI,I1,I2)

  integer,intent(in) :: XI
  integer,intent(out) :: I1
  integer,dimension(0:7),intent(out) :: I2

  I1 = XI/1716
  I2(0:7) = IndTable(0:7,mod(XI,1716) )

end subroutine decompInd

subroutine initHTable()

  integer :: I,J

  HTable(1,0:6)=1
  do I=2,8
  do J=0,6
    HTable(I,J)=sum(HTable(I-1,0:J))
  end do
  end do

end subroutine initHTable


function XInd(I1,I2)

  integer :: XInd
  integer,intent(in) :: I1
  integer,dimension(0:7),intent(in) :: I2

  integer :: R,I,J


  R=6
  XInd=0
  do I=0,6
    do J=0,I2(I)-1
      XInd = XInd + HTable(7-I,R-J)
    end do
    R=R-I2(I)
  end do

  XInd = XInd + I1*1716
end function XInd

!--------------------------------------------------------------------------------
! Runge-Kutta 4th
!--------------------------------------------------------------------------------
  subroutine R_Runge4thCustomized(X,XMat,H)

     real(8),dimension(-1:1),intent(inout) :: X
     real(8),dimension(1716,14),intent(inout) :: XMat
     real(8),intent(in) :: H

     real(8),dimension(-1:1) :: FA,FB,FC,FD
     real(8),dimension(1716,14) :: FMatA,FMatB,FMatC,FMatD

     call RHS_Model(X,XMat,FA,FMatA)
     call RHS_Model(X+0.5d0*H*FA,XMat+(0.5d0*H)*FMatA,FB,FMatB)
     call RHS_Model(X+0.5d0*H*FB,XMat+(0.5d0*H)*FMatB,FC,FMatC)
     call RHS_Model(X+H*FC,XMat+H*FMatC,FD,FMatD)
     X  = X + (H/6.0d0) * (FA+FD+FB+FB+FC+FC)
     XMat  = XMat + (H/6.0d0) * (FMatA+FMatD+FMatB+FMatB+FMatC+FMatC)

  end subroutine R_Runge4thCustomized
!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
! N, X  : First N elements in X is treated as variables to be optimized.
! Func  : Evaluation function.
! Table : N+1 square matrix, where Table(I,1:N) is the coordinate of the Ith
!         vertex Xi and Table(I,N+1) = Func(Xi). Table(I,N+1) are sorted in
!         ascending order (i.e. Table(1,N+1) is the minimal). 
!--------------------------------------------------------------------------------
  subroutine initNelderMead(N,X,Func,Table)

    integer,intent(in) :: N
    real(8),dimension(1:),intent(inout) :: X
    interface
      function Func(Y)
        real(8) :: Func
        real(8),dimension(1:),intent(in) :: Y
      end function Func
    end interface
    real(8),dimension(:,:),intent(inout) :: Table
    real(8),dimension(NVM+1) :: XTemp
    real(8) :: lam,R
    integer :: I,J,Swaps

    lam=-0.2d0
    XTemp=0.0d0
    Table=0.0d0

    do I=1,N+1
      call random_number(R)
      Table(I,1:N)=X(1:N)
      !if (I/=N+1) Table(I,I) = Table(I,I)+lam*(2.0d0*R-1.0d0)
      !if (I/=N+1) Table(I,I) = Table(I,I)+lam*(R-0.5d0+sign(0.5d0,R-0.5d0))
      if (I/=N+1) Table(I,I) = Table(I,I)+lam*(R+1.0d0)
      !if (I/=N+1) Table(I,I) = Table(I,I)+lam
      Table(I,N+1) = Func(Table(I,1:N))
      write(*,*) I, Table(I,N+1)
    end do

    do I=1,N
      Swaps=0
      do J=N+1,I+1,-1
        if (Table(J,N+1)<Table(J-1,N+1)) then
          XTemp(1:N+1)=Table(J,1:N+1)
          Table(J,1:N+1)=Table(J-1,1:N+1)
          Table(J-1,1:N+1)=XTemp(1:N+1)
          Swaps=Swaps+1
        end if
      end do
      if (Swaps==0) exit
    end do

    X(1:N)=Table(1,1:N)

  end subroutine initNelderMead
!--------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine forwardNelderMead(N,X,Func,Table)


    integer,intent(in) :: N
    real(8),dimension(1:),intent(inout) :: X
    interface
      function Func(Y)
        real(8) :: Func
        real(8),dimension(1:),intent(in) :: Y
      end function Func
    end interface
    real(8),dimension(:,:),intent(inout) :: Table
    real(8),dimension(NVM+1) :: XR,XE,XC,XTemp
    real(8),dimension(NVM) :: X0
    real(8),save :: alp=1.0d0,gam=2.0d0,rho=-0.5d0,sig=0.5d0
    integer :: I,J,Swaps


    X0=0.0d0
    do I=1,N
      X0(1:N)=X0(1:N)+Table(I,1:N)
    end do
    X0=X0/dble(N)

    XR=0.0d0
    XR(1:N)=X0(1:N)+alp*(X0(1:N)-Table(N+1,1:N))
    XR(N+1)=Func(XR(1:N))

    if ( (XR(N+1)>=Table(1,N+1)) .and. (XR(N+1)<Table(N+1,N+1)) ) then

      Table(N+1,1:N+1)=XR(1:N+1)

    else if (XR(N+1)<Table(1,N+1)) then

      call writeEvalDetails(6)

      XE=0.0d0
      XE(1:N)=X0(1:N)+gam*(X0(1:N)-Table(N+1,1:N))
      XE(N+1)=Func(XE(1:N))

      if (XE(N+1)<XR(N+1)) then
        call writeEvalDetails(6)
        Table(N+1,1:N+1)=XE(1:N+1)
      else
        Table(N+1,1:N+1)=XR(1:N+1)
      end if

    else

      XC=0.0d0
      XC(1:N)=X0(1:N)+rho*(X0(1:N)-Table(N+1,1:N))
      XC(N+1)=Func(XC(1:N))

      if (XC(N+1)<Table(N+1,N+1)) then
        Table(N+1,1:N+1)=XC(1:N+1)
      else
        write(*,*) "Shrink"
        do I=2,N+1
          Table(I,1:N)=Table(1,1:N)+sig*(Table(I,1:N)-Table(1,1:N))
          Table(I,N+1) = Func(Table(I,1:N))
        end do
      end if

    end if

    do I=1,N
      Swaps=0
      do J=N+1,I+1,-1
        if (Table(J,N+1)<Table(J-1,N+1)) then
          XTemp(1:N+1)=Table(J,1:N+1)
          Table(J,1:N+1)=Table(J-1,1:N+1)
          Table(J-1,1:N+1)=XTemp(1:N+1)
          Swaps=Swaps+1
        end if
      end do
      if (Swaps==0) exit
    end do

    X(1:N)=Table(1,1:N)

!    write(*,*) Ite, Table(1,N+1)

  end subroutine forwardNelderMead
!-------------------------------------------------------------------------------

  subroutine setOMPMKL(N,M)
    integer,intent(in) :: N
    integer,intent(in),optional :: M
    integer :: K
    if (present(M)) then
      K=M
    else
      !K=omp_get_num_procs()/N
      !write(*,*) omp_get_max_threads()
      K=omp_get_max_threads()/N
    end if
    write(*,*) "OMP:", N, ", MKL:",K
    !$ call omp_set_num_threads(N)
    !$ call mkl_set_num_threads(K)
  end subroutine setOMPMKL


end module model
