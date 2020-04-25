program main

  !$ use omp_lib
  use model
  implicit none

  real(8),dimension(NVM) :: E=0.0d0,Eptb=0.0d0
  real(8),dimension(0:7,0:1,0:1,0:1,0:1,0:1200) :: Phos,Conf
  real(8),dimension(0:7,0:1,0:1,2,0:1200) :: ATP
  real(8),dimension(0:7,0:1,0:1,0:6,0:1200) :: C1ATP
  real(8),dimension(0:7,0:1,0:1,6,0:1200) :: ATPase
  real(8),dimension(0:7,0:1,0:1,1,0:1200) :: Afs,TempAfs
  real(8),dimension(0:7,0:1,0:1,0:1,0:1200) :: Bfs
  integer :: I,J,ITime,IA,IB,NT,ITe
  character(len=3) :: CComp
  character(len=2) :: CTe
  character(len=100) :: Filename

  call setOMPMKL(8)
  !$ call mkl_set_dynamic(0)
  !$ call omp_set_nested(1)
  !$ call omp_set_max_active_levels(2)

  call initHTable()
  call initIndTable()
  call initRateMats()
  call initMasks()
  call initPrmset()

  Filename ="prmset.dat"
  call readPrmset(Filename,1)
  call initEqConst()


  ! Call subroutines in model.f90 as you want and print/analize data below.
  ! Folloing code is an example.

  call getDataPar_TOut(12000,10,E,Phos,ATP,ATPase,Afs,Bfs,Conf)

  do ITe=0,7
    IA=1
    IB=1
    NT=1200

    write(CTe,'(I2.2)') ITe*2+26
    CComp="c"
    if (IA==1) CComp = trim(CComp)//"a"
    if (IB==1) CComp = trim(CComp)//"b"

    open(10,file="phos"//CTe//"c.dat",status="unknown")
    do ITime=0,NT
      write(10,'(20(ES13.5e2))') 0.1d0*dble(ITime),Phos(ITe,IA,IB,0,0,ITime),Phos(ITe,IA,IB,0,1,ITime),&
        Phos(ITe,IA,IB,1,1,ITime),Phos(ITe,IA,IB,1,0,ITime)
    end do
    close(10)
  end do

end program main
