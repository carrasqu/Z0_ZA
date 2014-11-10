
! Kagome lattice XYZ model SSE program     
!-------------------------------------------------

!--------------------
module configuration
!--------------------

save

!lattice has L=3*lx*lysites


integer :: lx  
integer :: ly 
integer :: lz ! not used in this code
integer :: L  ! total number of sites or spins
integer :: nb ! number of plquettes
integer :: nver ! number of vertices
integer(4) :: nhl ! number of sites or spins again
integer(4) :: ndim ! number of sites in the bravais lattice of the kagome
integer(4) :: nbase ! number of lattice sites in the unit cell
integer(4) :: nindm,nindms,maxmulti,maxmultis ! dimensions of tables of neighbours
integer(4), allocatable :: imulti(:),ivicn(:,:,:),ivicns(:,:,:),imultis(:)

integer :: nh,nh2 ! order of the SSE expansion for the two replicas 
integer :: mm,mm2 ! max order considered in the SSE expansion (for each replica)
integer :: dd ! dimensionality
integer :: disseed ! seed for the disorder  configuration
integer :: nl ! number of loops per update
integer :: maxloop ! maximum length of the loop operators
integer :: mloopc ! coefficient for maximum loop length maxloop=mloopc*<n>
integer :: termal ! termalization call sign
integer :: totvisit
integer :: looplength
integer :: spino
integer :: nh0,nh20 ! for averages over the different processors
integer :: densiono ! measure density
integer :: info ! lapack info variable
integer :: inittype ! type of initialization (mf based=1 , random=0)
integer :: thermalization ! whether the simulation is thermalizing or not
integer :: full ! what part of nk to measure (0 only along x. 1 full nk)
integer :: directed
integer :: restart ! restart=1(restart) 0(start from scratch)
integer :: kstart ! number of times writeresults have been executed 

real(8) :: beta  ! inverse temperature 
real(8) :: aprob ! nb*beta
real(8) :: mu    ! chemical potential
real(8) :: cadd  ! sse constant makes all config postive definite 
real(8) :: delta ! bound of the disordered potential
real(8) :: t     ! "hopping" matrix element for S^+S^- +  S^-S^+
real(8) :: tp    ! "hopping" matrix element for S^+S^+ +  S^-S^-
real(8) :: Jz    ! SzSz term (or NN  interaction V for hardcore bosons )
real(8) :: sup   ! amplitude of the superlattice potential is present
real(8) :: vo    ! strength of the optical trap potential vo*(r)**2, r distance from the center of the trap
real(8) :: timeli,timelo,timeme,timew

real(8), allocatable :: epsdis(:) ! espdis(i)=mu+ epsilon_i : disorder configuration plus overall chemical potential
real(8), allocatable :: prob(:,:,:,:) ! probabilities for the loop updates
                       !prob(bond,vertex type(st bond p), entrance leg, exit leg) en principio.  

!real(8), allocatable :: wdatgf(:,:),wdat2gf(:,:)                        

real(8), allocatable :: vtex(:,:) ! weight of the vertices
real(8) :: NNN
real(8) :: z
real(8), allocatable :: gf(:) ! one body density matrix
!real(8), allocatable :: eigenval(:) ! eigenvalues of gf
!real(8), allocatable :: workgf(:)

real(kind=8), allocatable :: loc(:) ! local density
real(kind=8), allocatable :: loc0(:) ! local density reduced to the zeroth node


integer, allocatable :: tryloop(:,:)

real(8) :: densgf

integer, allocatable :: spin(:),spin2(:),tspin(:),tspin2(:) ! spin (-1 or 1 ) or (boson 0 or 1 respectively) configuration
integer, allocatable :: spin2AFTloop(:)
integer, allocatable :: bsites(:,:) ! bond sites which is the lattice information
integer, allocatable :: bbond(:,:) 
integer, allocatable :: sm(:),sm2(:),sml(:) ! operator string
integer, allocatable :: reA(:),reAi(:),reAipl(:),dl(:)
integer :: sizedl,whichreg ! whichreg=1 reA=reAi. whichreg=0 reA=reAipl 
LOGICAL            :: accepted 

integer, allocatable :: frstspin(:),frstspin2(:),frstspin3(:) ! temporary; helps constructing vtx
integer, allocatable :: lastspin(:),lastspin2(:),lastspin3(:) ! temporary; helps constructing vtx
integer, allocatable :: vtx(:) ! linked vertex storage
!integer, allocatable :: A(:,:,:)

integer, allocatable :: vtype(:),vtype2(:) ! type of vertex
integer, allocatable :: vtypel(:) ! temporary type of vertex during loop construction
 
integer, allocatable :: newvtx(:,:,:) ! newvtx(lentrance,lexit,vtype before the change)
integer, allocatable :: optype(:) ! 0 for diagonal 1 for off-diagonal operators,it is function of the vertex type vtype(:)
integer, allocatable :: legspin(:,:) ! value of the spin or boson density for a given leg and type of vertex


end module configuration

!----------------------!
 module measurementdata
!----------------------!
 save

 real(8) :: enrg1=0.d0
 real(8) :: enrg2=0.d0
 real(8) :: token(8)
 real(8) :: token0(8)
 real(8) :: data1(8)=0.d0
 real(8) :: data2(8)=0.d0
 real(kind=8) pi
 
 real(kind=8),allocatable :: kkx(:,:)     ! k values
 real(kind=8),allocatable :: ord(:,:)  ! coordinates
 real(kind=8),allocatable :: distances(:,:)  ! coordinates
 integer(4),allocatable :: cdistances(:)
 integer(4) indist                           ! number of distances considered
 integer(4) nqv                              ! number of momenta considered   
 integer(4),allocatable :: site_to_d(:,:) 
 real(kind=8),allocatable :: q(:,:)    !
 real(kind=8)ba1(2),ba2(2)
 integer :: lk
 real(kind=8) :: kmax
 complex(8) ima

 !renyi entropy measurement
 integer(4) Ni,Nipl


 end module measurementdata
!--------------------------!



!============================!
 program hcb_sse
!============================!
 use configuration
 use measurementdata
 use mpi
 implicit none
 
 !include 'mpif.h'
 integer(4) stats(MPI_STATUS_SIZE)

 integer :: nbins,msteps,isteps,iseedd,k,m,irstart,mstepsq

 !parallelization variables
 integer :: size,rank,ierr

 real(8) :: drand1,avt,avn, time1,time,avn0,avt0

 !parallel initialization


  call mpi_init(ierr)

  if (ierr.ne.mpi_success) write(6,*)'some problem initializing parallel'

  call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

  

 !input file reading. All processors read the file.
 open(unit=10,file='input.in',form='formatted',status='unknown')

 !read(10,*) dd
 dd=2 ! works only for the kagome
 
 read(10,*) lx 
 read(10,*) ly
 !read(10,*) lz 

 read(10,*) beta
 read(10,*) mu
 read(10,*) t
 read(10,*) tp
 read(10,*) Jz
 read(10,*) nbins
 read(10,*) msteps
 read(10,*) 
 read(10,*) isteps
 read(10,*) iseedd
 read(10,*) disseed
 read(10,*) delta
 read(10,*) cadd 
 read(10,*) sup
 read(10,*) vo
 read(10,*) densiono
 read(10,*) 
 read(10,*) inittype
 read(10,*) nqv
 read(10,*) kmax
 read(10,*) full
 read(10,*) directed
 read(10,*) mloopc
 read(10,*) restart 
 
 cadd=4.0+3.0*abs(mu)/4.0+3.0*abs(Jz)/4.0

 
 !sqmeas
 ! nbins: number of bins, msteps: monte carlo steps, isteps: thermalization steps,iseed: seed for the random generator,
 ! mu:chemical potential,disseed: seed for the disorder, sup=amplitude staggered potential
 ! vo: the strength of the parabolic potential. 
 ! densiono =1 measure loc density, 0 to avoid writing it
 ! inittype= 1 initialization based on the  meanfield density. 0 for random initialization ! lk: number of k points in the momentum distribution along x (lx/2 should be used if pbc are used and one wants the right quasimomenta)
 ! kmax : maximum kx value considered. (kmax should bi \pi if pbc are used to
 ! get the right quasimomenta)
 ! directed=0: directed loop updates with no bounces when possible. If not possible, minimal bounces. directed>0 : heat bath solution 
 ! mloopc coefficient for maximum loop length maxloop=mloopc*<n> 

 close(10)


! if(rank==0)then ! rank 0 writes the energy each msteps
! open(unit=13,file='fort.12',form='unformatted',status='unknown') ! energy at  each mc step to check binning
! end if

!initialization ------------------------------------------------------- 

 if(dd==1)then
  ly=1 
  lz=1
 elseif(dd==2)then
  lz=1   
 endif         

 z=2 

 pi=2.0d0*asin(1.0d0)
 ima=cmplx(0.0d0,1.0d0)

 call makelattice()  !construct the lattice
 
 call newvertex()  

 aprob=beta*nb

 if(directed==0)then
  call probupdates() ! construct the probabilities used in the loop updates. Directed loop updates
 elseif(directed>0)then
  call probupdates_heat_bath() ! heat bath solution for the updates
 end if


 call restarts(iseedd)

 ! this bit makes a different seed for each processor
  irstart=iseedd
  do k=1,rank
   irstart=ran(iseedd)*2**29
  enddo
  iseedd=2*irstart+1
 ! write(6,*) 'seed',  iseedd,rank

 call rand_init(iseedd) ! initializaation of random number 

 call initconf() ! initial random configuration


! write(6,*) 'configuration', spin, rank


! write(6,*)'epsdis',epsdis,rank
!----------------------------------------------------------------------

  nl=10
  maxloop=30000
  avn=0.0d0
  avt=0.0d0
  timeli=0.0d0
  timelo=0.0d0  
  timeme=0.0d0
  timew=0.0d0 

  allocate(tryloop(nl,4))
 ! Thermalization
 thermalization=1 

 do k=1,isteps

  do m=1,msteps

   call diagonalupdate()
  
   call loopupdate()
 
   !call getspin2()
   ! collects the maximum nh in order to enlarge the cutoff
   call mpi_allreduce(nh,nh0,1,mpi_integer4,mpi_max,mpi_comm_world,ierr)
   call mpi_allreduce(nh2,nh20,1,mpi_integer4,mpi_max,mpi_comm_world,ierr)
    
   call adjustcutoff(k)
 
   
   avn=avn+dble(nh+nh2)
   avt=avt+dble(totvisit)

  ! call check()

  end do 

  
    avn=avn/dble(msteps)
    avt=avt/dble(msteps)    

    call mpi_reduce(avn,avn0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(avt,avt0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

    if(rank==0)then
            
     avn0=avn0/dble(size)
     avt0=avt0/dble(size)
     write(6,*)'avn,avt', avn0,avt0,nl  

     if(2.0d0*avn0-avt0>0)then
      nl=nl+1
     else
      nl=nl-1      
      if(nl==0)nl=1  
     end if
     
     maxloop=int(mloopc*avn0)
     avn0=0.0d0
     avt0=0.0d0
    

    end if

    call mpi_bcast(nl,1,mpi_integer4,0,mpi_comm_world,ierr)
    call mpi_bcast(maxloop,1,mpi_integer4,0,mpi_comm_world,ierr) 

    avn=0.0d0
    avt=0.0d0

   
   ! call writeresults(msteps,k,size,rank,nbins)

   

   deallocate(tryloop)
   allocate(tryloop(nl,4))


 end do

 if(rank==0)then
  write(6,*)'Finished equilibration, M = ',mm, 'maxloop=',maxloop, 'number of loops',nl
 end if


! call restarts()

 thermalization=0
 
 !rN=0.0d0
 !iN=0.0d0
 !rN0=0.0d0
 !iN0=0.0d0 
 !cN=(0.0d0,0.0d0) 

 

 
! gfm=0.0d0
 whichreg=1 
 !------- mc run ------------
 do k=kstart,kstart-1+nbins
  
  
  do m=1,msteps

   call diagonalupdate()
  
   call loopupdate()

   call getspin2()

   ! whichreg=1 reA=reAi. whichreg=0 reA=reAipl

   if(whichreg==1)then
    call checkspins()   
    if(accepted)then
     whichreg=0
     reA=reAipl
    end if
   elseif(whichreg==0)then
    call checkspins()
    if(accepted)then
     whichreg=1
     reA=reAi
    end if 
   end if
  
   if(whichreg==1)then
    Ni=Ni+1
   elseif(whichreg==0)then
    Nipl=Nipl+1
   end if

!   write(6,*)'whicreg',whichreg,reA  

   if(looplength>=0)then
    call check(m)
   end if
 
   !!!!call mpi_reduce(nh,nh0,1,mpi_real8,mpi_sum,0,mpi_comm_world,ierr) 
 
   !!!!if(rank==0)then 
   !!!! write(13) 1.0d0,1.0d0, -(1.0d0/beta)*dble(nh0)/dble(size)+cadd*dble(nb)
   !!!!end if 

  end do
  
  
  call writeresults(msteps,k,size,rank,nbins,iseedd)

 
 end do


 write(6,*)'done'
 !---------------------------
 if(rank==0)then
  close(13)
 end if

 call mpi_finalize(ierr)
 stop

 end program hcb_sse
!================================!


subroutine makelattice() 
!------------------------------------------------------------------------------
! Constructs the list of sites bsites(1,b) and bsites(2,b)  bsites(3,b)  in plaquette b
!------------------------------------------------------------------------------

use configuration
use measurementdata
implicit none

integer is,x1,x2,y1,y2,z1,z2,i,bond,j,qv,typeb,wh,sub,k,kk
real(8) drand1,qx,qy,qz
complex(16) valuec




is=0 

! ndim number of lattice sites in the bravais lattice
! nbase number of sites in the unit cell.
! total number of sites L=ndim*nbase

read(9) ndim,nbase,nindm,nindms,maxmulti,maxmultis

L=nbase*ndim
nhl=L
nb=2*L/nbase ! number of triangular plaquettes ( works for the kagome lattice)



ALLOCATE(imulti(nindm))
ALLOCATE(imultis(nindms))
ALLOCATE(ivicn(nhl,maxmulti,nindm))
ALLOCATE(ivicns(nhl,maxmultis,nindms))
allocate(ord(L,dd))

read(9) imulti,imultis,ivicn,ivicns,ord


! momentum

if(full==1)nqv=ndim


allocate(bsites(3,nb)) ! 3 from having 3 spins per plaquette

! constructing the plquettes. Works for the kagome defined by the provided
! geometry.d but (hopefully later?) can be generalized.

bond=1 ! bond means plaquette in this code
j=1

do i=1,L
 
   
   if(mod(i,nbase)==1)then
    bsites(1,bond)=i
    bsites(2,bond)=ivicn(i,1,1)
    bsites(3,bond)=ivicn(i,2,1)   
    bond=bond+1
   elseif(mod(i,nbase)==2)then
    bsites(1,bond)=i
    bsites(2,bond)=ivicn(i,3,1)
    bsites(3,bond)=ivicn(i,4,1)   
    bond=bond+1 
   end if 

end do



call cdist() ! organizes the different distances in the torus (PBC)



allocate(epsdis(L))

      
 epsdis=mu




end subroutine makelattice 

!------------------------
!initial configuration
!-----------------------
 subroutine initconf()

 use configuration
 implicit none
 integer :: i,reading 
 real(8) :: drand1, deni,maxni

 allocate(spin(L),tspin(L),tspin2(L),spin2(L),reA(L),reAi(L),reAipl(L),spin2AFTloop(L))
 spin=0
 spin2=0
 tspin=0
 tspin2=0
 spin2AFTloop=0
 !open(unit=10,file='input.in',form='formatted',status='unknown')
 open(unit=39,file='bipart.dat',form='formatted',status='unknown')
 do i=1,L
  read(39,*)reading
  if(reading==0.or.reading==1)then
   reA(i)=reading
  else
   write(6,*) 'wrong input in bipartition file'
   stop
  end if

 end do 
 read(39,*)sizedl
 allocate(dl(sizedl))
 do i=1,sizedl
  read(39,*)reading
  if(reA(reading)/=0)then
   write(6,*) 'wrong input in bipartition file dl sites should be outside region A', reading,reA(reading)
   stop
  elseif(reA(reading)==0)then
   dl(i)=reading
  end if
 end do

 reAi=reA
 reAipl=reA
 do i=1,sizedl
  reAipl(dl(i))=1
 end do


 close(39) 
 
 do i=1,L
 ! random initialization
  spin(i)=2*int(2.*drand1())-1
 end do

 do i=1,L
  if(reA(i)>0)then
   spin2(i)=spin(i)
  elseif(reA(i)==0)then
   spin2(i)=2*int(2.*drand1())-1
  end if
 end do

 mm=5 !initial guess for mm the max order of the SSE expansion
 mm2=5 

 allocate(sm(mm),sm2(mm2),sml(mm+mm2))
 sm(:)=0 !all mm operators in the initial configuration are of the [0,0]  type => nh=0
 sm2(:)=0
 sml(:)=0 
 nh=0 
 nh2=0
 allocate(frstspin(L),frstspin2(L),frstspin3(L))
 allocate(lastspin(L),lastspin2(L),lastspin3(L))
 allocate(vtx(0:6*(mm+mm2)-1))
 !allocate(A(mm,L,2))
 allocate(vtype(mm),vtypel(mm+mm2),vtype2(mm2))
 end subroutine initconf
!-----------------------

 subroutine diagonalupdate()
 use configuration
 implicit none
 integer :: i, b,op,cc
 real(8) :: drand1, ratio,nran,error
 integer :: spint(3), typo,legt(0:5),vertex
 tspin=spin

 !vtype(:)=0

 do i=1,mm
  op=sm(i)
  if(op==0)then ! [0,0] operator found
   b=int(drand1()*nb+1)     !select a plaquette at random among the possible nb 

   spint(1)=spin(bsites(1,b))
   spint(2)=spin(bsites(2,b))
   spint(3)=spin(bsites(3,b))  

   ratio=cadd+ 0.5d0/z*(epsdis(bsites(1,b))*spin(bsites(1,b))+ &
             epsdis(bsites(2,b))*spin(bsites(2,b))+ epsdis(bsites(3,b))*spin(bsites(3,b))     )-&
             0.25*Jz*( spin(bsites(1,b))*spin(bsites(2,b))+ spin(bsites(1,b))*spin(bsites(3,b)) +&
              spin(bsites(2,b))*spin(bsites(3,b)))

   ratio=aprob*ratio/dble(mm-nh)   
 
  
   ! metropolis        
   nran=drand1()           
   if(ratio>=1.0d0 .or. ratio>=nran  )then
    
    sm(i)=2*b  !insert diagonal operator      
    nh=nh+1    ! increases the order of the SSE expansion 

   end if  

  elseif(mod(op,2)==0)then ! diagonal operator found
   !op=2*b
   b=op/2

   spint(1)=spin(bsites(1,b))
   spint(2)=spin(bsites(2,b))
   spint(3)=spin(bsites(3,b))  

   ratio=cadd+ 0.5d0/z*(epsdis(bsites(1,b))*spin(bsites(1,b))+ &
             epsdis(bsites(2,b))*spin(bsites(2,b))+ epsdis(bsites(3,b))*spin(bsites(3,b))     )-&
             0.25*Jz*(spin(bsites(1,b))*spin(bsites(2,b))+ spin(bsites(1,b))*spin(bsites(3,b)) +&
             spin(bsites(2,b))*spin(bsites(3,b)))

   ratio=dble(mm-nh+1)/(aprob*ratio)
        
       
         

   !metropolis      
   nran=drand1()
   if(ratio>=1.0d0 .or. ratio>=nran  )then

    sm(i)=0  !remove diagonal operator      
    nh=nh-1  ! decreases the order of the SSE expansion

   end if


  else ! off-diagonal operator found
          
   ! op=2*b+1  
   !propagate the state |alpha(p)> 

   b=op/2 ! (its integer part) 

   spint(1)=spin(bsites(1,b))
   spint(2)=spin(bsites(2,b))
   spint(3)=spin(bsites(3,b))
  
  ! write(6,*) 'should be zero', spint(1)-legspin(0,vtype(i))+ &
  !                             spint(2)-legspin(1,vtype(i))+ &
  !                             spint(3)-legspin(2,vtype(i))
  ! write(6,*) 'shoild be >8',vtype(i)

   spin(bsites(1,b))=legspin(3,vtype(i))
   spin(bsites(2,b))=legspin(4,vtype(i))
   spin(bsites(3,b))=legspin(5,vtype(i))

  end if 

  !construction of the types of vertices
  if(sm(i)>0)then
 
   cc=iabs(spint(1)-spin(bsites(1,b)))+iabs(spint(2)-spin(bsites(2,b)))+iabs(spint(3)-spin(bsites(3,b)))


   if(cc==0)then ! diagonal 
  
     legt(0)=spint(1)
     legt(1)=spint(2)
     legt(2)=spint(3)
     legt(3)=spint(1)
     legt(4)=spint(2)
     legt(5)=spint(3)
   
     call searchvtx(vertex,legt) 
     vtype(i)=vertex
  
   end if

  elseif(sm(i)==0)then

   vtype(i)=0

  end if 


 end do

!if(thermalization==0)then
!!=== handy for debugging=== 
! typo=0
! do i=1,mm
 
!  if(vtype(i).ne.0)then
! ! write(6,*)vtype(i) !,i

!  typo=typo+optype(vtype(i))
  
!  end if
! end do

!if(sum(tspin-spin)/=0)then 
! write(6,*) 'diagonal check', sum(tspin-spin), 'typo=',typo
!endif
!endif

!write(6,*)'   prop spin',spin

!vtype(:)=0

 do i=1,L
  if(reA(i)>0)then
  spin2(i)=spin(i)
  end if
 end do

 tspin2=spin2
 write(6,*)"key comparisons========================================" 
 write(6,*)'spin2 error get vs prop',sum(abs(spin2AFTloop-spin2))
 write(6,*)'propspin',spin
 write(6,*)'2AFTloop',spin2AFTloop
 write(6,*)'   spin2',spin2
 write(6,*)"*************************************************************************"

 do i=1,mm2
  op=sm2(i)
  if(op==0)then ! [0,0] operator found
   b=int(drand1()*nb+1)     !select a plaquette at random among the possible nb 

   spint(1)=spin2(bsites(1,b))
   spint(2)=spin2(bsites(2,b))
   spint(3)=spin2(bsites(3,b))  

   ratio=cadd+ 0.5d0/z*(epsdis(bsites(1,b))*spin2(bsites(1,b))+ &
             epsdis(bsites(2,b))*spin2(bsites(2,b))+ epsdis(bsites(3,b))*spin2(bsites(3,b))     )-&
             0.25*Jz*( spin2(bsites(1,b))*spin2(bsites(2,b))+ spin2(bsites(1,b))*spin2(bsites(3,b)) +&
              spin2(bsites(2,b))*spin2(bsites(3,b)))

   ratio=aprob*ratio/dble(mm2-nh2)   
 
  
   ! metropolis        
   nran=drand1()           
   if(ratio>=1.0d0 .or. ratio>=nran  )then
    
    sm2(i)=2*b  !insert diagonal operator      
    nh2=nh2+1    ! increases the order of the SSE expansion 

   end if  

  elseif(mod(op,2)==0)then ! diagonal operator found
   !op=2*b
   b=op/2

   spint(1)=spin2(bsites(1,b))
   spint(2)=spin2(bsites(2,b))
   spint(3)=spin2(bsites(3,b))  

   ratio=cadd+ 0.5d0/z*(epsdis(bsites(1,b))*spin2(bsites(1,b))+ &
             epsdis(bsites(2,b))*spin2(bsites(2,b))+ epsdis(bsites(3,b))*spin2(bsites(3,b))     )-&
             0.25*Jz*(spin2(bsites(1,b))*spin2(bsites(2,b))+ spin2(bsites(1,b))*spin2(bsites(3,b)) +&
             spin2(bsites(2,b))*spin2(bsites(3,b)))

   ratio=dble(mm2-nh2+1)/(aprob*ratio)

   !metropolis      
   nran=drand1()
   if(ratio>=1.0d0 .or. ratio>=nran  )then

    sm2(i)=0  !remove diagonal operator      
    nh2=nh2-1  ! decreases the order of the SSE expansion

   end if

  else ! off-diagonal operator found
          
   ! op=2*b+1  
   !propagate the state |alpha(p)> 

   b=op/2 ! (its integer part) 

   spint(1)=spin2(bsites(1,b))
   spint(2)=spin2(bsites(2,b))
   spint(3)=spin2(bsites(3,b))
  
  ! write(6,*) 'should be zero', spint(1)-legspin(0,vtype(i))+ &
  !                             spint(2)-legspin(1,vtype(i))+ &
  !                             spint(3)-legspin(2,vtype(i))
  ! write(6,*) 'shoild be >8',vtype(i)

   spin2(bsites(1,b))=legspin(3,vtype2(i))
   spin2(bsites(2,b))=legspin(4,vtype2(i))
   spin2(bsites(3,b))=legspin(5,vtype2(i))

  end if 

  !construction of the types of vertices
  if(sm2(i)>0)then
 
   cc=iabs(spint(1)-spin2(bsites(1,b)))+iabs(spint(2)-spin2(bsites(2,b)))+iabs(spint(3)-spin2(bsites(3,b)))


   if(cc==0)then ! diagonal 
  
     legt(0)=spint(1)
     legt(1)=spint(2)
     legt(2)=spint(3)
     legt(3)=spint(1)
     legt(4)=spint(2)
     legt(5)=spint(3)
   
     call searchvtx(vertex,legt) 
     vtype2(i)=vertex
  
   end if

  elseif(sm2(i)==0)then

  vtype2(i)=0

  end if 


 end do

 do i=1,L
  if(reA(i)>0)then
   spin(i)=spin2(i)
  end if 
 end do 

 error=0.0d0
 do i=1,L
  if(reA(i)==0)then
   error=error+abs(spin2(i)-tspin2(i))
  end if
 end do
 if(error>0.00000001)then
  write(6,*)'ERROR in spin2',error,whichreg
  write(6,*)' spin2',spin2
  write(6,*)'tspin2',tspin2
  write(6,*)'Pspin2',spin2AFTloop
  write(6,*)' spin',spin
  write(6,*)'tspin',tspin  
  write(6,*)'reA',reA 
  stop 
 end if
 spin2=tspin2



if(thermalization==0)then
!=== handy for debugging=== 
 typo=0
 do i=1,mm2
 
  if(vtype2(i).ne.0)then
 ! write(6,*)vtype(i) !,i

  typo=typo+optype(vtype2(i))
  
  end if
 end do

if(sum(tspin-spin)/=0)then 
 write(6,*) 'diagonal check', sum(tspin-spin), 'typo=',typo
endif
endif
! write(6,*)'order of exp n', nh

!write(6,*)'sm after',sm
! write(6,*)'type after',vtype
!================================= 
 end subroutine diagonalupdate
 !------------------------------------------------------------

 subroutine loopupdate()
 use configuration
 use measurementdata
 implicit none
 integer :: i,n,b,op,s1,s2,s3,v0,v1,v2,v3,p,ir,br,br1,k,aa,kk,cc,pbef,plast,io
 integer :: jo,j,li,le,ex,break,po,bo,igfo,lii,spinprop
 integer :: lf,pf,bf,pff,bff,lff, igff,igf,lo,typ,it,typoo
 integer :: crossed, go,lfi,jf,nobounce,dist,anylink
 real(8) drand1, ratio,nran !,time,time1



 densgf=0.0d0
 tryloop=-1
 
 do i=1,nl
   br=drand1()*(mm+mm2)+1
   br1=drand1()*L+1
   tryloop(i,1)=br1
   tryloop(i,2)=br 

   if(drand1()<0.5)then
    tryloop(i,3)=1 ! going up
   else
    tryloop(i,3)=0 ! going down         
   end if        

 end do
!------ construction of the linked vertex list -------

 !time=mclock()

 !A=-100 ! Given a random point in space-time, A tells you which is the closest leg going up or down in time
 ! if given the random point A is negative(-100), then that site has only trivial
 ! operators acting on it. 
 ! A is a huge memory waste, not used anymore.  

 anylink=-1
 frstspin(:)=-1
 lastspin(:)=-1
 frstspin2(:)=-1
 lastspin2(:)=-1
 frstspin3(:)=-1
 lastspin3(:)=-1

 vtx=-2
 do p=1,mm+mm2
 
  v0=6*(p-1)
  !p=v0/4+1 ! position of the operator
  
  if(p<=mm)then
   op=sm(p) ! type of operator
  elseif(p>mm)then
   op=sm2(p-mm)
  end if

  if(op/=0)then
   b=op/2  ! plaquette on which operator acts
   s1=bsites(1,b) 
   s2=bsites(2,b)
   s3=bsites(3,b)
   
   if(p<=mm)then
    if(reA(s1)>0)then
     v1=lastspin(s1)
    else
     v1=lastspin2(s1)
    end if 
    if(reA(s2)>0)then
     v2=lastspin(s2)
    else
     v2=lastspin2(s2)
    end if    
    if(reA(s3)>0)then
     v3=lastspin(s3)
    else
     v3=lastspin2(s3)
    end if  
   elseif(p>mm)then
    if(reA(s1)>0)then
     v1=lastspin(s1)
    else
     v1=lastspin3(s1)
    end if
    if(reA(s2)>0)then
     v2=lastspin(s2)
    else
     v2=lastspin3(s2)
    end if
    if(reA(s3)>0)then
     v3=lastspin(s3)
    else
     v3=lastspin3(s3)
    end if 
   end if
   !write(6,*) 'v1,s1,v2,s2,v0',v1,s1,v2,s2,v0,mm,mm2
   if(v1/=-1)then  ! found a link to v0 
    vtx(v1)=v0
    vtx(v0)=v1
    anylink=1 
   else
    if(p<=mm)then
     if(reA(s1)>0)then
      frstspin(s1)=v0  ! no link so it is first visit to that site
     else
      frstspin2(s1)=v0 
     endif
    elseif(p>mm)then
      if(reA(s1)>0)then
      frstspin(s1)=v0  ! no link so it is first visit to that site
     else
      frstspin3(s1)=v0
     endif 
    endif 
   end if 

   if(v2/=-1)then  ! found a link to v0+1 
    vtx(v2)=v0+1
    vtx(v0+1)=v2
    anylink=1
   else
    if(p<=mm)then
     if(reA(s2)>0)then
      frstspin(s2)=v0+1  ! no link so it is first visit to that site
     else
      frstspin2(s2)=v0+1
     endif
    elseif(p>mm)then
      if(reA(s2)>0)then
      frstspin(s2)=v0+1  ! no link so it is first visit to that site
     else
      frstspin3(s2)=v0+1
     endif
    endif
   end if
  
   if(v3/=-1)then  ! found a link to v0+1 
    pbef=v3/6+1
    vtx(v3)=v0+2
    vtx(v0+2)=v3
    anylink=1
   else
    if(p<=mm)then
     if(reA(s3)>0)then
      frstspin(s3)=v0+2  ! no link so it is first visit to that site
     else
      frstspin2(s3)=v0+2
     endif
    elseif(p>mm)then
      if(reA(s3)>0)then
      frstspin(s3)=v0+2  ! no link so it is first visit to that site
     else
      frstspin3(s3)=v0+2
     endif
    endif  
   end if
  
   ! new last visited spins in vertex p
   if(p<=mm)then
    if(reA(s1)>0)then
      lastspin(s1)=v0+3 
    else
       lastspin2(s1)=v0+3
    end if 
    if(reA(s2)>0)then
     lastspin(s2)=v0+4 
    else
     lastspin2(s2)=v0+4   
    end if
    if(reA(s3)>0)then
     lastspin(s3)=v0+5 
    else
     lastspin2(s3)=v0+5 
    end if
   elseif(p>mm)then
    if(reA(s1)>0)then
      lastspin(s1)=v0+3
    else
       lastspin3(s1)=v0+3
    end if
    if(reA(s2)>0)then
     lastspin(s2)=v0+4
    else
     lastspin3(s2)=v0+4
    end if
    if(reA(s3)>0)then
     lastspin(s3)=v0+5
    else
     lastspin3(s3)=v0+5
    end if
   
   end if   

  else 
  
   vtx(v0:v0+5)=-1 ! -1 indicates no links in operators.

  end if        

 end do

! links across the "imaginary time boundary" p=1 <--> p=mm

 do i=1,L
  it=frstspin(i)
  io=lastspin(i)
  if(it/=-1)then
   vtx(io)=it
   vtx(it)=io
   anylink=1 
  end if        
  it=frstspin2(i)
  io=lastspin2(i)
  if(it/=-1)then
   vtx(io)=it
   vtx(it)=io
   anylink=1
  end if
  it=frstspin3(i)
  io=lastspin3(i)
  if(it/=-1)then
   vtx(io)=it
   vtx(it)=io
   anylink=1 
  end if 
  

 end do


!  write(6,*)'touching right spins during link list const'
!  write(6,*) 'i, reA(i), frstspin1(i), frstspin2(i),frstspin3(i)'  
!  do i=1,L

!  write(6,*) i, reA(i), frstspin(i), frstspin2(i),frstspin3(i)
!  write(6,*) i, reA(i), lastspin(i), lastspin2(i),lastspin3(i)
!  end do
!  write(6,*) " "


  !time1=mclock()
  !timeli=timeli+time1-time
  !write(6,*)'linked vertex list', timeli,time1-time
  !time=mclock()

!  write(6,*) 'i, vtype,plaquette'
!  do i=1,mm
!   write(6,*)i,vtype(i),sm(i)/2
!  end do 

!  write(6,*) 'vtx table'
!  do i=0,6*mm-1

!    write(6,*)i, vtx(i)
!  end do

! write(6,*) 'proposals'
! do i=1,nl
!  write(6,*) i,tryloop(i,1),tryloop(i,2),tryloop(i,3),tryloop(i,4)
! end do
! write(6,*) 'spin',spin


!----- loop construction and update------------------------
 if(anylink==-1)return
 totvisit=0

 vtypel(1:mm)=vtype
 vtypel(mm+1:mm2+mm)=vtype2
 sml(1:mm)=sm
 sml(mm+1:mm+mm2)=sm2

! write(6,*) '======= start operator loops======'
loop0: do i=1,nl

  b=-1
  do while(b<0)
   br=drand1()*6*(mm+mm2)
   b=vtx(br)
  end do

  jo=br

  j=jo

  break=0
  looplength=0
loop1:  do while(break==0)  ! loop operator construction

  p=j/6+1 !vertex number
    
  b=sml(p)/2 ! plaquette where the operator is located

  li=mod(j,6) !entrance leg index
 

! **** selection of exit leg ******
  nran=drand1()
  le=0
  ex=0


loopleg:  do while(ex==0)

   !write(6,*)'??', p
   if(vtypel(p)<0) write(6,*) 'p,vtypel(p)', p, vtypel(p)


   if(nran<prob(b,vtypel(p),li,le))then
    ex=1
    exit loopleg
   end if 
   le=le+1

  end do loopleg


  ! le is the exit leg

      if(newvtx(li,le,vtypel(p))<0)then
       write(6,*) 'li,le,vtypel(p),nran', li,le,vtypel(p),nran
       
       write(6,*) 'nb,b,vtypel,li',nb,b,vtypel(p),li
       write(6,*) 'probs',prob(b,vtypel(p),li,0),prob(b,vtypel(p),li,1), &
       prob(b,vtypel(p),li,2),prob(b,vtypel(p),li,3)


       stop
      end if


   vtypel(p)=newvtx(li,le,vtypel(p)) ! update vertex   


   j=6*(p-1)+le
   jf=j

   if(j==jo)then ! head finds tail and closes the loop?
    exit loop1
   end if

   j=vtx(j)

   if(j==jo)then
    break=1
   end if   


   if(li/=le)then 
     looplength=looplength+1 
   end if
  
   ! if the loop does not close after a long construction  it should be killed.

   if(looplength>=maxloop)then 
    write(6,*)'broken===========**************************************'
     looplength=-1
    exit loop0
    exit loop1
  end if



  end do loop1
  


    totvisit=totvisit+looplength 

 end do loop0

 !time1=mclock()
 !timelo=timelo+time1-time
 !write(6,*) 'loop construction', timelo, time1-time

 if(looplength>=0)then
  
  vtype=vtypel(1:mm)
  vtype2=vtypel(mm+1:mm2+mm) 
  
  !update of the operator chain after the loops

!****************************************????????????????????????????????????*************************
  do i=1,mm+mm2

   if(i<=mm)then
    if(sm(i)>0)then 
     b=sm(i)/2
     sm(i)=2*b+optype(vtype(i)) 
    end if     
   elseif(i>mm)then
    if(sm2(i-mm)>0)then
     b=sm2(i-mm)/2
     sm2(i-mm)=2*b+optype(vtype2(i-mm))
    end if

   end if   

  end do 
 
  ! update of spin configuration (only spin and the region \bar{A} of spin2)

  do i=1,L
   if(reA(i)>0)then 
    if(frstspin(i)/=-1)then
     p=frstspin(i)/6+1
     le=mod(frstspin(i),6)  
     spin(i)=legspin(le,vtypel(p))
    else ! free spin
     !if(thermalization==0)then       
      if (drand1()<0.5) spin(i)=-spin(i)       
     !end if 
    end if
   else
    if(frstspin2(i)/=-1)then
     p=frstspin2(i)/6+1
     le=mod(frstspin2(i),6)
     spin(i)=legspin(le,vtypel(p))
    else ! free spin
     !if(thermalization==0)then
      if (drand1()<0.5) spin(i)=-spin(i)
     !end if
    end if 
    if(frstspin3(i)/=-1)then
     p=frstspin3(i)/6+1
     le=mod(frstspin3(i),6)
     spin2(i)=legspin(le,vtypel(p))
    else ! free spin
     !if(thermalization==0)then
      if (drand1()<0.5) spin2(i)=-spin2(i)
     !end if
    end if       
   end if
  end do 
 end if        

write(6,*) 'spin end of loop update',spin
write(6,*) 'spin2 end of loop update',spin2
write(6,*) 'after all loops whichreg',whichreg
do i=1,mm
write(6,*)  i,'vtp',vtype(i),'sm',sm(i), 'bond=',sm(i)/2, bsites(:,sm(i)/2)
write(6,*) 
end do

write(6,*)'replica'
do i=1,mm2
write(6,*)  i,'vtp',vtype2(i),'sm',sm2(i), 'bond=',sm2(i)/2, bsites(:,sm2(i)/2)
write(6,*)
end do

!write(6,*)  'states on the legs'
!do i=1,mm
!write(6,*)i, legspin(0:3,vtype(i))
!end do

 !----------------------------------------------------------
 end subroutine loopupdate 
                      
                     
 subroutine bouncefree(ma,i,vlist,ok)
 use configuration
 implicit none
 real(8) ma(6,6),w1,w2,w3,w4,w5,w6
 real(8)dw(6),dwr(6),ap(6,6),C,D
 integer(4)indx(6),ii,jj,i,vlist(6)
 integer :: ok,IER


 ! Syljuasen solution
 ! 
 ! 0 1 1 1 1 1 1
 ! 1 0 1 0 0 0 0
 ! 1 1 0 0 0 0 0 
 ! 1 0 0 0 0 1 0
 ! 1 0 0 0 1 0 1 
 ! 1 0 0 0 0 1 0 

 do ii=1,6
  dw(ii)=vtex(i,vlist(ii))
 end do

dwr=dw

indx=-1

! write(6,*)'bfore', dw

 
 call DPSORT(dw,6,indx,-2,IER) 
!  write(6,*)'after', dw
!  write(6,*) 'index',indx
   
  !DSORTX (x, incx, n, indx)
 ! call dsortx(dw,1,4,indx)
  

 ap=0.0d0

 ap(indx(1),indx(6))=dw(6)/2.0
 ap(indx(1),indx(5))=(dw(5)-dw(6))/2.0
 ap(indx(1),indx(4))=dw(4)-dw(5)/2.0
 ap(indx(4),indx(5))=dw(5)/2.0
 ap(indx(5),indx(6))=dw(6)/2.0
 ap(indx(1),indx(2))=(dw(1)+dw(2)-dw(3)-dw(4))/2.0d0
 ap(indx(1),indx(3))=(dw(1)-dw(2)+dw(3)-dw(4))/2.0d0
 ap(indx(2),indx(3))=(-dw(1)+dw(2)+dw(3)+dw(4))/2.0d0

 
 ap(indx(6),indx(1))=ap(indx(1),indx(6))
 ap(indx(5),indx(1))=ap(indx(1),indx(5))
 ap(indx(4),indx(1))=ap(indx(1),indx(4))
 ap(indx(5),indx(4))=ap(indx(4),indx(5))
 ap(indx(6),indx(5))=ap(indx(5),indx(6))
 ap(indx(2),indx(1))=ap(indx(1),indx(2))
 ap(indx(3),indx(1))=ap(indx(1),indx(3)) 
 ap(indx(3),indx(2))=ap(indx(2),indx(3))
 

 
 ma=ap
  
! write(6,*) dwr(1),sum(ap(1,:))
! write(6,*) dwr(2),sum(ap(2,:))
! write(6,*) dwr(3),sum(ap(3,:))
! write(6,*) dwr(4),sum(ap(4,:))
! write(6,*) dwr(5),sum(ap(5,:))
! write(6,*) dwr(6),sum(ap(6,:))
 

! write(6,*)'MATRIX',ap
 
 
 ok=1 

 do ii=1,6
  do jj=1,6
   if(ma(ii,jj)<0.0d0)ok=0
!    write(6,*) ii,jj,ap(ii,jj)
  enddo
 enddo

! write(6,*) 'a23,a13,a12',  ap(indx(2),indx(3)),ap(indx(1),indx(3)),ap(indx(1),indx(2))


 ! Bounce from the largest vertex only solution
 ! 
 ! 1 1 1 1 1 1 1
 ! 1 0 0 0 0 0 0
 ! 1 0 0 0 0 0 0 
 ! 1 0 0 0 0 0 0
 ! 1 0 0 0 0 0 0 
 ! 1 0 0 0 0 0 0 
 
 if(ok==0)then
  !write(6,*)'a23',ap(indx(2),indx(3))
  ap=0.0d0

  ap(indx(1),indx(1))=dw(1)-(dw(2)+dw(3)+dw(4)+dw(5)+dw(6))

  ap(indx(1),indx(2))=dw(2)
  ap(indx(1),indx(3))=dw(3)
  ap(indx(1),indx(4))=dw(4)
  ap(indx(1),indx(5))=dw(5)
  ap(indx(1),indx(6))=dw(6)
  
  ap(indx(2),indx(1))=dw(2)
  ap(indx(3),indx(1))=dw(3)
  ap(indx(4),indx(1))=dw(4)
  ap(indx(5),indx(1))=dw(5)
  ap(indx(6),indx(1))=dw(6) 

  ma=ap

! write(6,*)' W1 is large'
  
! write(6,*) dwr(1),sum(ap(1,:))

! write(6,*) dwr(2),sum(ap(2,:))

! write(6,*) dwr(3),sum(ap(3,:))

! write(6,*) dwr(4),sum(ap(4,:))

! write(6,*) dwr(5),sum(ap(5,:))

! write(6,*) dwr(6),sum(ap(6,:))

!write(6,*) 'MATRIX W1 is large'
! do ii=1,6
!  do jj=1,6
!    write(6,*) ii,jj,ap(ii,jj)
!  enddo
! enddo

 end if


 ! double check everything is OK, but it should already be. 
 ok=1

 do ii=1,6
  do jj=1,6
   if(ma(ii,jj)<0.0d0)ok=0
  enddo
 enddo
 if(ok==0)then
  write(6,*) 'wrong solutions'
  stop
 end if
 return
 end subroutine bouncefree
 
 subroutine assign(wma,ii,vlist)
 use configuration
 implicit none 
 integer(4)ii,i1,i2,i3,i4,i5,i6,en,ex,vt(0:5),li,vlist(6)
 real(8)tot,wma(6,6)  

  do  li=1,6
   vt(li-1)=vlist(li) 
  end do
 
  do en=0,5
   
   tot=sum(wma(en+1,:)) 
   ex=0
   prob(ii,vt(en),en,ex)=wma(en+1,ex+1)/tot 
   !write(6,*)ii,vt(en),en,ex,wma(en+1,ex+1)/tot 
   do ex=1,5     
     prob(ii,vt(en),en,ex)=prob(ii,vt(en),en,ex-1)+wma(en+1,ex+1)/tot
    ! write(6,*)ii,vt(en),en,ex,wma(en+1,ex+1)/tot
   enddo 
  end do

 end subroutine assign
subroutine probupdates !_loopupdate_XXZ
 
 use configuration
 implicit none
 integer :: i,j,k,ok,update,el,vlist(6)
 real(kind=8)ma(6,6),mach_prec,tot,smallest
 

 allocate(vtex(nb,nver),prob(nb,nver,0:5,0:5))

 smallest=10000000000.0
! prob(bond, vertex type( bond p), entrance leg, exit leg)

 do i=1,nb
   

  ! calculation of the vertices 

  ! diagonal vertices

  do j=1,8
   vtex(i,j)=cadd+ 0.5d0/z*(epsdis(bsites(1,i))*legspin(0,j)+ &
             epsdis(bsites(2,i))*legspin(1,j)+ epsdis(bsites(3,i))*legspin(2,j)     )-&
             0.25*Jz*( legspin(0,j)*legspin(1,j)+legspin(0,j)*legspin(2,j) + legspin(1,j)*legspin(2,j) )
  end do

  ! hopping vertices
  do j=9,20
   vtex(i,j)=t/2.0
  end do

  ! S+S+ S-S- vertices
  do j=21,32
   vtex(i,j)=tp/2.0
  end do

  do j=1,nver
 
   if(j<9)then 
    if(vtex(i,j)<smallest)then
     smallest=vtex(i,j) 
    end if        
   end if
   if(vtex(i,j)<0)then
    write(6,*) 'negative vertex i, j , V(i,j)', i,j, vtex(i,j)       
    stop 
   end if        

  end do

  
  ! Solution of directed loop equations
  do j=1,nver
   
   ! finds out the closed sets of vertices
   do el=0,5
    vlist(el+1)=newvtx(0,el,j)     
   end do  
   
   ! determine the solutions 
   call bouncefree(ma,i,vlist,ok)
   
   ! assign the solutions to the probabilities
   call assign(ma,i,vlist)
 
  end do

 
 end do

 write(6,*) '+++++++++++++++ Smallest diagonal VTEX=====',smallest
 
 end subroutine probupdates !_loopupdate_XXZ


! ----------- probabilities according to the loop update (not directed loop
! update)
 subroutine probupdates_heat_bath
 
 use configuration
 implicit none
 integer :: i,j,k,m
 real(8) deno
 allocate(vtex(nb,nver),prob(nb,nver,0:5,0:5))

 prob=0.0d0
! prob(bond, vertex type( bond p), entrance leg, exit leg)

 do i=1,nb

  
 ! diagonal vertices
  do j=1,8
   vtex(i,j)=cadd+ 0.5d0/z*(epsdis(bsites(1,i))*legspin(0,j)+ &
             epsdis(bsites(2,i))*legspin(1,j)+ epsdis(bsites(3,i))*legspin(2,j)     )-&
             0.25*Jz*( legspin(0,j)*legspin(1,j)+legspin(0,j)*legspin(2,j) + legspin(1,j)*legspin(2,j) )
  
  end do

  ! hopping vertices
  do j=9,20
   vtex(i,j)=t/2.0  
  end do

  ! S+S+ S-S- vertices
  do j=21,32
   vtex(i,j)=tp/2.0
  end do 

 end do
 

 do i=1,nb
     ! generation of the probabilities
 
  do j=1,nver
   do k=0,5 ! loop entrance leg 

    !write(6,*)'evertex,entrance k, exit m, overtex' 
    deno=0.0d0
    do m=0,5 ! loop exit leg
     !newvtx(en,ex,i)
     deno=deno+vtex(i,newvtx(k,m,j))
     !write(6,*)'evertex,entrance k, exit m, overtex'
    ! write(6,*)j,k,m,newvtx(k,m,j)
    end do

    prob(i,j,k,0)=vtex(i,newvtx(k,0,j))/deno
    do m=1,5
     prob(i,j,k,m)=prob(i,j,k,m-1)+vtex(i,newvtx(k,m,j))/deno
    end do
    
   end do

  end do

 end do

 end subroutine probupdates_heat_bath



 subroutine newvertex
 
 use configuration
 implicit none
 integer i,j,k,li,le,vertex,en,ex
 integer,allocatable :: legt(:)

 nver=32
 allocate(newvtx(0:5,0:5,nver),optype(nver),legspin(0:5,nver),legt(0:5))


 ! diagonal vertices

 !vertex type 1    0 0 0             3 4 5   (names of the legs)
 !                 =====             =====
 !                 0 0 0             0 1 2 

 optype(1)=0 

 legspin(0,1)=-1
 legspin(1,1)=-1
 legspin(2,1)=-1
 legspin(3,1)=-1
 legspin(4,1)=-1
 legspin(5,1)=-1
 


!vertex type 2    X 0 0
!                 =====
!                 X 0 0

 optype(2)=0

 legspin(0,2)=1
 legspin(1,2)=-1
 legspin(2,2)=-1
 legspin(3,2)=1
 legspin(4,2)=-1
 legspin(5,2)=-1


!vertex type 3   0 X 0 
!                =====
!                0 X 0 

 optype(3)=0

 legspin(0,3)=-1
 legspin(1,3)=1
 legspin(2,3)=-1
 legspin(3,3)=-1
 legspin(4,3)=1
 legspin(5,3)=-1


!vertex type 4   0 0 X 
!                =====
!                0 0 X 

 optype(4)=0

 legspin(0,4)=-1
 legspin(1,4)=-1
 legspin(2,4)=1
 legspin(3,4)=-1
 legspin(4,4)=-1
 legspin(5,4)=1
 

!vertex type 5   X X 0 
!                =====
!                X X 0 

 optype(5)=0

 legspin(0,5)=1
 legspin(1,5)=1
 legspin(2,5)=-1
 legspin(3,5)=1
 legspin(4,5)=1
 legspin(5,5)=-1


!vertex type 6   X 0 X 
!                =====
!                X 0 X 

 optype(6)=0

 legspin(0,6)=1
 legspin(1,6)=-1
 legspin(2,6)=1
 legspin(3,6)=1
 legspin(4,6)=-1
 legspin(5,6)=1


!vertex type 7   0 X X 
!                =====
!                0 X X 

 optype(7)=0

 legspin(0,7)=-1
 legspin(1,7)=1
 legspin(2,7)=1
 legspin(3,7)=-1
 legspin(4,7)=1
 legspin(5,7)=1


!vertex type 8   X X X 
!                =====
!                X X X 

 optype(8)=0

 legspin(0,8)=1
 legspin(1,8)=1
 legspin(2,8)=1
 legspin(3,8)=1
 legspin(4,8)=1
 legspin(5,8)=1



! S+S- S-S+ vertices


!vertex type 9   0 X 0  
!                =====
!                X 0 0 

 optype(9)=1

 legspin(0,9)=1
 legspin(1,9)=-1 
 legspin(2,9)=-1
 legspin(3,9)=-1
 legspin(4,9)=1
 legspin(5,9)=-1 


!vertex type 10  X 0 0  
!                =====
!                0 X 0 
      
 optype(10)=1
 
 legspin(0,10)=-1
 legspin(1,10)=1
 legspin(2,10)=-1
 legspin(3,10)=1
 legspin(4,10)=-1
 legspin(5,10)=-1



!vertex type 11  0 0 X  
!                =====
!                0 X 0 
      
 optype(11)=1
 
 legspin(0,11)=-1
 legspin(1,11)=1
 legspin(2,11)=-1
 legspin(3,11)=-1
 legspin(4,11)=-1
 legspin(5,11)=1


!vertex type 12  0 X 0  
!                =====
!                0 0 X 

 optype(12)=1
 
 legspin(0,12)=-1
 legspin(1,12)=-1
 legspin(2,12)=1
 legspin(3,12)=-1
 legspin(4,12)=1
 legspin(5,12)=-1


!vertex type 13  0 0 X  
!                =====
!                X 0 0 

 optype(13)=1

 legspin(0,13)=1
 legspin(1,13)=-1
 legspin(2,13)=-1
 legspin(3,13)=-1
 legspin(4,13)=-1
 legspin(5,13)=1


!vertex type 14  X 0 0  
!                =====
!                0 0 X 

 optype(14)=1

 legspin(0,14)=-1
 legspin(1,14)=-1
 legspin(2,14)=1
 legspin(3,14)=1
 legspin(4,14)=-1
 legspin(5,14)=-1


!vertex type 15  0 X X  
!                =====
!                X 0 X 

 optype(15)=1

 legspin(0,15)=1
 legspin(1,15)=-1
 legspin(2,15)=1
 legspin(3,15)=-1
 legspin(4,15)=1
 legspin(5,15)=1


!vertex type 16  X 0 X  
!                =====
!                0 X X 

 optype(16)=1

 legspin(0,16)=-1
 legspin(1,16)=1
 legspin(2,16)=1
 legspin(3,16)=1
 legspin(4,16)=-1
 legspin(5,16)=1


!vertex type 17  X 0 X  
!                =====
!                X X 0 

 optype(17)=1

 legspin(0,17)=1
 legspin(1,17)=1
 legspin(2,17)=-1
 legspin(3,17)=1
 legspin(4,17)=-1
 legspin(5,17)=1


!vertex type 18  X X 0  
!                =====
!                X 0 X 

 optype(18)=1

 legspin(0,18)=1
 legspin(1,18)=-1
 legspin(2,18)=1
 legspin(3,18)=1
 legspin(4,18)=1
 legspin(5,18)=-1


!vertex type 19  X X 0  
!                =====
!                0 X X 

 optype(19)=1

 legspin(0,19)=-1
 legspin(1,19)=1
 legspin(2,19)=1
 legspin(3,19)=1
 legspin(4,19)=1
 legspin(5,19)=-1


!vertex type 20  0 X X  
!                =====
!                X X 0 

 optype(20)=1

 legspin(0,20)=1
 legspin(1,20)=1
 legspin(2,20)=-1
 legspin(3,20)=-1
 legspin(4,20)=1
 legspin(5,20)=1


! S+S+ S-S- vertices

!vertex type 21  0 0 0  
!                =====
!                X X 0 

 optype(21)=1

 legspin(0,21)=1
 legspin(1,21)=1
 legspin(2,21)=-1
 legspin(3,21)=-1
 legspin(4,21)=-1
 legspin(5,21)=-1

 
!vertex type 22  0 0 0  
!                =====
!                X 0 X 

 optype(22)=1

 legspin(0,22)=1
 legspin(1,22)=-1
 legspin(2,22)=1
 legspin(3,22)=-1
 legspin(4,22)=-1
 legspin(5,22)=-1


!vertex type 23  0 0 0  
!                =====
!                0 X X 

 optype(23)=1

 legspin(0,23)=-1
 legspin(1,23)=1
 legspin(2,23)=1
 legspin(3,23)=-1
 legspin(4,23)=-1
 legspin(5,23)=-1



!vertex type 24  0 0 X  
!                =====
!                X X X 

 optype(24)=1

 legspin(0,24)=1
 legspin(1,24)=1
 legspin(2,24)=1
 legspin(3,24)=-1
 legspin(4,24)=-1
 legspin(5,24)=1


!vertex type 25  0 X 0  
!                =====
!                X X X 

 optype(25)=1

 legspin(0,25)=1
 legspin(1,25)=1
 legspin(2,25)=1
 legspin(3,25)=-1
 legspin(4,25)=1
 legspin(5,25)=-1


!vertex type 26  X 0 0  
!                =====
!                X X X 

 optype(26)=1

 legspin(0,26)=1
 legspin(1,26)=1
 legspin(2,26)=1
 legspin(3,26)=1
 legspin(4,26)=-1
 legspin(5,26)=-1



!vertex type 27  X X X  
!                =====
!                0 0 X 

 optype(27)=1

 legspin(0,27)=-1
 legspin(1,27)=-1
 legspin(2,27)=1
 legspin(3,27)=1
 legspin(4,27)=1
 legspin(5,27)=1


!vertex type 28  X X X  
!                =====
!                0 X 0 

 optype(28)=1

 legspin(0,28)=-1
 legspin(1,28)=1
 legspin(2,28)=-1
 legspin(3,28)=1
 legspin(4,28)=1
 legspin(5,28)=1


!vertex type 29  X X X  
!                =====
!                X 0 0 

 optype(29)=1

 legspin(0,29)=1
 legspin(1,29)=-1
 legspin(2,29)=-1
 legspin(3,29)=1
 legspin(4,29)=1
 legspin(5,29)=1


!vertex type 30  X X 0  
!                =====
!                0 0 0 

 optype(30)=1

 legspin(0,30)=-1
 legspin(1,30)=-1
 legspin(2,30)=-1
 legspin(3,30)=1
 legspin(4,30)=1
 legspin(5,30)=-1



!vertex type 31  X 0 X  
!                =====
!                0 0 0 

 optype(31)=1

 legspin(0,31)=-1
 legspin(1,31)=-1
 legspin(2,31)=-1
 legspin(3,31)=1
 legspin(4,31)=-1
 legspin(5,31)=1


!vertex type 32  0 X X  
!                =====
!                0 0 0 

 optype(32)=1

 legspin(0,32)=-1
 legspin(1,32)=-1
 legspin(2,32)=-1
 legspin(3,32)=-1
 legspin(4,32)=1
 legspin(5,32)=1


!write(6,*) 'newvertex table'
!write(6,*) 'entrance leg, exit leg, input vertex, output vertex '
do i=1,nver
 do en=0,5
   do ex=0,5
    legt=legspin(:,i)
    legt(en)=-legt(en)
    legt(ex)=-legt(ex)
    call searchvtx(vertex,legt)
    newvtx(en,ex,i)=vertex
    !write(6,*) en,ex,i,vertex
   end do 
 end do 
end do

end subroutine newvertex

subroutine searchvtx(vertex,legt)
use configuration
implicit none
integer(4) vertex, i ,j,k,legt(0:5),sumi

search: do i=1,nver

 sumi=0
 do j=0,5
  sumi=sumi+abs(legt(j)-legspin(j,i))
 end do
 if(sumi==0)then
  vertex=i 
  exit search
 end if
end do search

end subroutine searchvtx


!-----------------------------!
 subroutine adjustcutoff(step)
!-----------------------------!
 use configuration
 implicit none

 integer, allocatable :: stringcopy(:),vtypecopy(:),stringcopy2(:),vtypecopy2(:),vtypelcopy(:)
 integer :: mmnew,mmnew2,step

! mmnew=nh+nh/3 !serial version

  mmnew=nh0+nh0/3
  mmnew2=nh20+nh20/3

 if(mmnew>=mm)then 
  allocate(stringcopy(mm))
  allocate(vtypecopy(mm))
  stringcopy(:)=sm(:)
  vtypecopy(:)=vtype(:)
  deallocate(sm,vtype)
  allocate(sm(mmnew),vtype(mmnew))
  sm(1:mm)=stringcopy(:)
  sm(mm+1:mmnew)=0
  vtype(1:mm)=vtypecopy(:)
  vtype(mm+1:mmnew)=0
  deallocate(stringcopy,vtypecopy)
  mm=mmnew
 endif
 if(mmnew2>=mm2)then
  allocate(stringcopy(mm2))
  allocate(vtypecopy(mm2))
  stringcopy(:)=sm2(:)
  vtypecopy(:)=vtype2(:)
  deallocate(sm2,vtype2)
  allocate(sm2(mmnew2),vtype2(mmnew2))
  sm2(1:mm2)=stringcopy(:)
  sm2(mm2+1:mmnew2)=0
  vtype2(1:mm2)=vtypecopy(:)
  vtype2(mm2+1:mmnew2)=0
  deallocate(stringcopy,vtypecopy)
  mm2=mmnew2
 endif

 deallocate(vtypel,vtx,sml)
 allocate(vtypel(mm+mm2),vtx(0:6*(mm+mm2)-1),sml(mm+mm2))


 end subroutine adjustcutoff

!------------------------------------!
 subroutine writeresults(msteps,bins,size,rank,nbins,iseedd)
!------------------------------------!
 use configuration;use mpi; use measurementdata; implicit none

 integer :: i,msteps,bins,shift,size,rank,ierr,j,nbins,qv,k,iseedd
 integer(4) stats(MPI_STATUS_SIZE)

 real(8) :: wdata1(8),wdata2(8), nk0,nk01,errnk0, dii,value ! ,wdatgf(L,L),wdat2gf(L,L)
 complex(8)valuec

 enrg1=enrg1/msteps
 enrg2=enrg2/msteps

 enrg2=(enrg2-enrg1*(enrg1+1.d0))/dble(L)
 enrg1=-enrg1/(2.0*beta*dble(L))+cadd*dble(nb)/dble(L)
 !write(6,*) 'sqmsteps',sqmsteps 
 
 
 token(1)=enrg1
 token(2)=enrg2

 call mpi_reduce(token,token0,8,mpi_real8,mpi_sum,0,mpi_comm_world,ierr)

 if(rank==0)then ! rank 0

  token0=token0/dble(size)
  ! accumulate the measurements
  data1(1)=data1(1)+token0(1)
  data1(2)=data1(2)+token0(2)
  data1(3)=data1(3)+token0(3)
  data1(4)=data1(4)+token0(4)
  data1(6)=data1(6)+token0(6)
  data1(5)=data1(5)+token0(5)
  data1(7)=data1(7)+token0(7) 
  data1(8)=data1(8)+token0(8)
  
  
  ! and the square of them
  data2(1)=data2(1)+token0(1)**2
  data2(2)=data2(2)+token0(2)**2
  data2(3)=data2(3)+token0(3)**2
  data2(4)=data2(4)+token0(4)**2
  data2(6)=data2(6)+token0(6)**2
  data2(5)=data2(5)+token0(5)**2
  data2(7)=data2(7)+token0(7)**2
  data2(8)=data2(8)+token0(8)**2
  

  do i=1,8
    wdata1(i)=data1(i)/bins
    wdata2(i)=data2(i)/bins
    wdata2(i)=sqrt(abs(wdata2(i)-wdata1(i)**2)/bins)
  end do



!************************************************************
  open(10,file='results.txt',status='replace')
  write(10,*)' Cut-off L : ',mm
  write(10,*)' Number of bins completed : ',bins
  write(10,*)' ========================================='
  write(10,10)' -E/N       : ',wdata1(1),wdata2(1)
  write(10,*)' ========================================='
  10 format(1x,a,2f14.8)
  close(10)

open(unit=55,file='data.dat',form='unformatted',status='unknown')
rewind(55)
write(55)bins,data1,data2,lx,ly,beta,mu,t,tp,Jz,iseedd

close(55)
 
token=0.0d0
token0=0.0d0


end if ! rank 0 

 enrg1=0.d0
 enrg2=0.d0
  
 end subroutine writeresults
!---------------------------!

subroutine check(mc)
use configuration
use measurementdata
implicit none


integer :: i,j,b,op,s1,s2,am,jj(0:2),opt,kkk,kk,nt,kin,btype,i1,i2,a1,a2,ax1,ax2,ay1,ay2,counter,mc,modi,mea
real(8)current,compi,comp2i

 enrg1=enrg1+dfloat(nh+nh2)
 enrg2=enrg2+dfloat(nh+nh2)**2
 
end subroutine check

subroutine coord(i,a1,a2)
use configuration
use measurementdata
implicit none
integer(4)s1,a1,a2,k,kk,i


if(mod(i,3*ly)==0)then
  k=3*ly
 else
  k=mod(i,3*ly)
 end if

 !write(*,*)i, ceiling(dble(k)/dble(3))-1,mod(i,3) !! with ceiling... and mod i,3 one can get the a2 coordinate

 if(mod(i,3)==0)then
  kk=3
 else
  kk=mod(i,3)
 end if


! write(*,*)i, ceiling(dble(k)/dble(3))-1, ceiling(dble(i)/dble(3*ly))-1,kk

 if(kk==1)then
  a2=2*(ceiling(dble(k)/dble(3))-1)
  a1=2*(ceiling(dble(i)/dble(3*ly))-1)
 elseif(kk==2)then
  a2=2*(ceiling(dble(k)/dble(3))-1)
  a1=2*(ceiling(dble(i)/dble(3*ly))-1) +1
 elseif(kk==3)then
  a2=2*(ceiling(dble(k)/dble(3))-1)+1
  a1=2*(ceiling(dble(i)/dble(3*ly))-1)
 end if

end subroutine coord

subroutine cdist()
use configuration
use measurementdata
implicit none
integer(4) i,j,k,d,ci,cit(1),cnt,modi,found
real(8) dist9(9),aa1,aa2,candidate(2),diff,eps


eps=1.0d-12

allocate(distances(nbase*L,dd),cdistances(nbase*L))

! kagome 
ba1(1)=2.0d0*dble(lx)
ba1(2)=0.0d0
ba2(1)=2.0d0*dble(ly)*dcos(pi/3.0d0)
ba2(2)=2.0d0*dble(ly)*dsin(pi/3.0d0)


d=1

do i=1,nbase !,L
 do j=1,L ! ,L
   
   aa1=0.0d0
   aa2=0.0d0 
   dist9(1)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ) )**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=1.0d0
   aa2=0.0d0 
   dist9(2)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  ) 
  
   aa1=1.0d0
   aa2=1.0d0 
   dist9(3)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=0.0d0
   aa2=1.0d0  
   dist9(4)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=-1.0d0
   aa2=1.0d0
   dist9(5)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=-1.0d0
   aa2=0.0d0
   dist9(6)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
  
   aa1=-1.0d0
   aa2=-1.0d0
   dist9(7)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
   
   aa1=0.0d0
   aa2=-1.0d0
   dist9(8)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )

   aa1=1.0d0
   aa2=-1.0d0
   dist9(9)=sqrt( (ord(i,1)-( ord(j,1) +aa1*ba1(1)+aa2*ba2(1) ))**2+(ord(i,2)-( ord(j,2)+aa1*ba1(2)+aa2*ba2(2)  ))**2  )
   
   
   cit=minloc(dist9)
   ci=cit(1)
   

   if(ci==1)then
    aa1=0.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)1,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==2)then
    aa1=1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)2,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==3)then
    aa1=1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)3,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==4)then
    aa1=0.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)4,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==5)then
    aa1=-1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)5,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==6)then
    aa1=-1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)6,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==7)then
    aa1=-1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)7,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==8)then
    aa1=0.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)8,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   elseif(ci==9)then
    aa1=1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !do k=1,ndim
    !write(6,*)9,exp(2*pi*ima*sum(candidate(:)*q(k,:)))
    !enddo
   endif 
   !write(6,*) '=================================' 
   !write(6,*)i,j,ci,dist9  
   !write(6,*) 'distance',candidate 
   !write(6,*) '================================'
    
   
    distances(d,:)=candidate
    d=d+1 
   
 end do
end do



indist=nbase*L
!write(6,*)'indist',indist,d-1 
!do i=1,nbase*L
! write(6,*)i,distances(i,1), distances(i,2),sqrt(distances(i,1)**2+ distances(i,2)**2)
!end do

allocate(site_to_d(L,L))

cdistances=0

do i=1,L
 do j=1,L

   modi=mod(i,nbase)
   if(modi==0)modi=nbase
   found=0 
   loopkk: do k=1,L

    aa1=0.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
    !write(6,*)'00', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk 
    end if  

    aa1=1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'10', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if
  
    aa1=1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'11', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if

    aa1=0.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'01', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if

    aa1=-1.0d0
    aa2=1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'-11', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:) 
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if

    aa1=-1.0d0
    aa2=0.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'-10', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k) +1
     found=1
     exit loopkk
    end if     
 
    aa1=-1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))

     !write(6,*)'-1-1', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k) +1
     found=1
     exit loopkk
    end if

    aa1=0.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'0-1', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)= cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if
    aa1=1.0d0
    aa2=-1.0d0
    candidate(1)=-(ord(i,1)-(ord(j,1)+aa1*ba1(1)+aa2*ba2(1)))
    candidate(2)=-(ord(i,2)-(ord(j,2)+aa1*ba1(2)+aa2*ba2(2)))
     !write(6,*)'11', k,modi, (modi-1)*L+k, candidate,-distances((modi-1)*L+k,:)
    diff=(candidate(1)-distances((modi-1)*L+k,1) )**2+(candidate(2)-distances((modi-1)*L+k,2) )**2
    if(diff<eps)then
     site_to_d(i,j)=(modi-1)*L+k
     cdistances((modi-1)*L+k)=cdistances((modi-1)*L+k)+1
     found=1
     exit loopkk
    end if
 
   end do loopkk
  ! write(6,*) 'found', i,j,found
 
 end do
end do  


!write(6,*) 'distances'
!do i=1,L
! do j=1,L
!  write(6,*)i,j,site_to_d(i,j),distances(site_to_d(i,j),:)
! end do
!enddo


!write(6,*)'cdistances'
!do i=1,nbase*L
!write(6,*) i, cdistances(i),distances(i,:)
!end do

end subroutine cdist
subroutine restarts(iseedd)
use configuration
use measurementdata
implicit none
integer lxt,lyt,iseedd,iseeddr
real(8)betat,mut,tt,tpt,Jzt,eps

eps=0.0000000001
if(restart==0)then
  kstart=1
 elseif(restart>0)then

  open(unit=55,file='data.dat',form='unformatted',status='old')
  rewind(55) 
  read(55)kstart,data1,data2,lxt,lyt,betat,mut,tt,tpt,Jzt,iseeddr
         kstart=kstart+1
  if(lxt.ne.lx)then
   write(6,*) 'check restart parameters lx',lxt,lx
   stop
  end if
  if(lyt.ne.ly)then
   write(6,*) 'check restart parameters ly',lyt,ly
   stop
  end if
  if(abs(betat-beta).ge.eps)then
   write(6,*) 'check restart parameters beta',betat,beta
   stop
  end if
  if(abs(mut-mu).ge.eps)then
   write(6,*) 'check restart parameters mu',mut,mu
   stop
  end if
  if(abs(tt-t).ge.eps)then
   write(6,*) 'check restart parameters t',tt,t
   stop
  end if
  if(abs(tpt-tp).ge.eps)then
   write(6,*) 'check restart parameters tp',tpt,tp
   stop
  end if
  if(abs(Jzt-Jz).ge.eps)then
   write(6,*) 'check restart parameters Jz',Jzt,Jz
   stop
  end if
  if(iseeddr==iseedd)then
   write(6,*) 'when restarting you need a different seed for the QMC',iseedd,iseeddr
   write(6,*) 'code selects one for you'
   CALL SYSTEM('echo $RANDOM>random.dat')
   open(unit=100,file='random.dat',form='formatted',status='unknown')
   read(100,*) iseedd
   close(100)
   CALL SYSTEM('rm random.dat')
  endif
 close(55)
 end if
end subroutine restarts


 subroutine getspin2()
 use configuration
 implicit none
 integer :: i, b,op,cc
 real(8) :: drand1, ratio,nran,error
 integer :: spint(3), typo,legt(0:5),vertex
 tspin=spin

 do i=1,mm
  op=sm(i)
  if(mod(op,2)==1)then ! 0ff-diagonal operator found
   b=op/2 ! (its integer part) 
   spin(bsites(1,b))=legspin(3,vtype(i))
   spin(bsites(2,b))=legspin(4,vtype(i))
   spin(bsites(3,b))=legspin(5,vtype(i))
  end if 
 end do


spin2AFTloop=spin
do i=1,L
 if(reA(i)==0)then
  spin2AFTloop(i)=spin2(i)
 end if
end do

tspin2=spin
spin=tspin
tspin=tspin2
spin2=spin2AFTloop
 end subroutine getspin2
 !------------------------------------------------------------

 subroutine checkspins()
 use configuration
 implicit none
 integer(4)i,j,k
 ! whichreg=1 reA=reAi. whichreg=0 reA=reAipl

 
  accepted=.true. 
  do i=1,sizedl 
   if(spin(dl(i))/=spin2AFTloop(dl(i)))then
    accepted=.false.
   end if 
  end do

!write(6,*)'configuration=================================================='
! do i=1,mm
!write(6,*)  i,'vtp',vtype(i),'sm',sm(i), 'bond=',sm(i)/2, bsites(:,sm(i)/2)
!write(6,*)
!end do

!write(6,*)'replica'
!do i=1,mm2
!write(6,*)  i,'vtp',vtype2(i),'sm',sm2(i), 'bond=',sm2(i)/2, bsites(:,sm2(i)/2)
!write(6,*)
!end do
  
! write(6,*)'inside check spins accepted??',accepted, 'whichreg=', whichreg
! write(6,*)'        spin',spin
! write(6,*)'spin2AFTloop',spin2AFTloop
 
 
 end subroutine checkspins
