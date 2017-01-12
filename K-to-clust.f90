program K-to-clust

! program to generate mean concentrations of clusters for a finite-N system
! at a given total concentration of monomers, given equilibrium association 
! constants for clusters, using explicit generation of integer partitions.
!  J. T. Kindt, Emory University, 2012

! Please cite "Accounting for Finite-Number Effects on Cluster Size Distributions
! in Simulations of Equilibrium Aggregation", J. T. Kindt, J. Chem. Theor.
! Comp., volume 9 pp 147-152 (2013)

implicit none

integer, parameter:: Nmax=170
       ! to go beyond Nmax=170, will need to 
       ! do something about overflow of factorials

integer:: comp(0:Nmax+1,3)  ! comps with linked list attached
!         comp(i,1) contains the number of i-mers in the composition/partition
!         comp(i,2) contains the next-highest value of j with comp(j,1)>0
!         comp(i,3) contains the next-lowest value of j with comp(j,1) >0
!         comp(Nsim+1,3) = the largest aggregate present in the composition
!         comp(0,2) = the smallest aggregate (incl. free monomers) present
!         comp(i,3) =0 if i is the smallest aggregate
!         comp(i,2) =Nsim+1 if i is the greatest aggregate
integer:: most_prob_comp(0:Nmax+1)
real*8:: factorial(0:Nmax)
real*8:: facinv(0:Nmax)
real*8::  conc(3*Nmax),q(0:Nmax),Qtot, P ,Psave(0:Nmax+1)
real*8:: bulk_limit_conc(Nmax)
real*8:: aveNout(0:Nmax),totconc,V,rhotot  
real*8:: mserror,sdstep,steepmag
real*8:: Kapp(Nmax),Kactual(Nmax),testC1,testCtot,error,dCdC1
real*8:: Pcalc,biggest_P

integer:: ilist(0:Nmax+1),revlist(0:Nmax+1)
integer:: i,ncomp,j,k,l,nsum,iter,niter,nci_gt_0,ni_in,Nsim

integer::  inext,irest,ibig,nsplit
integer:: ipoint,iloop,inext_minus_1,irestnext,inextpoint,prestart
logical:: ci_gt_0(Nmax)

100 format (20d8.1)

! initialize factorial array
! note: floating-point error for Nmax > 170
factorial(0)=1.0
facinv(0)=1.0d0
do i=1,Nmax
 factorial(i)=factorial(i-1)*i
 facinv(i)=1.0d0/factorial(i)
 !write (*,*) i,factorial(i),facinv(i)
enddo 


nci_gt_0=0
ci_gt_0=.FALSE.
ilist=0
conc=0.0d0
read (*,*) ni_in,Nsim,rhotot
                              ! ni_in: the total number of i-mers whose bulk equilibrium concentration
                              !        will be given
                              ! Nsim: total number of monomers for the simulation to be modeled.
                              ! rhotot: desired total # of monomers/volume.  If rhotot
                              !         set to a number <= 0, then the total
                              !         concentration will be the same as the
                              !         total concentration of the reference
                              !         system. 
                              
                              
if (Nsim.gt.Nmax) then
       write (*,*) "Nsim is greater than Nmax- adjust the Nmax parameter to equal Nsim."       
       stop
endif
if (ni_in.gt.3*Nmax) then
        write (*,*) "# Warning: Will only use concentration values up to ",3*Nmax
        ni_in=3*Nmax
endif
read (*,*) i,conc(i)     ! these can be equilibrium association constants if conc(1)=1.
if (conc(1).le.0.0d0) then
        write (*,*) "Concentration of free monomer must be > 0.0d0."
        stop
endif
totconc=conc(1)
ci_gt_0(1)=.TRUE.
nci_gt_0=1


!    ilist contains the list of all i-mers present; ilist(nci_gt_0) is the biggest aggregate observed.
ilist(0)=0
ilist(1)=1
Kactual=0.0d0
Kactual(1)=1.0d0
write (*,*) "# i     K_assoc(i)"
write (*,*) "#", 1,     Kactual(1)
do j=2,ni_in
 read (*,*) i,conc(i)          ! input bulk system concentration of i-mers
  if (conc(i).gt.0.0d0) then ! keep track of what aggregate sizes are present
       if (i.le.Nsim) then
         ci_gt_0(i)=.TRUE.
         nci_gt_0=nci_gt_0+1
         ilist(nci_gt_0)=i
       endif
  endif
  totconc=totconc+i*conc(i)
  Kactual(i)=conc(i)/conc(1)**i  ! equilibrium constant in bulk
  write (*,*) "#",i,Kactual(i)
enddo 
ilist(nci_gt_0+1)=Nsim+1

! make revlist: keeps track of index of equal or next-lowest aggregate size present;
        ! for any aggregate size i, equal or next-lowest aggregate size is ilist(revlist(i))
revlist=nci_gt_0
do j=0,Nmax+1
  do i=nci_gt_0, 0, -1
    if (ilist(i).gt.j) then
      revlist(j)=i-1
    endif
  enddo
  !write (*,*) j,revlist(j)
enddo


write (*,*) "##########"
!conc_scale_factor=1.0d0/conc(1)
!conc=conc*conc_scale_factor


if (rhotot.le.0.0d0) then
    V=dfloat(Nsim)/totconc   ! calculate volume of Nsim,V,T system at same total
                             ! monomer conc
    write (*,*) "# total monomer concentration is:",totconc
    write (*,*) "# volume is ",Nsim,"/",totconc," = ",V
    rhotot=totconc
else
    V=dfloat(Nsim)/rhotot
    write (*,*) "# total monomer concentration is:",rhotot
    write (*,*) "# volume is ",Nsim,"/",rhotot," = ",V
endif

q(0)=1.0d0  !  no physical meaning; will not affect results; just needs to be defined
            !  as a placeholder in calculating probabilities for case where irest=0 
q(1)=1.0d0
! find single-aggregate partition functions relative to monomer
write (*,*) "# i     q(i)/q(1)"
do i=1,Nsim
  q(i)=Kactual(i)*V**(1-i)
  write (*,*) "# ",i,q(i)
enddo


! generate all possible compositions of a system with Nsim monomers
! start off with composition containing maximum number of biggest defined clusters

ibig=ilist(nci_gt_0)
inext=ibig

aveNout=0.0d0
comp=0
comp(inext,1) = Nsim/inext
comp(inext,2) = Nsim+1  ! linked list pointing upwards; 
                           ! index corresponds to ilist
comp(Nsim+1,3)=inext
irest=mod(Nsim,inext)
write (*,*) "# starting with highest monomer = ",ibig
do while (irest.ne.0)
   irestnext=ilist(revlist(irest))  ! skip down in case irest-mers are not present
   comp(irestnext,1)=irest/irestnext
   comp(irestnext,2)=inext
   comp(inext,3)=irestnext
   !write (*,*) inext,irest,irestnext
   irest=mod(irest,irestnext)
   inext=irestnext
enddo   

comp(0,2)=inext
comp(inext,3)=0
inext=max( comp(1,2),comp(0,2))  

! count probability weight for initial composition
Qtot=0.0d0
Psave=1.0d0

comp(0,1)=0

prestart=Nsim+1  

ncomp=0
biggest_P=0.0d0
do while (inext.le.Nsim)
  ! calculate new probability
     ! weight contributions from i-mers with i >= prestart are
     ! unchanged, save time by not recalculating them
   P=Pcalc(Psave,facinv,q,comp,Nsim,Nmax,prestart)
   if (P.gt.biggest_P) then
        biggest_P=P
        most_prob_comp=comp(:,1)
   endif
     !write (*,*) "# back from Pcalc ",P,ncomp,Qtot

  ! count current composition
  !  
  ncomp=ncomp+1
  Qtot=Qtot+P
  ipoint=comp(Nsim+1,3)
  do while (ipoint.gt.0)
     aveNout(ipoint)=aveNout(ipoint)+comp(ipoint,1)*P
     ipoint=comp(ipoint,3)
  enddo
  
     
  !
  ! determine lowest aggregate size greater than a monomer
    
    inext=max( comp(1,2),comp(0,2))  ! 
   !write (*,*) "biggest and smallest:", comp(Nsim+1,3),inext,comp(0,2)
   !write (*,*)  comp(1:10,1)
  !      explanation: either 0 points to 1 and 1 points up to the next-highest size,
  !                   or there are no free monomers and 0 points directly to inext
  ! divide up all the monomers belonging to one inext-mer, with free monomers added on, ! into the greatest possible aggregate size (i.e., ilist(revlist(inextpoint-1))- mers )

     prestart=comp(inext,2)
     inextpoint=revlist(inext)
     inext_minus_1=ilist(inextpoint-1)  ! note that this will be inext-1 or, if conc(inext-1)=0, 
                                          ! the next lower aggregate size
     nsplit=(inext+comp(1,1))/inext_minus_1
     irest=mod(inext+comp(1,1),inext_minus_1)
     comp(1,:)=0
     comp(0,:)=0
     comp(inext,1)=comp(inext,1)-1
     comp(inext_minus_1,1)=nsplit
     if (comp(inext,1).gt.0) then  ! link up and down to inext_minus_1
         comp(inext,3)=inext_minus_1
         comp(inext_minus_1,2)=inext
     else                          ! if inext is cleared out, shift links down
         if (comp(inext,2).eq.Nsim+1) then
                 ! keeping track of overall progress
                 write (*,*) "# shifting highest i-mer from ",inext,"to ",inext_minus_1," at ",ncomp,Qtot
         endif
         comp(comp(inext,2),3)=inext_minus_1
         comp(inext_minus_1,2)=comp(inext,2)
         comp(inext,2)=0; comp(inext,3)=0
     endif
      
     ! deal with remainder
     iloop=inext_minus_1
     do while (irest.ne.0)
        irestnext=ilist(revlist(irest))  ! skip down in case irest-mers are not present
        comp(irestnext,1)=irest/irestnext
        comp(irestnext,2)=iloop
        comp(iloop,3)=irestnext
   !write (*,*) inext,irest,irestnext
        irest=mod(irest,irestnext)
        iloop=irestnext
     enddo
     comp(iloop,3)=0
     comp(0,2)=iloop
    
 
enddo


write (*,*) "# ncomp", ncomp


aveNout=aveNout/Qtot
write (*,*) "# Qtot",Qtot," most probable composition has likelihood ",biggest_P/Qtot

write (*,*)"# calculating high-N limit concentration distribution (tolerance:0.1%)"

! find bulk-limit distribution; start with a guess at the monomer concentration,
!  generate all i-mer concentrations, iterate to get the correct total
!  concentration.
testC1=aveNout(i)/V
error=10.0
bulk_limit_conc=0.0d0
do while (error.gt.0.001d0)
  testCtot=0.0d0
  dCdC1=0.0d0
  do i=1,nsim
     testCtot=testCtot+i*Kactual(i)*testC1**i
     dCdC1=dCdC1+i*i*Kactual(i)*testC1**(i-1)
  enddo
  error = abs(testCtot/rhotot-1.0)
  testC1=testC1+0.2d0*(rhotot-testCtot)/dCdC1
  if (testC1.lt.0.0d0) then  ! numerical instability; still should print out main results
        write (*,*) "# calculation failed - disregard column 3"
        bulk_limit_conc=-1.0d0
        exit
  endif
enddo
if (error.le.0.001d0) then  
     do i=1,nsim
        bulk_limit_conc(i)=Kactual(i)*testC1**i
     enddo
endif

write (*,*) "# i  <Ni/V>,Nsim   <Ni/V>,bulk   K_assoc  K_assoc(apparent) Ni_at_most_probable_comp"
do i=1,nsim
   if (aveNout(i).gt.0.0d0 .or. Kactual(i).gt.0.0d0) then
      Kapp(i)=(aveNout(i)/V)/( aveNout(1)/V )**i 
      write (*,*) i,aveNout(i)/V,bulk_limit_conc(i),Kactual(i),Kapp(i),most_prob_comp(i)
   endif
enddo



end

real*8 function Pcalc(Psave,facinv,q,comp,Nsim,Nmax,prestart)

implicit none

real*8:: facinv(0:Nmax),q(0:Nmax),Psave(0:Nmax+1)
integer:: comp(0:Nmax+1,3)
integer:: Nsim, ipoint,Nmax,prestart

Pcalc=Psave(prestart)
!write (*,*) "prestart",prestart,Psave(prestart)
ipoint=comp(prestart,3)  ! working down from largest unchanged aggregate 
do while (ipoint.gt.0)
   !write (*,*) ipoint,Psave(ipoint)
   Psave(ipoint)=Pcalc*(q(ipoint)**(comp(ipoint,1)))*facinv(comp(ipoint,1))
   Pcalc=Psave(ipoint)
   ipoint=comp(ipoint,3)
enddo
Psave(0)=Pcalc
!write (*,*) "p = ",Pcalc
if (Pcalc.lt.0.0d0) then
  Pcalc=1.0d0
  ipoint=comp(Nsim+1,3)
  do while (ipoint.gt.0)
   write (*,*) ipoint,q(ipoint),comp(ipoint,1),Pcalc,facinv(comp(ipoint,1))
   Pcalc=Pcalc*(q(ipoint)**(comp(ipoint,1)))*facinv(comp(ipoint,1))
   ipoint=comp(ipoint,3)
  enddo
  stop
endif
return
end
