program md
implicit none
integer(4)L,ndim,nbase,indist,i,modj,modk

ndim=16
nbase=3
L=ndim*nbase

indist=nbase*L

modj=1
do i=1,indist
 
   modk=mod(i,nbase)
   if(modk==0)modk=nbase 
   write(6,*) i,modj,modk
   if(mod(i,L)==0)then
    modj=modj+1
   endif 

end do

write(6,*) 'compare',indist,L*L

end program md
