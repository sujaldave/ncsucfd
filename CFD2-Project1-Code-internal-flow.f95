program potential2d
implicit none
integer :: i,j,k,nelem,nface,npoin,ndimn,ntype,ngeel,ip1,ip2,ip3,ngefa,vx,vy,itermax,iter
real :: junk,a3,b3,total,tol,res,res1,sum1,sum2,sum3
real, allocatable :: inpoel(:,:),coord(:,:),bface(:,:),geoel(:,:),geofa(:,:),lhspo(:,:),rhspo(:)
real, allocatable :: phiold(:),phinew(:),vx_ele(:),vy_ele(:),vx_poin(:),vy_poin(:),totvel(:)
junk=0.0
ndimn=2
ntype=3
ngeel=5
ngefa=2
vx=1
vy=0

!-------------READ THE GRID FILE----------
 
!SKIP THE FIRST SEVEN LINES
open(12,file='grid file for channel.txt')
do i=1,7
  read(12,*)
end do

read(12,*) nelem,npoin,nface
read(12,*)

allocate(coord(ndimn,npoin),inpoel(ntype,nelem),bface(ntype,nface),geoel(ngeel,nelem),geofa(ngefa,nface))
allocate(lhspo(npoin,npoin),rhspo(npoin),phiold(npoin),phinew(npoin),vx_ele(nelem),vy_ele(nelem))
allocate(vx_poin(npoin),vy_poin(npoin),totvel(npoin))
do i=1,nelem
  read(12,*) junk,inpoel(1,i),inpoel(2,i),inpoel(3,i)
end do

read(12,*)

do i=1,npoin
  read(12,*) junk,coord(1,i),coord(2,i)
end do

do i=1,npoin+2
  read(12,*)
end do

do i=1,nface
  read(12,*) junk,bface(1,i),bface(2,i),bface(3,i)
end do

close(12)

!-------------SETTING UP DATA STRUCTURE----------

do i=1,nelem
  ip1=inpoel(1,i)
  ip2=inpoel(2,i)
  ip3=inpoel(3,i)
  geoel(5,i)=((coord(1,ip1)*(coord(2,ip2)-coord(2,ip3))) &
  			& -(coord(2,ip1)*(coord(1,ip2)-coord(1,ip3)))&
            & +((coord(1,ip2)*coord(2,ip3))-(coord(1,ip3)*coord(2,ip2))))
  geoel(5,i)=(abs(geoel(5,i))) !D
  geoel(1,i)=(coord(2,ip2)-coord(2,ip3))/geoel(5,i) !a1
  geoel(2,i)=(coord(2,ip3)-coord(2,ip1))/geoel(5,i) !a2
  geoel(3,i)=(coord(1,ip3)-coord(1,ip2))/geoel(5,i) !b1
  geoel(4,i)=(coord(1,ip1)-coord(1,ip3))/geoel(5,i) !b2
end do

do i=1,nface
  ip1=bface(1,i) !left point
  ip2=bface(2,i) !right point
  geofa(1,i)=coord(2,ip2)-coord(2,ip1) !nx-comp
  geofa(2,i)=coord(1,ip1)-coord(1,ip2) !ny-comp
end do

!-------------ASSEMBLING GLOBAL MATRIX----------

lhspo=0.0
rhspo=0.0

do i=1,nelem
  ip1=inpoel(1,i)
  ip2=inpoel(2,i)
  ip3=inpoel(3,i)

  a3=-geoel(1,i)-geoel(2,i)
  b3=-geoel(3,i)-geoel(4,i)

  lhspo(ip1,ip1)=lhspo(ip1,ip1)+(((geoel(1,i)*geoel(1,i))+(geoel(3,i)*geoel(3,i)))*0.5*geoel(5,i))
  lhspo(ip1,ip2)=lhspo(ip1,ip2)+(((geoel(1,i)*geoel(2,i))+(geoel(3,i)*geoel(4,i)))*0.5*geoel(5,i))
  lhspo(ip1,ip3)=lhspo(ip1,ip3)+(((geoel(1,i)*a3)+(geoel(3,i)*b3))*0.5*geoel(5,i))
  lhspo(ip2,ip1)=lhspo(ip2,ip1)+(((geoel(2,i)*geoel(1,i))+(geoel(4,i)*geoel(3,i)))*0.5*geoel(5,i))
  lhspo(ip2,ip2)=lhspo(ip2,ip2)+(((geoel(2,i)*geoel(2,i))+(geoel(4,i)*geoel(4,i)))*0.5*geoel(5,i))
  lhspo(ip2,ip3)=lhspo(ip2,ip3)+(((geoel(2,i)*a3)+(geoel(4,i)*b3))*0.5*geoel(5,i))
  lhspo(ip3,ip1)=lhspo(ip3,ip1)+(((a3*geoel(1,i))+(b3*geoel(3,i)))*0.5*geoel(5,i))
  lhspo(ip3,ip2)=lhspo(ip3,ip2)+(((a3*geoel(2,i))+(b3*geoel(4,i)))*0.5*geoel(5,i))
  lhspo(ip3,ip3)=lhspo(ip3,ip3)+((((a3*a3)+(b3*b3)))*0.5*geoel(5,i))
 
end do


do i=1,nface
  ip1=bface(1,i)
  ip2=bface(2,i)
  if(bface(3,i)==4) then
    rhspo(ip1)=rhspo(ip1)+(((vx*geofa(1,i))+(vy*geofa(2,i)))*0.5)
    rhspo(ip2)=rhspo(ip2)+(((vx*geofa(1,i))+(vy*geofa(2,i)))*0.5)
  end if
  if(bface(3,i)==2) then
    rhspo(ip1)=rhspo(ip1)+0.0
    rhspo(ip2)=rhspo(ip2)+0.0
  end if
end do

!-------------GAUSS SIEDEL SOLVER----------

phiold=0.0
phinew=0.0
itermax=10000
tol=0.0001


do iter=1,itermax
  
  do i=1,npoin
    total=0.0
    
        
    do j=1,i-1
      total=total+(lhspo(i,j)*phinew(j))
    end do
    
    do j=i+1,npoin
      total=total+(lhspo(i,j)*phiold(j))
    end do
    
    phinew(i)=(rhspo(i)-total)/lhspo(i,i)

  end do

  res=0.0
  do k=1,npoin
    res = res + ((phinew(k)-phiold(k))*(phinew(k)-phiold(k)))
  end do


  res = sqrt(res)
  if (iter==1) then
  	res1 = res
  end if
  
  res=res/res1
  
  if (res<=tol) then
    exit
  end if
  
  phiold=phinew
end do
   
!-------------POST PROCESSING----------

do i=1,nelem
  ip1=inpoel(1,i)
  ip2=inpoel(2,i)
  ip3=inpoel(3,i)

  a3=-geoel(1,i)-geoel(2,i)
  b3=-geoel(3,i)-geoel(4,i)

  vx_ele(i)=geoel(1,i)*phinew(ip1)+geoel(2,i)*phinew(ip2)+a3*phinew(ip3)
  vy_ele(i)=geoel(3,i)*phinew(ip1)+geoel(4,i)*phinew(ip2)+b3*phinew(ip3)
end do

do i=1,npoin
  sum1=0.0
  sum2=0.0
  sum3=0.0

  do j=1,nelem
    ip1=inpoel(1,j)
    ip2=inpoel(2,j)
    ip3=inpoel(3,j)

    if(ip1==i .or. ip2==i .or. ip3==i) then
      sum1=sum1+(vx_ele(j)*geoel(5,j)*0.5)
      sum2=sum2+(vy_ele(j)*geoel(5,j)*0.5)
      sum3=sum3+(geoel(5,j)*0.5)
    end if
  end do

  vx_poin(i)=sum1/sum3
  vy_poin(i)=sum2/sum3
  totvel(i)=sqrt((vx_poin(i)*vx_poin(i))+(vy_poin(i)*vy_poin(i)))
end do 	

open(20, file = 'Test1_final.txt')
write(20,*) 'VARIABLES = "X","Y","phi","vx","vy","v"'
write(20,*) 'ZONE F = FEPOINT, N=',npoin,', E=',nelem,', ET= TRIANGLE'

do i=1,npoin
  write(20,*) coord(1,i),coord(2,i),phinew(i),vx_poin(i),vy_poin(i),totvel(i)
end do

do i=1,nelem
  write(20,*) inpoel(1,i),inpoel(2,i),inpoel(3,i)
end do
close(20)

open(22, file= 'topfacevelocity.txt')
open(24, file= 'bottomfacevelocity.txt')

do i=1,nface
  if(bface(3,i)==2)then
    ip1=bface(1,i)
    ip2=bface(2,i)
    if(coord(2,ip1)<0.5)then
      write(24,*) coord(1,ip1),totvel(ip1)
    end if
    if(coord(2,ip1)>0.5)then
      write(22,*) coord(1,ip1),totvel(ip1)
    end if
    if(coord(2,ip2)<0.5)then
      write(24,*) coord(1,ip2),totvel(ip2)
    end if
    if(coord(2,ip2)>0.5)then
      write(22,*) coord(1,ip2),totvel(ip2)
    end if
  end if
end do
close(22)
close(24)

end program potential2d
