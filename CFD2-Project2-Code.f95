program project2
implicit none
integer :: ndimn,npoin,nnode,nelem,neqns,nconi,nbfac,nfael,i
integer :: mesup,ielem,junk,ndegr,iface,mafac,mbfac,nafac,inode,istor
integer,allocatable :: inpoel(:,:)
REAL*8,allocatable :: coord(:,:),unkno0(:,:),unknoN(:,:),dt(:),geoel(:,:),geofa(:,:),resi(:)&
                    &,rhsel(:,:),rho(:),u(:),v(:),p(:),dens(:),vx(:),vy(:),press(:)
integer,allocatable :: esup1(:),esup2(:),psup1(:),psup2(:),esuel(:,:)
integer,allocatable :: intfac(:,:),bface(:,:)
REAL*8 :: x1,x2,x3,y1,y2,y3,length,xpoin,ypoin,cfl
REAL*8 :: pi,attack,gamma,rho0,mach0,vel0,csou0,p0,u0,v0,tol,entr0,delt
integer :: testcase,flag,iflag,poin1,poin2,accuracy,iteration,ipoin,ie,itermax,iter
REAL*8 :: sumrho,sumu,sumv,sump,sumarea,gam1,sumg


!------------------------------------------------------------------------------------------------------- READING THE GRID FILE



                                                                                                 !!!--------CHANGE THIS FOR DIFFERENT QUESTION
open(12,file='Grid.txt')
read(12,*) junk                         !Skipping intial lines
		do i=1,junk
          read(12,*)
        end do
read(12,*)
read(12,*) ndimn,nnode
read(12,*)
read(12,*) nelem,npoin,nbfac
read(12,*)

itermax=100000
iteration=0
ndegr=1
nconi=3
nfael=3
mafac=3*nelem
mbfac=3*nelem
neqns=3

allocate(inpoel(1:nnode,1:nelem),coord(1:ndimn,1:npoin))
allocate(bface(1:nconi,1:nbfac),esup2(1:npoin+1))
allocate(esuel(1:nnode,1:nelem),intfac(4,1:mafac))

!                               coord  ------> Grid Points
!                               inpoel ------> Element Point Connectivity Matrix
!                               bface  ------> Boundary Face
!                               esup2  ------> Element Surrounding Point
!                               esuel ------>  Element Surrounding Element
!                               intfac ------> Element Face Connectivity Matrix

do ielem=1,nelem
  read(12,*) junk,inpoel(1,ielem),inpoel(2,ielem),inpoel(3,ielem)
end do
read(12,*)

do ipoin=1,npoin
  read(12,*) junk,coord(1,ipoin),coord(2,ipoin)
end do
read(12,*)

do ipoin=1,npoin
  read(12,*)
end do
read(12,*)

do iface=1,nbfac
  read(12,*) junk,bface(1,iface),bface(2,iface),bface(3,iface)
end do

close(12)

!------------------------------------------------------------------------------------------------------- PRE PROCESSING







!------------------------------------------------------------------------------------------------------- PRE PROCESSING

call elesurpoi                                                                                 		 !Element Surrounding Points

call elesurele(esuel,esup1,esup2,inpoel,nfael,nelem,mesup,npoin,nnode)                               !Elements surrouding elements

call findfaces(nnode,mbfac,nelem,nbfac,nafac,inpoel,esuel,intfac,nfael,bface,nconi,mafac)            !Faces

print*,'PRE PROCESSING SUCCESS !!!!!!!'




!------------------------------------------------------------------------------------------------------- GEOMETRIC VARIABLES
allocate(geofa(3,nafac),geoel(3,nelem+nbfac))     
            
!                               geofa  ------> Face normals and length
!                               geoel ------>  Centroids and Area

do iface=1,nafac
  x1 = coord(1,intfac(3,iface))
  y1 = coord(2,intfac(3,iface))
  x2 = coord(1,intfac(4,iface))
  y2 = coord(2,intfac(4,iface))
  length = ((x2-x1)**2 + (y2-y1)**2)**0.5
  geofa(1,iface) = (y2-y1)/length                                          ! x-coordinate of Normal
  geofa(2,iface) = (x1-x2)/length                                          ! y-coordinate of Normal
  geofa(3,iface) = length                                                  ! Length of the face
end do

do ielem=1,nelem
  x1 = coord(1,inpoel(1,ielem))
  y1 = coord(2,inpoel(1,ielem))
  x2 = coord(1,inpoel(2,ielem))
  y2 = coord(2,inpoel(2,ielem))
  x3 = coord(1,inpoel(3,ielem))
  y3 = coord(2,inpoel(3,ielem))
  geoel(1,ielem) = (x1 + x2 + x3)/3.0d0                                          !x-coordinate of Centroid
  geoel(2,ielem) = (y1 + y2 + y3)/3.0d0                                          !y-coordinate of Centroid
  geoel(3,ielem) = (-(x3-x1)*(y2-y3) - (-(x2-x3))*(y3-y1))/2.0d0                 !Area of the element
end do
                                                                                 !Ghost cells
do iface = 1,nbfac
    
    x1 = coord(1,intfac(3,iface))
    y1 = coord(2,intfac(3,iface))
    x2 = coord(1,intfac(4,iface))
    y2 = coord(2,intfac(4,iface))
    
    poin1 = geoel(1,intfac(1,iface))
    poin2 = geoel(2,intfac(2,iface))
    xpoin = x1 + x2 - poin1
    ypoin = y1 + y2 - poin2
    geoel(1,intfac(2,iface)) = xpoin                                                          
    geoel(2,intfac(2,iface)) = ypoin
    geoel(3,intfac(2,iface)) = geoel(3,intfac(1,iface))
end do

print*,'GEOMETRIC MATRICES SUCCESS !!!!!!!'

!!---------Calculate averaged grid size-----------------
sumg=0.0d0
do ielem=1,nelem
  sumg=sumg+geoel(3,ielem)
end do
sumg=sumg/nelem
sumg=sqrt(sumg)
print *,sumg
!stop
!!-------------------------------------------------------
!-------------------------------------------------------------------------------------------------------SETTING UP INITIAL CONDITIONS AND INITIALIZATION

pi=0.25d0*atan(1.0d0)
attack=1.0d0											 !Angle of attack
gamma=1.4d0											     !Heat Capacity Ratio
rho0  = 1.0d0                                            !Density
mach0 = 0.85d0                                  		 !Mach Number
vel0 = 1.0d0                               				 !Velocity Assumed
csou0  = vel0/mach0                              	     !Speed of Sound
u0 = vel0*(cos(attack*(pi/180.0d0)))          	   	     !Velocity in the x-direction
v0 = vel0*(sin(attack*(pi/180.0d0)))          			 !Velocity in the y-direction
p0 =  (rho0*(u0*u0+v0*v0))/(gamma*mach0*mach0)                  		 !Pressure
tol = 1.0e-6                              			     !Tolerance
entr0 = p0/(rho0**gamma)                                 !Entropy
gam1=gamma - 1.0d0
allocate(unkno0(4,nelem+nbfac),unknoN(4,nelem+nbfac),resi(4),rhsel(4,nelem+nbfac))
allocate(dt(nelem))

!                               unkno0  ------> Unknown Conservative variables at Initial time
!                               unknoN ------>  Unknown Conservative variables at Final time

unkno0(1,:) = rho0
unkno0(2,:) = rho0*u0
unkno0(3,:) = rho0*v0
unkno0(4,:) = (p0/(gamma-1.0)) + (rho0*(u0*u0+v0*v0))*0.500

cfl=0.1d0
testcase = 1                   ! 1 -----> NACA 0012 and Cylinder
                               ! 2 -----> Channel bump flow

                                                                                                 !!!--------CHANGE THIS FOR DIFFERENT QUESTION


accuracy=1

                                                                                                 !!!--------CHANGE THIS FOR DIFFERENT QUESTION

                                                                                                 
							   ! Accuracy = 1   -----> FIRST ORDER
							   ! Accuracy = 2   -----> SECOND ORDER        !!!!!---Almost done----Debugging ongoing---!!!---Subroutines at the end---!!!
resi(:)=30.0d0
!-------------------------------------------------------------------------------------------------------RUNGE-KUTTA 3

call apply_boundcon(testcase,unkno0,intfac,bface,geofa,&
                   &nbfac,nafac,nelem,gamma,p0,u0,v0,rho0)

call calc_timestep(dt,delt,coord,unkno0,esuel,inpoel,geoel,gamma,&
                         &nnode,nelem,ndimn,npoin,nbfac,nfael)

open(22,file='Residual.txt')
write(22,*) 'VARIABLES="ITERATIONS","DENSITY","X-MOMENTUM DENSITY","Y-MOMENTUM DENSITY","ENERGY DENSITY"'
do while( (resi(1)>tol .or. resi(2)>tol .or. resi(3)> tol .or. resi(4)>tol))                                 !!!!VERY STRICT TOLERANCE CONDITION
  resi(:)=0.0d0

  call calc_timestep(dt,delt,coord,unkno0,esuel,inpoel,geoel,gamma,&
                         &nnode,nelem,ndimn,npoin,nbfac,nfael)

  call rk3(testcase,accuracy,unkno0,intfac,bface,geofa,geoel,&
           &rhsel,nafac,nbfac,nelem,gamma,p0,u0,v0,rho0,unknoN, &
           & delt)

  

  !-------------------------------------------------------------------------------------------------------Calculate Residual

  do ielem=1,nelem
    resi(:) = resi(:) + geoel(3,ielem)*(unknoN(:,ielem) - unkno0(:,ielem))*&
                                     &(unknoN(:,ielem) - unkno0(:,ielem))
  end do

  resi(:) = sqrt(resi(:))
  

  iteration=iteration+1
  print *, iteration,resi(1),resi(2),resi(3),resi(4)

  write(22,*) iteration,resi(1),resi(2),resi(3),resi(4)


  unkno0=unknoN
  
  if(iteration.gt.itermax) then
    goto 66
  end if
end do

close(22)
print*,'CONVERGENCE SUCCESS !!!!!!!'

66 continue
!-------------------------------------------------------------------------------------------------------Post Processing
allocate(rho(1:nelem),u(1:nelem),v(1:nelem),p(1:nelem))
allocate(dens(1:nelem),vx(1:nelem),vy(1:nelem),press(1:nelem))

rho(:)=unknoN(1,1:nelem)
u(:)=unknoN(2,1:nelem)/unknoN(1,1:nelem)
v(:)=unknoN(3,1:nelem)/unknoN(1,1:nelem)
p(:)=gam1*(unknoN(4,1:nelem)-0.5d0*rho(:)*(u(:)*u(:)+ v(:)*v(:)))

do ipoin=1,npoin
  
  sumrho=0.0d0
  sumu=0.0d0
  sumv=0.0d0
  sump=0.0d0
  sumarea=0.0d0

  do ie=esup2(ipoin)+1,esup2(ipoin+1)
    sumrho = sumrho + rho(esup1(ie))*geoel(3,esup1(ie))
    sumu = sumu + u(esup1(ie))*geoel(3,esup1(ie))
    sumv = sumv + v(esup1(ie))*geoel(3,esup1(ie))
    sump = sump + p(esup1(ie))*geoel(3,esup1(ie))
    sumarea = sumarea + geoel(3,esup1(ie))
  end do

  dens(ipoin)=sumrho/sumarea
  vx(ipoin)=sumu/sumarea
  vy(ipoin)=sumv/sumarea
  press(ipoin)=sump/sumarea

end do

!-------------------------------------------------------------------------------------------------------Write to FILE

open(88,file='results.txt')
write(88,*)'zone f=fepoint, N=',npoin,',E=',nelem,',et=triangle'
!write(88,*)'VARIABLES="X","Y","U","V","PRESSURE","DENSITY"'
do ipoin=1,npoin
  write(88,*) coord(1,ipoin),coord(2,ipoin),vx(ipoin),vy(ipoin),press(ipoin),dens(ipoin)
end do
do ielem=1,nelem
  write(88,*) inpoel(1,ielem),inpoel(2,ielem),inpoel(3,ielem)
enddo

close(88)
  
contains
!-------------------------------------------------------------------------------------------------------ELEMENT SURROUNDING POINTS
subroutine elesurpoi                                                                !This subroutine is contained in the program due to passing of an
implicit none                                                                       !allocatable array (esup1)
esup2=0

do ielem=1,nelem
  do inode=1,nnode
     
    esup2(inpoel(inode,ielem)+1)=esup2(inpoel(inode,ielem)+1)+1

  end do
end do


do ipoin=2,npoin+1

  esup2(ipoin)=esup2(ipoin)+esup2(ipoin-1)

end do


mesup = esup2(npoin+1)
allocate(esup1(mesup))

do ielem=1,nelem
  do inode=1,nnode

    ipoin=inpoel(inode,ielem)
    istor=esup2(ipoin)+1
    esup2(ipoin)=istor
    esup1(istor)=ielem

  end do
end do

do ipoin=npoin+1,2,-1

  esup2(ipoin)=esup2(ipoin-1)

end do

esup2(1)=0

end subroutine 


end program

!=====================================================================================ELEMENT SURROUNDING ELEMENTS=================

subroutine elesurele(esuel,esup1,esup2,inpoel,nfael,nelem,mesup,npoin,nnode)
implicit none
integer :: mesup,npoin,nelem,nnode,ielem,istor,ifael,nfael
integer :: ip1,ip2,jelem,icoun,ieadj
integer :: esup1(1:mesup),esup2(1:npoin+1),inpoel(1:nnode,1:nelem)
integer :: esuel(1:nfael,1:nelem),lpoin(1:npoin),lhelp(2,3)
lhelp=reshape((/2,3,3,1,1,2/),shape(lhelp))

lpoin=0
esuel=0
nfael=3

do ielem=1,nelem
  do ifael=1,nfael
    ieadj=0

    ip1=inpoel(lhelp(1,ifael),ielem)
    ip2=inpoel(lhelp(2,ifael),ielem)
    lpoin(ip1)=1
    lpoin(ip2)=1

    do istor=esup2(ip1)+1,esup2(ip1+1)
      jelem=esup1(istor)
      if(jelem.ne.ielem) then
        icoun=lpoin(inpoel(1,jelem))+lpoin(inpoel(2,jelem))+&
               &lpoin(inpoel(nnode,jelem))
               if(icoun.eq.nnode-1) then
                 ieadj=jelem
                 goto 10
               end if
      end if
    end do
    10 continue
    esuel(ifael,ielem)=ieadj
    lpoin(ip1)=0
    lpoin(ip2)=0
  end do
end do

end subroutine

!====================================================================================================FIND FACES=================
subroutine findfaces(nnode,mbfac,nelem,nbfac,nafac,inpoel,esuel,intfac,nfael,bface,nconi,mafac)
implicit none 
integer :: nnode,mbfac,nelem,nbfac,nafac,nfael,nconi,mafac
integer :: iface,ielem,inode,inode1,inode2,jelem,ip1,ip2
integer :: inpoel(1:nnode,1:nelem),intfac(4,mafac)
integer :: esuel(1:nfael,1:nelem),bface(1:nconi,1:nbfac)

nbfac = 0

do ielem  = 1,nelem
  do inode = 1,nnode
     inode1 = mod(inode,nnode) + 1
     inode2 = mod(inode1,nnode) + 1
     jelem = esuel(inode,ielem)
     if(jelem.eq.0)then
       nbfac = nbfac + 1
       if(nbfac.gt.mbfac)then
         print*, 'Please Increase mbfac =',nbfac
         stop
       end if
       intfac(1,nbfac) = ielem
       intfac(2,nbfac) = nelem + nbfac
       intfac(3,nbfac) = inpoel(inode1,ielem)
       intfac(4,nbfac) = inpoel(inode2,ielem)
       esuel(inode,ielem) = nelem + nbfac
     end if
   end do
end do


do iface = 1,nbfac

  ip1 = bface(1,iface)
  ip2 = bface(2,iface)
  do ielem = 1,nelem
    do inode = 1,nnode
      inode1 = mod(inode, nnode)+1
      inode2 = mod(inode1,nnode)+1
      if ((inpoel(inode1,ielem).eq.ip1).and.(inpoel(inode2,ielem).eq.ip2)) then
        intfac(1,iface) = ielem
        intfac(2,iface) = nelem + iface
        intfac(3,iface) = ip1
        intfac(4,iface) = ip2
      endif

    enddo
  enddo

enddo

nafac = nbfac


do ielem = 1,nelem
  do inode = 1,nnode
    inode1 = mod(inode, nnode)+1
    inode2 = mod(inode1,nnode)+1
    jelem = esuel(inode,ielem)
    if (jelem > ielem .and. jelem <= nelem) then
      nafac = nafac + 1
      intfac(1,nafac) = ielem
      intfac(2,nafac) = jelem
      intfac(3,nafac) = inpoel(inode1,ielem)
      intfac(4,nafac) = inpoel(inode2,ielem)
    endif
  end do
end do

print*,'   nafac=',nafac,'    nbfac=',nbfac
!stop

end subroutine

!=====================================================================================TIME STEP CALCULATIONS=================

subroutine calc_timestep(dt,delt,coord,unkno0,esuel,inpoel,geoel,gamma,&
                         &nnode,nelem,ndimn,npoin,nbfac,nfael)
implicit none
integer :: nnode,nelem,ndimn,npoin,nbfac,nfael
REAL*8 :: unkno0(4,nelem+nbfac),geoel(3,nelem+nbfac),coord(1:ndimn,1:npoin)
integer :: inpoel(1:nnode,1:nelem),esuel(nfael,nelem)
integer :: ielem,iface,poin1,poin2,icell,jcell
REAL*8 :: x1,x2,y1,y2,nx,ny,l,surflen,rho1,rho2,u1,u2,v1,v2,rhoe1,rhoe2,p1,p2
REAL*8 :: csou1,csou2,area,gamma,gam1,dt(nelem),eps,small,delt,cfl

small=10.0d0
gam1=gamma-1.0d0
eps=1.0e-6
cfl=0.1d0

do ielem = 1,nelem
  surflen = 0.0
  do iface=1,3
          
    poin1 = mod(iface,3)+1
    poin2 = mod(poin1,3)+1
        
    x1 = coord(1,inpoel(poin1,ielem))
    y1 = coord(2,inpoel(poin1,ielem))
    x2 = coord(1,inpoel(poin2,ielem))
    y2 = coord(2,inpoel(poin2,ielem))
           
    l = sqrt((x2-x1)**2 + (y2-y1)**2)
    nx = (y2-y1)/l
    ny = (x1-x2)/l
         
    icell = ielem
    jcell = esuel(iface,ielem)
         
    rho1 = unkno0(1,icell)
    rho2 = unkno0(1,jcell)
         
    u1 = unkno0(2,icell)/unkno0(1,icell)
    u2 = unkno0(2,jcell)/unkno0(1,jcell)
            
    v1 = unkno0(3,icell)/unkno0(1,icell)
    v2 = unkno0(3,jcell)/unkno0(1,jcell)

    rhoe1=unkno0(4,icell)
    rhoe2=unkno0(4,jcell)
           
    p1 = gam1* (rhoe1 - 0.5*rho1*(u1*u1 + v1*v1))
    p2 = gam1* (rhoe2 - 0.5*rho2*(u2*u2 + v2*v2))

    p1=max(eps,p1)
    p2=max(eps,p2)
    
    csou1 = sqrt(gamma*p1/rho1)
    csou2 = sqrt(gamma*p2/rho2)
          
    surflen = surflen + 0.5*(csou1+csou2  + abs((u1+u2)*nx + (v1+v2)*ny))
  end do
  area  = geoel(3,ielem)
  dt(ielem)  = cfl*area/surflen

end do
   
do ielem=1,nelem
  if(dt(ielem).lt.small)then
   small=dt(ielem)
  end if
end do

delt=small
print *,delt
end subroutine

!=====================================================================================BOUNDARY CONDITIONS=================

subroutine apply_boundcon(testcase,unkno0,intfac,bface,geofa,&
                         &nbfac,nafac,nelem,gamma,p0,u0,v0,rho0)
implicit none
integer ::testcase,nbfac,nafac,nelem
REAL*8 :: p0,u0,v0,rho0,csou0
integer :: intfac(4,nafac),bface(3,nbfac)
REAL*8 :: geofa(3,nafac),unkno0(4,nelem+nbfac)
integer :: bcell,gcell,iface
REAL*8 :: rho_g,rhoe_b,rho_b,rhou_b,rhov_b,gamma,gam1,mach_b,vt_b,vt0,vn0
REAL*8 :: u_g,v_g,p_g,u_b,v_b,p_b,csou_b,nx,ny,length,vn_b
csou0=sqrt(gamma*p0/rho0)
gam1=gamma-1.0
							   ! bcell -----> Boundary Cell
							   ! gcell -----> Ghost Cell
do iface=1,nbfac

  bcell = intfac(1,iface)
  gcell = intfac(2,iface)

  rho_b = unkno0(1,bcell)
  rhou_b = unkno0(2,bcell)
  rhov_b = unkno0(3,bcell)
  rhoe_b = unkno0(4,bcell)

  u_b=rhou_b/rho_b
  v_b=rhov_b/rho_b
  p_b=gam1*(rhoe_b - 0.5*rho_b*(u_b*u_b + v_b*v_b))
  csou_b=sqrt(gamma*p_b/rho_b)

  nx = geofa(1,iface)
  ny = geofa(2,iface)
  length = geofa(3,iface)

  if(testcase.eq.1)then
    if(bface(3,iface).eq.2)then
      rho_g = rho_b                                               !Solid Wall
      p_g = p_b
      u_g = u_b - 2.0*(u_b*nx + v_b*ny)*nx
      v_g = v_b - 2.0*(u_b*nx + v_b*ny)*ny
    else
      if(bface(3,iface).eq.4)then
      rho_g = rho0                                                 !Farfield
      p_g = p0
      u_g = u0
      v_g = v0
      end if
    end if
    
  else
    if(bface(3,iface).eq.2)then
      rho_g = rho_b                                               !Solid Wall
      p_g = p_b
      u_g = u_b - 2.0*(u_b*nx + v_b*ny)*nx
      v_g = v_b - 2.0*(u_b*nx + v_b*ny)*ny
    else if(bface(3,iface).eq.4)then
      vn_b  = u_b*nx + v_b *ny
      vt_b  = -u_b*ny + v_b*nx
      vn0 = u0*nx + v0*ny
      vt0 = -u0*ny + v0*nx
      mach_b  = vn_b/csou_b
      if(mach_b .gt.1.0)then                                          !Supersonic OUTflow
        p_g = p_b
        u_g = u_b
        v_g = v_b
        rho_g = rho_b
        elseif((mach_b.lt.1.0).and.(mach_b.gt.0.0))then                  !Subsonic OUTflow
          u_g = u_b - vn_b*nx + (0.5d0*(vn0+vn_b) + (csou_b-csou0)/gam1)*nx
          v_g = v_b - vn_b*ny + (0.5d0*(vn0+vn_b) + (csou_b-csou0)/gam1)*ny
          rho_g = rho_b*((0.25d0*gam1*(vn_b-vn0) + 0.5d0*(csou_b+csou0))/csou_b)**(2.0/gam1)
          p_g = ((0.25d0*gam1*(vn_b-vn0) + 0.5d0*(csou_b+csou0))**2)*rho_g/gamma
          elseif((mach_b.lt.0.0).and.(mach_b.gt.-1.0))then                  !Subsonic INflow
            u_g = u0 - vn0*nx + (0.5d0*(vn0+vn_b) + (csou_b-csou0)/gam1)*nx
            v_g = v0 - vn0*ny + (0.5d0*(vn0+vn_b) + (csou_b-csou0)/gam1)*ny
            rho_g = rho0*((0.25d0*gam1*(vn_b-vn0) + 0.5d0*(csou_b+csou0))/csou0)**(2.0d0/gam1)
            p_g =((0.25d0*gam1*(vn_b-vn0) + 0.5d0*(csou_b+csou0))**2)*rho_g/gamma
            else
              rho_g = rho0                                             !Supersonic INflow
              p_g = p0
              u_g = u0
              v_g = v0
          end if
        end if

    
  end if

  unkno0(1,gcell) = rho_g
  unkno0(2,gcell) = rho_g*u_g
  unkno0(3,gcell) = rho_g*v_g
  unkno0(4,gcell) = (p_g/gam1) + 0.5*rho_g*(u_g*u_g + v_g*v_g)

end do

end subroutine

!=====================================================================================Van Leer Flux Calculation=================

subroutine calc_flux(accuracy,unkno0,rhsel,intfac,geofa,nafac,&
                     &nelem,nbfac,gamma)
implicit none
integer :: accuracy,nafac,nelem,nbfac
REAL*8 :: gamma
REAL*8 :: unkno0(4,nelem+nbfac),rhsel(4,nelem+nbfac),geofa(3,nafac)
integer :: intfac(4,nafac)
integer :: iface,icell,jcell
REAL*8 :: rho1,rhou1,rhov1,rhoe1,rho2,rhou2,rhov2,rhoe2,u1,u2,v1,v2,p1,p2
REAL*8 :: csou1,csou2,mach1,mach2,vn1,vn2,flux(4),flux1(4),flux2(4)
REAL*8 :: eps,gam1,nx,ny

rhsel=0.0d0
gam1=gamma- 1.0d0
eps=1.0e-8

if(accuracy.eq.1) then

  do iface=1,nafac

    icell=intfac(1,iface)
    jcell=intfac(2,iface)

    nx=geofa(1,iface)
    ny=geofa(2,iface)
        
    rho1=unkno0(1,icell)
    rhou1=unkno0(2,icell)
    rhov1=unkno0(3,icell)
    rhoe1=unkno0(4,icell)

    rho2=unkno0(1,jcell)
    rhou2=unkno0(2,jcell)
    rhov2=unkno0(3,jcell)
    rhoe2=unkno0(4,jcell)

    u1=rhou1/rho1
    v1=rhov1/rho1
    p1=max(eps,gam1*(rhoe1- 0.5d0*rho1*(u1*u1 + v1*v1)))
    csou1=sqrt(gamma*p1/rho1)
    vn1=u1*nx + v1*ny
    mach1=vn1/csou1

    u2=rhou2/rho2
    v2=rhov2/rho2
    p2=max(eps,gam1*(rhoe2- 0.5d0*rho2*(u2*u2 + v2*v2)))
    csou2=sqrt(gamma*p2/rho2)
    vn2=u2*nx + v2*ny
    mach2=vn2/csou2

    if(mach1.lt.-1.0d0)then
      flux1(:)=0.0d0
      else
        if(mach1.gt.1.0d0)then
          flux1(1)=rho1*vn1
          flux1(2)=rhou1*vn1 + p1*nx
          flux1(3)=rhov1*vn1 + p1*ny
          flux1(4)=vn1*(rhoe1+p1)
          else
            flux1(1)=0.25d0*csou1*rho1*(mach1+1.0d0)*(mach1+1.0d0)
            flux1(2)=flux1(1) * (u1 + nx*(csou1+csou1-vn1)/gamma)
            flux1(3)=flux1(1) * (v1 + ny*(csou1+csou1-vn1)/gamma)
            flux1(4)=flux1(1) * (0.5d0*(u1*u1 + v1*v1 - vn1*vn1) &
                        &       + 0.5d0*((gam1*vn1 + csou1+csou1)* &
                        &   (gam1*vn1 + csou1+csou1))/(gamma*gamma - 1.0d0))
        end if
    end if

    if(mach2.gt.1.0d0)then
      flux2(:)=0.0
      else
        if(mach2.lt.-1.0d0)then
          flux2(1)=rho2*vn2
          flux2(2)=rhou2*vn2 + p2*nx
          flux2(3)=rhov2*vn2 + p2*ny
          flux2(4)=vn2*(rhoe2+p2)
          else
            flux2(1)=-0.25d0*csou2*rho2*(mach2-1.0d0)*(mach2-1.0d0)
            flux2(2)=flux2(1) * (u2 + nx*(-csou2-csou2-vn2)/gamma)
            flux2(3)=flux2(1) * (v2 + ny*(-csou2-csou2-vn2)/gamma)
            flux2(4)=flux2(1) * (0.5d0*(u2*u2 + v2*v2 - vn2*vn2) &
                        &       + 0.5d0*((gam1*vn2 - csou2-csou2)* &
                        &   (gam1*vn2 - csou2-csou2))/(gamma*gamma - 1.0d0))
        end if
    end if
  
  flux(:)=(flux1(:)+flux2(:))*geofa(3,iface)
  
  rhsel(:,icell)=rhsel(:,icell)-flux(:)
  rhsel(:,jcell)=rhsel(:,jcell)+flux(:)

    
end do
end if

end subroutine


!=====================================================================================RUNGE KUTTA 3 TVD=================

subroutine rk3(testcase,accuracy,unkno0,intfac,bface,geofa,geoel,&
               &rhsel,nafac,nbfac,nelem,gamma,p0,u0,v0,rho0,unknoN, &
               & delt)
implicit none
integer :: accuracy,nafac,nelem,nbfac,testcase,ielem
integer :: intfac(4,nafac),bface(3,nbfac)
REAL*8 :: geoel(3,nelem+nbfac),unkno0(4,nelem+nbfac),rhsel(4,nelem+nbfac)
REAL*8 :: geofa(3,nafac),Ustg1(4,nelem+nbfac),Ustg2(4,nelem+nbfac)
REAL*8 :: unknoN(4,nelem+nbfac),delt,gamma,p0,u0,v0,rho0



                   
call apply_boundcon(testcase,unkno0,intfac,bface,geofa,&
                   &nbfac,nafac,nelem,gamma,p0,u0,v0,rho0)
                   
call calc_flux(accuracy,unkno0,rhsel,intfac,geofa,nafac,&
              &nelem,nbfac,gamma)
do ielem=1,nelem
Ustg1(:,ielem) = unkno0(:,ielem) + &
                     &(delt*rhsel(:,ielem)/geoel(3,ielem))                   ! Stage 1
end do



call apply_boundcon(testcase,Ustg1,intfac,bface,geofa,&
                   &nbfac,nafac,nelem,gamma,p0,u0,v0,rho0)
                   
call calc_flux(accuracy,Ustg1,rhsel,intfac,geofa,nafac,&
              &nelem,nbfac,gamma)
              
do ielem=1,nelem                   
Ustg2(:,ielem) = (0.75d0*unkno0(:,ielem)) + &
                     &(0.25d0*(Ustg1(:,ielem) + &
                    &(delt*rhsel(:,ielem)/geoel(3,ielem))))! Stage 2
end do

call apply_boundcon(testcase,Ustg2,intfac,bface,geofa,&
                   &nbfac,nafac,nelem,gamma,p0,u0,v0,rho0)
                   
call calc_flux(accuracy,Ustg2,rhsel,intfac,geofa,nafac,&
              &nelem,nbfac,gamma)

do ielem=1,nelem              
unknoN(:,ielem) = (unkno0(:,ielem)) + &
                      &(2.0d0*(Ustg2(:,ielem) + &
                      &(delt*rhsel(:,ielem)/geoel(3,ielem))))! Stage 3 
unknoN(:,ielem) = unknoN(:,ielem)/3.0d0
end do
call apply_boundcon(testcase,unknoN,intfac,bface,geofa,&
                   &nbfac,nafac,nelem,gamma,p0,u0,v0,rho0)

end subroutine rk3

!=================================================================    
!===============================================================

!subroutine poisurpoi(psup1,psup2,esup1,esup2,inpoel,mesup,npoin,nnode,nelem)
!implicit none
!integer :: mesup,npoin,nnode,inode,ielem,ipoin,istor,iesup,nelem
!integer :: esup1(1:mesup),esup2(1:npoin+1)
!integer :: jpoin,inpoel(1:nnode,1:nelem)
!integer :: psup1(1:mesup),psup2(1:npoin+1),lpoin(1:npoin)

!lpoin=0
!psup2(1)=0
!istor=0

!do ipoin=1,npoin
!  lpoin(ipoin)=ipoin
!  do iesup=esup2(ipoin)+1,esup2(ipoin+1)

!    ielem=esup1(iesup)

!    do inode=1,nnode
!      jpoin=inpoel(inode,ielem)
!      if(lpoin(jpoin).ne.ipoin) then
!        istor=istor+1
!        psup1(istor)=jpoin
!        lpoin(jpoin)=ipoin
!      end if
!    end do
!  end do
  
!  psup2(ipoin+1)=istor
!end do

!end subroutine
!=================================================================    
!===============================================================
subroutine reconstruct(lesq,geoel,esuel,nelem,nbfac,nfael)
integer :: nelem,nbfac,nfael
REAL*8 :: geoel(3,nelem+nbfac),lesq(3,nelem)
REAL*8 :: wtx2,wtxy,wty2,dx,dy,inv,x1,x2,y1,y2
integer :: ielem,inode,esuel(nfael,nelem)

do ielem=1,nelem
  wtx2=0.0d0
  wty2=0.0d0
  wtxy=0.0d0
  do inode=1,3
    x1=geoel(1,ielem)
    x2 = geoel(1,esuel(inode,ielem))
    y1 = geoel(2,ielem)
    y2 = geoel(2,esuel(inode,ielem))
    length = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
    inv = 1/length                                      ! inverse length
    dx = x1-x2
    dy = y1-y2
    wtx2 = wtx2 + (inv**2)*(dx*dx)                     !weight x^2
    wty2 = wty2 + (inv**2)*(dy*dy)                     !weight y^2
    wtxy = wtxy + (inv**2)*(dx*dy)                     !weight xy
    end do
  lesq(1,ielem) = wtx2
  lesq(2,ielem) = wty2
  lesq(3,ielem) = wtxy
  end do
end subroutine

!=================================================================                                     !!!___FIND GRADIENTS AND VAN ALBADA____!!!
!===============================================================

!!!!!-------THIS VAN ALBADA RECONSTRUCTION IS SUPPOSED TO BE INSIDE THE FLUXING SUBROUTINE IN THE IF (ACCURACY=2 ) CONDITION



  
  



