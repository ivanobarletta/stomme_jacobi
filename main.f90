program stommel_jacobi
! Use finite difference method to solve the stommel problem
! See Stommel, 1948 (Westward intensification of wind-driven Ocean Currents)
! Solution is reached iteratively with Jacobi method

USE netcdf

implicit none

integer         :: nx,ny
integer         :: ji,jj,jj_mid

double precision                :: zLx,zLy
double precision                :: zdeltax, zdeltay
double precision                :: zbeta, zalpha

double precision, allocatable   :: x(:),y(:),xv(:)

double precision, allocatable   :: psi(:,:), psin(:,:)
double precision, allocatable   :: vpsi(:)

double precision, allocatable   :: zf(:,:)

double precision                :: zpi

double precision                :: ztau0,zH,zrho0,zr 

double precision                :: zdeltaxy2
double precision                :: za1,za2,za3,za4,za5,zka 

integer                         :: iter, niters
double precision                :: zres

integer                         :: ierr
integer                         :: ncid, xdimid, ydimid
integer                         :: xvarid,yvarid,psivarid,vpsivarid,zfvarid

character(len=80)               :: outname

integer                         :: iunit, ios, rc,fu

iunit   = 777

NAMELIST/params/zH,zr,zbeta,ztau0,zrho0,zLx,zLy,zdeltax,zdeltay,niters,outname
open (action='read', file='namelist.str', iostat=rc, newunit=fu)
read (nml=params, iostat=rc, unit=fu)

write(6,*) 'zh',zH

zpi     = 4.0D0 * ATAN(1.0D0)

write(6,*) 'zLx',zLx
write(6,*) 'zLy',zLy

nx      = zLx / zdeltax + 1
ny      = zLy / zdeltay + 1

zalpha  = (zpi*ztau0) / (zH*zrho0*zLy) 

! allocate variables
allocate(x   (1:nx),y(1:ny) )
allocate(xv  (1:nx-1) )
allocate(vpsi(1:nx))
allocate(zf  (1:nx   ,1:ny) )
allocate(psi (1:nx   ,1:ny) )
allocate(psin(1:nx   ,1:ny) )

psi(:,:)        = 0.0D0 
psin(:,:)       = 0.0D0 

! build coordinates
do ji = 1,nx
   x(ji) = (ji-1)*zdeltax
end do
do jj = 1,ny
   y(jj) = (jj-1)*zdeltay
end do
do ji = 1,nx-1
   xv(ji) = 0.5*zdeltax + (ji-1)*zdeltax
end do

write(6,*) 'x'
write(6,*) x
write(6,*) 'y'
write(6,*) y

! find jj of y=Ly/2
jj_mid  = ny / 2 

write(6,*) 'zr   ', zr
write(6,*) 'zH   ', zH
write(6,*) 'zbeta', zbeta
write(6,*) 'zrho0', zrho0
write(6,*) 'ztau0', ztau0


! initialize zf(x,y)
do jj = 1,ny
   do ji = 1,nx
      zf(ji,jj) = SIN((zpi*y(jj))/(zLy)) 
   end do
end do

za1     = 2*( zr / zdeltax**2 + zr / zdeltay**2 )
za2     = zr / zdeltax**2 + 0.5*zbeta / zdeltax
za3     = zr / zdeltax**2 - 0.5*zbeta / zdeltax
za4     = zr / zdeltay**2

do iter = 1,niters
   zres = 0.0D0
   do jj = 2,ny-1
      do ji = 2,nx-1
         psin(ji,jj) = (za2/za1)*psi(ji+1,jj) + (za3/za1)*psi(ji-1,jj) + (za4/za1)*psi(ji,jj+1) + (za4/za1)*psi(ji,jj-1) - (zalpha/za1)*zf(ji,jj)      
      end do
   end do

   zres = SUM( (psi-psin)**2)
   write(6,*) 'iter,zres', iter,zres 
   psi(:,:) = psin(:,:)
   if ( zres .lt. 1D-12 ) exit
end do

! calculate v=-d psi / dx
vpsi(1)         = -1.0D20
vpsi(nx)        = -1.0D20
do ji=2,nx-1
   vpsi(ji) = -(psi(ji+1,jj_mid)-psi(ji-1,jj_mid)) / (2*zdeltax)
end do

! write netcdf
ierr    = nf90_create(outname,NF90_NETCDF4,ncid)

ierr    = nf90_def_dim(ncid,'x',nx,xdimid)
ierr    = nf90_def_dim(ncid,'y',ny,ydimid)

ierr    = nf90_def_var(ncid,'x',NF90_DOUBLE,xdimid,xvarid)
ierr    = nf90_def_var(ncid,'y',NF90_DOUBLE,ydimid,yvarid)

ierr    = nf90_def_var(ncid,'psi',NF90_DOUBLE,(/xdimid,ydimid/),psivarid)
ierr    = nf90_def_var(ncid,'zf' ,NF90_DOUBLE,(/xdimid,ydimid/),zfvarid)

ierr    = nf90_def_var(ncid,'vpsi',NF90_DOUBLE,xdimid,vpsivarid)
ierr    = nf90_put_att(ncid, vpsivarid, "_FillValue", -1.0D20)

ierr    = nf90_enddef(ncid)

ierr    = nf90_put_var(ncid,xvarid,x,start=(/1/),count=(/nx/))
ierr    = nf90_put_var(ncid,yvarid,y,start=(/1/),count=(/ny/))

ierr    = nf90_put_var(ncid,zfvarid  ,zf  ,start=(/1,1/),count=(/nx,ny/))

ierr    = nf90_put_var(ncid,psivarid ,psi ,start=(/1,1/),count=(/nx,ny/))

ierr    = nf90_put_var(ncid,vpsivarid,vpsi,start=(/1/)  ,count=(/nx/))

ierr    = nf90_close(ncid)

end program stommel_jacobi
