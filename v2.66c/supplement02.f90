module mPlot
character*10 :: file_ext='png'
integer :: num_pix_x=900, num_pnt_x=3000
!margin in unit of pixel
integer :: margin_left=80, margin_right=120, margin_bot=50, margin_top=100
integer :: num_pix_y, num_pnt_y
real :: pnt_pix
integer :: curg_grp=1
character(len=5) :: color_set(0:8)= &
    (/'SMALL', 'VGA', 'RAIN', 'SPEC', 'GREY', 'RRAIN', 'RSPEC', 'RGREY','TEMP' /)
end module mPlot

subroutine OutputPlot
use mPlot
use paraset2    , only :num_cmesh
use paraset1
use paraset4 , only: num_group

integer i,j,k, msf(3), ipx
integer IsLoc
integer :: set_number=1
integer ing, ing_t

if( PngSizeFactor .gt. 9.0) PngSizeFactor=9
if( PngSizeFactor .lt. 0.1) PngSizeFactor=0.1

num_pix_x=int(num_pix_x*PngSizeFactor)
margin_left=margin_left*PngSizeFactor
margin_right=margin_right*PngSizeFactor
margin_bot=margin_bot*PngSizeFactor
margin_top=margin_top*PngSizeFactor

if(num_pix_x .gt. num_pnt_x) num_pnt_x=num_pix_x

msf=0
ing_t=1
!mat
ipx=abs(MatSrcFlx)
if(ipx .eq. 1 .or. ipx .eq. 3 .or. ipx .eq. 5 .or. ipx .ge. 7) then
msf(1)=1
endif

!src
if(ipx .eq. 2 .or. ipx .eq. 3 .or.  ipx .ge. 6) then
msf(2)=1
endif

!flux
if(ipx .ge. 4) then
msf(3)=1
endif

msf_loop : do k=1 ,3 

if (msf(k) .eq. 0) cycle
ing_t=1

if(k .eq. 3) then
  if(num_group_eff .eq. 0) then
    write(*,*) 'WARNING, Failed to plot flux because no flux files loaded'
    return
  else
    ing_t=num_group
  endif
endif

do ing=1, ing_t

if(k .eq. 3) then
if (grp_flag(ing,set_number) .eq. 1) then
 curg_grp=ing
else
 cycle
endif
endif

if (MatSrcFlx .lt. 0) then
ipx=-k
else
ipx=k
endif

IsLoc=0
!z cut plot
if(IsZPlot .eq. 1 .and. num_zloc .ge. 1) then
  do i=1, num_zloc
    Call PlotLoc(i,3,ipx)
  enddo
  IsLoc=IsLoc+1
endif

!y cut plot
if(IsYPlot .eq. 1 .and. num_yloc .ge. 1) then
  do i=1, num_yloc
    Call PlotLoc(i,2,ipx)
  enddo
  IsLoc=IsLoc+1
endif

!x cut plot
if(IsXPlot .eq. 1 .and. num_xloc .ge. 1) then
  do i=1, num_xloc
    Call PlotLoc(i,1,ipx)
  enddo
  IsLoc=IsLoc+1
endif

if(IsLoc .ne. 0) cycle
!z mid-level plot
if(MidPlot .eq. 1 .or. MidPlot .eq. 3 .or. MidPlot .eq. 5 .or. MidPlot .ge. 7) then
 if(IsZSingle .eq. 0) then
  do j=1, num_cmesh(3)
   call PlotLoc(-j, 3,ipx)
  enddo
 else
  do j=z_start, z_end
    call PlotLoc(-j, 3,ipx)
  enddo
 endif

endif
!y mid-level plot
if(MidPlot .eq. 2 .or. MidPlot .eq. 3 .or.  MidPlot .ge. 6) then
do j=1, num_cmesh(2)
call PlotLoc(-j, 2,ipx)
enddo
endif
!x mid-level plot
if(MidPlot .ge. 4) then
do j=1, num_cmesh(1)
call PlotLoc(-j, 1,ipx)
enddo
endif

enddo !ing
enddo msf_loop

end subroutine


!plot an X, Y, or Z cut
!plotnum : if >0, is the plot number of a series of locations
!          if <0, plot at corarse mesh mid-plane
!xyzplot : axis index
!          =1, plot at x=...
!          =2, plot at y=...
!          =3, plot at z=...
!msfplot:  plot mat, sre or flx, if <0 , will output a data file
!          =1, plot mat. distribution
!          =2, plot fixed src dist.
!          =3, plot flux dist.
subroutine PlotLoc(plotnum, xyzplot, msfplot)
use DISLIN
use mPlot
use files
use funs
use ErrControl
use paraset1
use paraset2
use paraset4
use constants

implicit none
integer plotnum, xyzplot, msfplot

integer i,j,k,ini,inj
integer idx, ipx
character*80 cstr, cstr_t1, cstr_t2
character*3 :: xyzset='xyz' 
character(len=1) :: bigxyz(3)=(/'X','Y','Z'/)
real rtemp(4),x_size(2), y_size(2), xy_ratio, fin_size
integer itemp(4), org_pos(2) 
real, dimension(:,:), allocatable :: zmat,  z_d
real, dimension(:), allocatable :: xray, yray, zray

real plot_loc 
integer idx_i, idx_j, cm_ijk(3)
integer :: set_number=1
real :: fm_xyz(3), lim_val(2)=0, x_d(2), y_d(2), viewpnt(3)

character(LEN=20) :: what2plot(3)=(/'Mat Num', 'Fixed Src', 'Flux'/) 
logical ex
real :: fseg

!plot mat, src, or flx
idx=abs(xyzplot)
ipx=abs(msfplot)
! setting up png file name
!if(Is3DPlot .eq. 0) then
curfile%ext=file_ext
!else
!curfile%ext='tif'
!endif

curfile%fullname=''
curfile%dir='.'

!truncate proname to the first 3 letters
if(prbname_len .le. 3) then
write(form,"('(A',I0,',A2,I0, A1)')") prbname_len 
else
write(form,"('(A',I0,',A2,I0, A1)')") 3
endif

!adding axis index, plotnum, and 'i' (specified pos by -plotx) 'c' default mid-plane pos
if(plotnum .gt. 0) then
  plot_loc=xyzloc(plotnum, idx)
  write(cstr,form ) trim(prbname), '_'//xyzset(idx:idx), plotnum, 'i'
elseif(plotnum .lt. 0) then  !mid-level
  plot_loc=(cgdim(idx)%boundary(abs(plotnum)-1)+cgdim(idx)%boundary(abs(plotnum)) )/2.0
  write(cstr,form ) trim(prbname), '_'//xyzset(idx:idx), abs(plotnum), 'c'
endif

!adding mat, src, or flux
select case( ipx)
case (1)
  write(curfile%name, "(A,A1)") trim(cstr),'m'
case (2)
  write(curfile%name, "(A,A1)") trim(cstr),'s'
case (3)
  write(curfile%name, "(A,A1,I0)") trim(cstr),'f', curg_grp   
end select

!locate the chosen axis CM index 
cm_ijk(idx)=0
do k=num_cmesh(idx), 0, -1
 if(plot_loc .gt. cgdim(idx)%boundary(k)) then
   cm_ijk(idx)=k+1
   exit
 endif 
enddo
if(cm_ijk(idx) .gt. num_cmesh(idx) .or. cm_ijk(idx) .eq. 0) then
write(*, "('Warning: Plot position out of boundary at ', A1, '=', f10.2, ' plotting skipped')") &
   xyzset(idx:idx), plot_loc
return
endif
!change to relative position
!plot_loc=plot_loc-cgdim(idx)%boundary(cm_ijk(idx)-1)

!selecting the remaining two axses 
select case(idx)
case (1)
idx_i=2
idx_j=3
case (2)
idx_i=1
idx_j=3
case (3)
idx_i=1
idx_j=2
end select


cur_filename=Getfullname(curfile,0)
inquire(file=cur_filename, exist=ex)
!if file exist, delete it
if(ex) then
 open(unit=211, file=cur_filename)
 close(211, status='DELETE')
endif
write(*,"('Plotting slice at ', A1, '=', f10.2)") xyzset(idx:idx), plot_loc
call SETFIL(trim(cur_filename))

!if(Is3DPlot .eq. 0) then:
Call METAFL ('PNG')
!else
!Call METAFL ('TIFF')
!endif


CALL SCRMOD ('REVERS')

x_size(1)=cgdim(idx_i)%boundary(0)
x_size(2)=cgdim(idx_i)%boundary(num_cmesh(idx_i))
y_size(1)=cgdim(idx_j)%boundary(0)
y_size(2)=cgdim(idx_j)%boundary(num_cmesh(idx_j))
xy_ratio=(y_size(2)-y_size(1))/(x_size(2)-x_size(1))
num_pnt_y=int(num_pnt_x*xy_ratio)

!add margin
pnt_pix=1.0*num_pnt_x/(num_pix_x - margin_left -margin_right)
itemp(3)=int(margin_left*pnt_pix)+num_pnt_x + int(margin_right*pnt_pix)
itemp(4)=int(margin_bot*pnt_pix)+num_pnt_y + int(margin_top*pnt_pix)
itemp(1)=int(margin_left*pnt_pix)
itemp(2)=int(margin_bot*pnt_pix)
!figure real size num_pnt_x, num_pnt_y
!itemp(1)=int(num_pnt_x/6)
!itemp(2)=int(num_pnt_y/8)
!itemp(3)=itemp(1)*8+itemp(1)
!add y margin itemp(2)
!itemp(4)=itemp(2)*10

xy_ratio=1.0*itemp(4)/itemp(3)
num_pix_y=int(num_pix_x*xy_ratio)

!if(num_pix_y .lt. 400) then 
!  num_pix_y=400 
!  xy_ratio=1.0*num_pix_y/num_pix_x
!  itemp(4)=int(itemp(3)*xy_ratio)
!endif

call WINSIZ(num_pix_x, num_pix_y)
call page(itemp(3), itemp(4))

CALL DISINI()
!plot a page border
!     CALL PAGERA()
!Font
!CALL BMPFNT ('COMPLEX')
CALL COMPLX()
!CALL SIMPLX()
!CALL DUPLX()
!CALL HELVES()

!org_pos(1)=itemp(1)
!org_pos(2)=itemp(2)+num_pnt_y
org_pos(1)=int(margin_left*pnt_pix)
org_pos(2)=int(margin_top*pnt_pix) + num_pnt_y
Call AXSPOS(org_pos(1),org_pos(2))
Call AX3LEN(num_pnt_x,num_pnt_y,num_pnt_y)


CALL NAME(what2plot(ipx),'Z')
CALL NAME(bigxyz(idx_i)//' axis','X')
CALL NAME(bigxyz(idx_j)//' axis','Y')
if(Is3DPlot .eq. 0) then
CALL LABELS ('NONE', 'XY')
else
CALL LABELS ('FLOAT', 'XY')
endif
if(ipx .eq. 2 .or. ipx .eq. 3) then
 CALL LABELS ('FEXP', 'Z')
else
 CALL LABDIG(-1,'Z')
endif
if(IsLogPlot .eq. 1 .and. ipx .eq. 3) then
  CALL AXSSCL('LOG', 'Z')
  CALL LABELS ('LOG', 'Z')
  CALL LABDIG(-1,'Z')
endif

!CALL NAMDIS (int(itemp(1)*0.3), 'X')
!CALL NAMDIS (int(itemp(2)*1.0), 'Y')
CALL NAMDIS (int(itemp(2)*0.7), 'X')
CALL NAMDIS (int(itemp(1)*0.6), 'Y')
      
!CALL SETVLT ('SPEC') 
CALL SETVLT (trim(color_set(ColorMap))) 
CALL TICKS (0, 'XYZ')
if (IsLogPlot .eq. 1) then
CALL TICKS(2, 'Z')
endif
!viewpnt for 3-d plot
viewpnt(1)=4
viewpnt(2)=4
viewpnt(3)=4
! CALL NOBAR
select case(ipx)
case (1) 
 if(num_material .gt. 1) then	
 i=int((num_material-1)/20)+1
   if(Is3DPlot .eq. 0) then
     Call GRAF3(x_size(1),x_size(2),x_size(1),x_size(2),&
       y_size(1), y_size(2), y_size(1), y_size(2),&
	   1.0, num_material*1.0,1.0,i*1.0 )
   else
     CALL VIEW3D(viewpnt(1),viewpnt(2), viewpnt(3),'ABS')
     Call GRAF3D(x_size(1),x_size(2),x_size(1),(x_size(2)-x_size(1))/5,&
       y_size(1), y_size(2), y_size(1), (y_size(2)-y_size(1))/5,&
	   1.0, num_material*1.0,1.0,i*1.0 )
   endif
 
 else
   if(Is3DPlot .eq. 0) then
     Call GRAF3(x_size(1),x_size(2),x_size(1),x_size(2),&
      y_size(1), y_size(2), y_size(1), y_size(2),&
	  1.0, 2.0,1.0,1.0 )
   else
     CALL VIEW3D(viewpnt(1),viewpnt(2), viewpnt(3),'ABS')
	 Call GRAF3D(x_size(1),x_size(2),x_size(1),(x_size(2)-x_size(1))/5,&
      y_size(1), y_size(2), y_size(1), (y_size(2)-y_size(1))/5,&
	  1.0, 2.0,1.0,1.0 )
   endif
 
 endif

case(2)
call GetCutLimit(plot_loc, idx, ipx, lim_val)
if(lim_val(1) .eq. 0 ) then
  if(max_src .gt. 0) then
    lim_val(1)=max_src
  else
    lim_val(1)=1.0
  endif
endif
if(Is3DPlot .eq. 0) then
 Call GRAF3(x_size(1),x_size(2),x_size(1),x_size(2),&
      y_size(1), y_size(2), y_size(1), y_size(2),&
	  0.0, lim_val(1),0.0,lim_val(1)/20 )
else
  CALL VIEW3D(viewpnt(1),viewpnt(2), viewpnt(3),'ABS')
  Call GRAF3D(x_size(1),x_size(2),x_size(1),(x_size(2)-x_size(1))/5,&
      y_size(1), y_size(2), y_size(1), (y_size(2)-y_size(1))/5,&
	  0.0, lim_val(1),0.0,lim_val(1)/20 )
endif
 
case (3)
call GetCutLimit(plot_loc, idx, ipx, lim_val)
if(lim_val(1) .eq. 0 ) then
  if(glb_fluxset(set_number)%max_glb_flux_grp(curg_grp) .gt. 0) then
    lim_val(1)=glb_fluxset(set_number)%max_glb_flux_grp(curg_grp)
  else
   lim_val(1)=1.0
  endif
endif

if(flx_max_flag .ne. 0) lim_val(1)=flx_max
if(flx_min_flag .ne. 0) lim_val(2)=flx_min

if(IsLogPlot .eq. 1) then
 lim_val(1)=log10(lim_val(1))
 if(lim_val(1) .le. 0) then
 lim_val(1)=int(lim_val(1))*1.0
 else
 lim_val(1)=int(lim_val(1))+1.0
 endif

 if(lim_val(2) .gt. 0) then
   lim_val(2)=log10(lim_val(2))
   if(lim_val(2) .le. 0) then
     lim_val(2)=int(lim_val(2))-1.0
   else
    lim_val(2)=int(lim_val(2))*1.0
   endif
 else
   lim_val(2)=log10(tiny)
 endif
if(lim_val(1) - lim_val(2) .le. 0) lim_val(2)=lim_val(1)-1
 fseg=1
else
 fseg=lim_val(1)/20
endif


if(Is3DPlot .eq. 0) then
  Call GRAF3(x_size(1),x_size(2),x_size(1),x_size(2),&
      y_size(1), y_size(2), y_size(1), y_size(2),&
	  lim_val(2), lim_val(1),lim_val(2),fseg )
 else
  CALL VIEW3D(viewpnt(1),viewpnt(2), viewpnt(3),'ABS')
  Call GRAF3D(x_size(1),x_size(2),x_size(1),(x_size(2)-x_size(1))/5,&
      y_size(1), y_size(2), y_size(1), (y_size(2)-y_size(1))/5,&
	  lim_val(2), lim_val(1),lim_val(2),fseg )
endif

!Call GRAF3(x_size(1),x_size(2),x_size(1),x_size(2),&
!  y_size(1), y_size(2), y_size(1), y_size(2),&
! 0.0, glb_fluxset(set_number)%max_glb_flux_grp(curg_grp),0.0,glb_fluxset(set_number)%max_glb_flux_grp(curg_grp)/20 )

end select


!title text 
if(ipx .eq. 3) then
write(cstr_t1, "(A, I0)") trim(what2plot(ipx))//' distribution for group=', curg_grp
else if(ipx .eq. 2) then
write(cstr_t1, "(A)") trim(what2plot(ipx))//' distribution'
else
write(cstr_t1, "(A, I0, A1)") trim(what2plot(ipx))//' distribution (num_material= ', num_material, ')'
endif
CALL TITLIN(cstr_t1,1)

if(IsLogPlot .eq. 0) then
write(cstr_t2, "(A1,'-', A1,' plot at ',A1,'=', f10.2, 2x, '(', A1,'-level: ', I0, ')' )") &
  xyzset(idx_i:idx_i), xyzset(idx_j:idx_j), xyzset(idx:idx), plot_loc,xyzset(idx:idx), cm_ijk(idx)
else
write(cstr_t2, "(A1,'-', A1,' LOG plot at ',A1,'=', f10.2, 2x, '(', A1,'-level: ', I0, ')' )") &
  xyzset(idx_i:idx_i), xyzset(idx_j:idx_j), xyzset(idx:idx), plot_loc,xyzset(idx:idx), cm_ijk(idx)
endif

CALL TITLIN(cstr_t2,2)

CALL TITLE()

if(Is3DPlot .eq. 0) then
do i=0, num_cmesh(idx_i)
  write(cstr, "(f10.2)") cgdim(idx_i)%boundary(i)
  Call ADDLAB (adjustl(cstr), cgdim(idx_i)%boundary(i), 2, 'X')
enddo
      
do j=0, num_cmesh(idx_j)
  write(cstr, "(f10.2)") cgdim(idx_j)%boundary(j)
  Call ADDLAB (adjustl(cstr), cgdim(idx_j)%boundary(j), 2, 'Y')
enddo

endif
!  CALL SETRGB(0.7,0.7,0.7)
!  CALL GRID(4,4)
fm_xyz(idx)=plot_loc 
if(msfplot .lt. 0) then
 open (unit=67, file=trim(curfile%name)//'.cat')
 write(67, "(A)") trim(cstr_t1)
 write(67, "(A)") trim(cstr_t2)
 write(67, "('  x-pos         y-pos         z-pos           ', A)") trim(what2plot(ipx))
endif

! ploting CMs
do j=1, num_cmesh(idx_j)
do i=1, num_cmesh(idx_i)
 cm_ijk(idx_i)=i
 cm_ijk(idx_j)=j
      
 rtemp(1)=cgdim(idx_i)%boundary(i-1) +zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_i)/2.01
 rtemp(2)=cgdim(idx_i)%boundary(i) -zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_i)/2.01
 rtemp(3)=cgdim(idx_j)%boundary(j-1) +zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_j)/2.01
 rtemp(4)=cgdim(idx_j)%boundary(j) -zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_j)/2.01
if(Is3DPlot .eq. 0) then	  
  call SURSZE (rtemp(1), rtemp(2), rtemp(3), rtemp(4))
endif 
itemp(1)=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx_i)
itemp(2)=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx_j)      
allocate(zmat(itemp(1), itemp(2)))
allocate(xray(itemp(1)), yray(itemp(2)) )

!k=(zlevel(znum)%cm_zlev(i,j)%fine(3)-1)/2+1
!if(IsZPlot .eq. 1) then
fin_size=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx)
k=int((plot_loc-cgdim(idx)%boundary(cm_ijk(idx)-1))/fin_size)+1
if(k .gt. zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx)) &
  k=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx)

select case (idx)
case (3)
 select case (ipx)
  case (1)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(:,:,k)
  case (2)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(:,:,k)
  case (3)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(:,:,k) 
   if(IsLogPlot .eq. 1) then
   
   do inj=1, itemp(2)
   do ini=1, itemp(1)
     if(  zmat(ini,inj) .le. 0) then
	   zmat(ini,inj)=10**lim_val(2)
	 endif
    enddo
    enddo
   endif
 end select

 xray=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%cm_x
 yray=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%cm_y

case (2)

select case (ipx)
  case (1)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(:,k,:)
  case (2)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(:,k,:)
  case (3)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(:,k,:) 
   if(IsLogPlot .eq. 1) then
   
   do inj=1, itemp(2)
   do ini=1, itemp(1)
     if(  zmat(ini,inj) .le. 0) then
	   zmat(ini,inj)=10**lim_val(2)
      endif

    enddo
    enddo
   endif
 end select

xray=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%cm_x
yray=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%cm_z

case (1)

select case (ipx)
  case (1)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(k,:,:)
  case (2)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(k,:,:)
  case (3)
   zmat=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(k,:,:) 
   if(IsLogPlot .eq. 1) then
   
   do inj=1, itemp(2)
   do ini=1, itemp(1)
     if(  zmat(ini,inj) .le. 0) then
	   zmat(ini,inj)=10**lim_val(2)
     
	 endif 
    enddo
    enddo
   endif

 end select

xray=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%cm_y
yray=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%cm_z

end select

!output data file
if(msfplot .lt. 0) then

 do inj=1, itemp(2)
 do ini=1, itemp(1)
   fm_xyz(idx_i)=xray(ini)
   fm_xyz(idx_j)=yray(inj)
   write(67, "(3(ES12.5,2x),2x, ES14.6)") fm_xyz, zmat(ini,inj)
 enddo
 enddo
 
endif 


x_size(1)=0.0
x_size(2)=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_i)
y_size(1)=0.0
y_size(2)=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_j)

if(Is3DPlot .eq. 0) then

CALL TRFREL (x_size, y_size, 2)
itemp(3)=abs(int(x_size(2)-x_size(1)))
itemp(4)=abs(int(y_size(2)-y_size(1)))


CALL SETRES(itemp(3),itemp(4))

if(itemp(1) .eq. 1) then
 allocate(zray(itemp(2)))
 zray(:)=zmat(1,:)
 call CURVY3(rtemp(1),yray,zray, itemp(2) )
 deallocate(zray)
elseif (itemp(2) .eq. 1) then
 allocate(zray(itemp(1)))
 zray(:)=zmat(:,1)
 call CURVX3(xray,rtemp(3),zray,itemp(1) )
 deallocate(zray)
else
 CALL CRVMAT(zmat,itemp(1),itemp(2),1,1)
endif

else  !plot 3D
!  CALL SHDMOD('SMOOTH','SURFACE')
  if(itemp(1) .eq. 1 .and. itemp(2) .eq. 1) then
    x_d(1)=xray(1)-0.25*x_size(2)
	x_d(2)=xray(1)+0.25*x_size(2)
	y_d(1)=yray(1)-0.25*y_size(2)
	y_d(2)=yray(1)+0.25*y_size(2)
	allocate( z_d(2,2) )
	z_d=zmat(1,1)
	CALL SURSHD(x_d, 2, y_d, 2, z_d)
    deallocate(z_d)
  elseif (itemp(1) .eq. 1) then
    x_d(1)=xray(1)-0.25*x_size(2)
	x_d(2)=xray(1)+0.25*x_size(2)
    allocate( z_d( 2,itemp(2) ) )
    z_d(1,:)=zmat(1,:)
	z_d(2,:)=z_d(1,:)
	CALL SURSHD(x_d, 2, yray, itemp(2), z_d) 
    deallocate(z_d)
  elseif (itemp(2) .eq. 1) then
    y_d(1)=yray(1)-0.25*y_size(2)
	y_d(2)=yray(1)+0.25*y_size(2)
    allocate( z_d( itemp(1),2 ) )
    z_d(:,1)=zmat(:,1)
	z_d(:,2)=z_d(:,1)
	CALL SURSHD(xray, itemp(1), y_d, 2, z_d) 
    deallocate(z_d)
  
  else
    CALL SURSHD(xray,itemp(1),yray,itemp(2),zmat)
  endif

endif

!draw lines
!vertical 
if(Is3DPlot .eq. 0 .and. IsDrawFm .eq. 1) then
if(itemp(3) .gt. 8 .and. itemp(4) .gt. 8) then
rtemp(1)=cgdim(idx_i)%boundary(i-1)+zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_i)
rtemp(2)=cgdim(idx_j)%boundary(j-1)
rtemp(3)=rtemp(1)
rtemp(4)=cgdim(idx_j)%boundary(j)
do k=1, zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx_i)-1
 call rline(rtemp(1),rtemp(2),rtemp(1), rtemp(4))
 rtemp(1)=rtemp(1)+zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_i) 
enddo

! horizon
rtemp(1)=cgdim(idx_i)%boundary(i-1)
rtemp(2)=cgdim(idx_j)%boundary(j-1)+zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_j)
rtemp(3)=cgdim(idx_i)%boundary(i)
rtemp(4)=rtemp(idx_j)
do k=1, zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx_j)-1
 call rline(rtemp(1),rtemp(2),rtemp(3), rtemp(2))
 rtemp(2)=rtemp(2)+zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx_j) 
enddo
endif

endif !Is3DPlot

deallocate(zmat, xray, yray)
enddo
enddo

if(msfplot .lt. 0) close (67)

!draw coarse mesh line
!CALL LINWID (2)
if(Is3DPlot .eq. 0) then
call SETRGB (1.0, 1.0, 1.0)
!vertical
!rtemp(1)=cgdim(1)%boundary(1)
rtemp(2)=cgdim(idx_j)%boundary(0)
!rtemp(3)=rtemp(1)
rtemp(4)=cgdim(idx_j)%boundary(num_cmesh(idx_j))
do k=1, num_cmesh(idx_i)-1
 rtemp(1)=cgdim(idx_i)%boundary(k)
 call rline(rtemp(1),rtemp(2),rtemp(1), rtemp(4))
enddo

! horizon
rtemp(1)=cgdim(idx_i)%boundary(0)
!rtemp(2)=cgdim(2)%boundary(1)
rtemp(3)=cgdim(idx_i)%boundary(num_cmesh(idx_i))
!rtemp(4)=rtemp(2)
do k=1, num_cmesh(idx_j)-1
 rtemp(2)=cgdim(idx_j)%boundary(k)
 call rline(rtemp(1),rtemp(2),rtemp(3), rtemp(2))
enddo

call LINWID (1)
CALL COLOR('FORE')
endif

CALL DISFIN()

write(READLOG,"('Plotting slice at ', A1, '=', f10.2,', ', A, ', file: ', A)") &
  xyzset(idx:idx), plot_loc, trim(cstr_t1), trim(cur_filename)

end subroutine 

!find the max and min for a plane
!max: lim_out(1),  min: lim_out(2)
subroutine GetCutLimit(cut_loc, icut, msfcut,lim_out)
use paraset4
use paraset2
use mPlot , only : curg_grp

real cut_loc, lim_out(2)
integer icut, msfcut

integer i,j,k
integer cm_ijk(3), idx, idx_i, idx_j, ipx
real fin_size, max_cut, min_cut
integer :: set_number=1

!locate the chosen axis CM index 
cm_ijk(icut)=0
idx=abs(icut)
ipx=abs(msfcut)
lim_out(1)=0
lim_out(2)=0

do k=num_cmesh(idx), 0, -1
 if(cut_loc .gt. cgdim(idx)%boundary(k)) then
   cm_ijk(idx)=k+1
   exit
 endif 
enddo
if(cm_ijk(idx) .gt. num_cmesh(idx) .or. cm_ijk(idx) .eq. 0) then
write(*, "('Warning: cut position out of boundary ')") 
return
endif

!selecting the remaining two axses 
select case(idx)
case (1)
idx_i=2
idx_j=3
case (2)
idx_i=1
idx_j=3
case (3)
idx_i=1
idx_j=2
end select
do j=1, num_cmesh(idx_j)
do i=1, num_cmesh(idx_i)

cm_ijk(idx_i)=i
cm_ijk(idx_j)=j
 

fin_size=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%finsize(idx)
k=int((cut_loc-cgdim(idx)%boundary(cm_ijk(idx)-1))/fin_size)+1
if(k .gt. zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx)) &
  k=zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fine(idx)


select case (idx)
case (3)
 select case (ipx)
  case (1)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(:,:,k))
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(:,:,k))
  case (2)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(:,:,k))
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(:,:,k))
  case (3)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(:,:,k)) 
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(:,:,k)) 
 end select



case (2)

select case (ipx)
  case (1)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(:,k,:))
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(:,k,:))
  case (2)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(:,k,:))
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(:,k,:))
  case (3)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(:,k,:) )
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(:,k,:) )
 end select


case (1)

select case (ipx)
  case (1)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(k,:,:))
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%mat_matrix(k,:,:))
  case (2)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(k,:,:))
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%src_matrix(k,:,:))
  case (3)
   max_cut=maxval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(k,:,:)) 
   min_cut=minval(zlevel(cm_ijk(3))%cm_zlev(cm_ijk(1),cm_ijk(2))%fluxset(set_number)%locg_set(curg_grp)%locg_flux(k,:,:)) 
 end select

end select
if (i .eq. 1 .and. j .eq. 1) then
lim_out(1)=max_cut
lim_out(2)=min_cut
else
if(max_cut .gt. lim_out(1)) lim_out(1)=max_cut
if(min_cut .lt. lim_out(2)) lim_out(2)=min_cut
endif

enddo
enddo

return
end subroutine
