!subcode3.f90   assignment of petran input files varibles
!subroutines included in this file
!          AssignPentranVaribles
!          TransCmIndex:   transfer coarse mesh nuber to x,y,z index number.  
subroutine AssignPentranVaribles
use paraset1
use paraset2
use paraset3
use paraset4
use paraset5
use files
use ErrControl
use funs
use fido

integer i,j,k,m,i_var, mi
integer temp(2)
integer cm_i,cm_j,cm_k

integer mfront, mback, mleft, mright
!real*8 :: src_total=0
!real*8 t
logical ex

real sca_guess

!for generating the Gemetry block
cgdim(1)%var1='x'
cgdim(1)%var2='xmesh='
cgdim(1)%var3='ixfine='
cgdim(1)%var4='ixmed='

cgdim(1)%int3=7
cgdim(1)%int4=6

cgdim(2)%var1='y'
cgdim(2)%var2='ymesh='
cgdim(2)%var3='jyfine='
cgdim(2)%var4='jymed='
cgdim(2)%int3=7
cgdim(2)%int4=6

cgdim(3)%var1='z'
cgdim(3)%var2='zmesh='
cgdim(3)%var3='kzfine='
cgdim(3)%var4='kzmed='
cgdim(3)%int3=7
cgdim(3)%int4=6


sn=sn_order
group=num_group
num_mat=num_material

num_dir=sn*(sn+2)
maxmem_byte=tot_num_fm*num_dir*group*8
maxmem=maxmem_byte/1e6+10
maxlin=0
maxarr=0
nctlim=0
max_fido_char=0
max_input_len=0

maxleg=pn_order
maxcrs=0
maxmed=0
maxfin=0

maxmmc=1
maxfmc=1

do i=1,3
  if(num_cmesh(i) .gt. maxcrs) then
    maxcrs=num_cmesh(i)
  endif
enddo


!refined maxfin maxmed maxmmc maxfmc
do k=1, num_cmesh(3)
do j=1, num_cmesh(2)
do i=1, num_cmesh(1)
   temp(1)=zlevel(k)%cm_zlev(i,j)%fine(1)*zlevel(k)%cm_zlev(i,j)%fine(2)*&
        zlevel(k)%cm_zlev(i,j)%fine(3)   
   temp(2)=zlevel(k)%cm_zlev(i,j)%med(1)*zlevel(k)%cm_zlev(i,j)%med(2)*&
          zlevel(k)%cm_zlev(i,j)%med(3)
   if(temp(1) .gt. maxfmc)  maxfmc=temp(1) 
   if(temp(2) .gt. maxmmc)  maxmmc=temp(2)
   do mi=1,3
     if(maxfin .lt. zlevel(k)%cm_zlev(i,j)%fine(mi) ) maxfin=zlevel(k)%cm_zlev(i,j)%fine(mi)
     if(maxmed .lt. zlevel(k)%cm_zlev(i,j)%med(mi) )  maxmed=zlevel(k)%cm_zlev(i,j)%med(mi)
   enddo
enddo
enddo
enddo

!Upon Request maxmmc=maxfmc
maxmmc=maxfmc
maxmed=maxfin


ngroup(1)=group
ngroup(2)=1
ngroup(3)=0

maxleg=max(leg_order,pn_order)
!block 2	
allocate(flxini(num_cm))
if(s_format .eq. 2) then
 sca_guess=1.0
else
 sca_guess=0.0
endif

do i=1 ,num_cm
   flxini(i)=sca_guess
enddo

allocate(mathmg(num_cm))
mathmg=0

! block 3
legord=pn_order
legoxs=leg_order
nxtyp=xs_format
ihm=len_table
iht=3

if(ihm .eq. num_group+3) then
ihs=4
elseif(ihm .eq. 2*num_group+2) then
ihs=iht+num_group
elseif(ihm  .gt. num_group+3) then
ihs=ihm-num_group+1
else
 ihs=4
 write(warning_message,"('ihm is invalid' )") 
 call DisplayMsg
endif

ihng=0

!allocate(chig(num_mat*group))
!chig=0

! block4
ncoupl=1
nrdblk=0

!tolin_len=1 ! or =num_cm
!allocate(tolin(tolin_len))
if (sec_flag (3) .eq. 0) then
tolin(1)=1E-3
! tolout
tolout(1)=1E-5
!maxitr
maxitr(1)=70  !max_no_of_iter=50
maxitr(2)=10  !criticality_inner_limit=10
rkdef=1.0
endif

!med_grid_multiplier
tolout(2)=1E-3


!maxitr=50
methit=1
methac=1
!nzonrb
number_of_zones=num_cm
damping_fact=0.999
skip_iter=0

dtwmxw=0.97
nquit=4
num_source=0
do k=1, num_cmesh(3)
  do j=1, num_cmesh(2)
    do i=1, num_cmesh(1)
     if(cmdiff_global .ne. 0) zlevel(k)%cm_zlev(i,j)%ndmeth=ndmeth_global
     if(zlevel(k)%cm_zlev(i,j)%s_flag .ge. 1) num_source=num_source+1
    enddo
  enddo
enddo
nprtyp=num_source
! block 5


allocate(nsdef(num_source))
nsdef=0   ! volume source
allocate(nscmsh(num_source))
m=0
do k=1, num_cmesh(3)
  do j=1, num_cmesh(2)
    do i=1, num_cmesh(1)
    if(zlevel(k)%cm_zlev(i,j)%s_flag .ge. 1) then 
      m=m+1
      nscmsh(m)=(k-1)*num_cmesh(1)*num_cmesh(2)+(j-1)*num_cmesh(1)+i
!    else
!     if(allocated(zlevel(k)%cm_zlev(i,j)%src_matrix)) &
!	 deallocate(zlevel(k)%cm_zlev(i,j)%src_matrix)
	endif
    enddo
  enddo  
enddo

allocate(sref(3*num_source))
sref=0.0

!do i=1, num_source
!  do j=1,3
!     sref((i-1)*3+j)=0.0
!  enddo  
!enddo

allocate(serg(group*num_source))
if(num_source .gt. 0) then
curfile%name=prbname
curfile%ext='spc'
curfile%fullname=''
curfile%dir=input_dir
cur_filename=Getfullname(curfile,1)

inquire(file=cur_filename,exist=ex)
if(ex ) then
 write(warning_message,"('Reading Src spectrum(.spc) file: ',A)") trim(cur_filename)

 call DisplayMsg
 open(unit=44,file=cur_filename)
 
 read(44,*,iostat=i_var) (serg(i),i=1,group)
 close(44)
 
 if (i_var .ne. 0) then
  write(warning_message,"('read error, skiped file: ',A)") trim(cur_filename)
  call TrapReadError(0)  
  do j=1,num_source
  serg((j-1)*group+1)=1.0
  do i=2, group
   serg((j-1)*group+i)=0
  enddo
 enddo
 else
  write(warning_message,"('.spc file done: ',A)" ) trim(cur_filename)
  call DisplayMsg
 do i=2,num_source
 do j=1,group
   serg((i-1)*group+j)=serg(j)
 enddo
 enddo
 endif
else 
 do j=1,num_source
  serg((j-1)*group+1)=1.0
  do i=2, group
   serg((j-1)*group+i)=0
  enddo
 enddo
endif
endif

allocate(chig(num_mat*group))
chig=0
do i=1, num_mat
  chig((i-1)*group+1)=1
enddo

curfile%name=prbname
curfile%ext='chi'
curfile%fullname=''
curfile%dir=input_dir
cur_filename=Getfullname(curfile,1)

inquire(file=cur_filename,exist=ex)
if(ex ) then
 
 write(warning_message,"('Reading fission spectrum(.chi) file:' ,A)") trim(cur_filename)
 call DisplayMsg

 open(unit=45,file=cur_filename)
 read(45,*,iostat=i_var) (chig(i),i=1,group)
 close(45)
 if(i_var .ne. 0) then
   write(warning_message,"('read error, skiped file: ',A)") trim(cur_filename)
   call TrapReadError(0)
   chig=0
 else
 write(warning_message,"('.chi file done: ',A)" ) trim(cur_filename)
 call DisplayMsg
 do i=2,num_mat
 do j=1,group
   chig((i-1)*group+j)=chig(j)
 enddo
 enddo
endif
endif

if(num_source .gt. 0) then
allocate(smag(num_source))
do i=1, num_source
  call TransCmIndex(nscmsh(i),cm_i,cm_j,cm_k)
  smag(i)=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%sum_act
!  src_total=src_total+smag(i)
enddo 
endif

!read in fission src
call ReadInFis

! block  6 BCs
!Penmsh:  ibback, ibfrnt, jbright,jbleft,kbottom,ktop
!pentran: ibback, ibfrnt, jbeast,jbwest,kbsout,kbnort
bc_char(1)='ibback='
bc_char(2)='ibfrnt='
bc_char(3)='jbeast='
bc_char(4)='jbwest='
bc_char(5)='kbsout='
bc_char(6)='kbnort='

do i=1,6
 if(bcs(i) .eq. 1) then
  allocate(albedo(i)%albedo_grp(group))    
  do j=1, group
   albedo(i)%albedo_grp(j)=1
  enddo
 endif
enddo

!block 7

allocate (nadump(1))
nadump=0
!dump angular flux at boundary
mfront=0
mback=0
mleft=0
mright=0
!allocate( nadump(2*num_cmesh(3)*(num_cmesh(1)+num_cmesh(2)) ) 
!allocate (front(num_cmesh(3)*num_cmesh(1)) )
!allocate (back(num_cmesh(3)*num_cmesh(1)) )
!allocate (left(num_cmesh(3)*num_cmesh(2)) )
!allocate (right(num_cmesh(3)*num_cmesh(2)) )
!do i=1, num_cm
!  if(cmesh(i)%index_cm(2) .eq. 1) then
!    mfront=mfront+1
!	front(mfront)=i
!  endif
!  if(cmesh(i)%index_cm(2) .eq. num_cmesh(2)) then
!    mback=mback+1
!	back(mback)=i
!  endif
!  if(cmesh(i)%index_cm(1) .eq. 1) then
!    mleft=mleft+1
!	left(mleft)=i
!  endif
!  if(cmesh(i)%index_cm(1) .eq. num_cmesh(1)) then
!    mright=mright+1
!	right(mright)=i
!  endif
!enddo

return
end subroutine

!transfer cm number to i,j, k index

subroutine TransCmIndex(cm_t, cm_i, cm_j, cm_k)
USE paraset1
USE paraset2
USE paraset3
USE paraset4

integer cm_t, cm_i,cm_j,cm_k

integer remain1
if(cm_t .le. 0) stop 'wrong cm number in calling trans_cm_index'

cm_k=int((cm_t-1)/num_cm_zlev)+1
remain1=cm_t-(cm_k-1)*num_cm_zlev
cm_j=int((remain1-1)/num_cmesh(1))+1
cm_i=remain1-(cm_j-1)*num_cmesh(1)

return
end subroutine 

!read in penmsh format prbname.src, only take each point src to
!one fine mesh
subroutine ReadInSrc
use paraset1
use paraset2
use paraset3
use paraset4
use files
use lineread
use ErrControl
use funs

implicit  none
integer i,j,k
integer ink,inj,ini,si,sj,sk
real fi,fj,fk,finvol
integer i_var,IsPart
logical ex
integer mesh(6), cm_start(3), cm_end(3) !, fm_start(3), fm_end(3)
integer::  cen_lft(3)=0, cen_rgt(3)=0 !, cen_src(3)=0
real ::fra_z(2)=1.0,fra_y(2)=1.0,fra_x(2) 
real x,y,z
real ::tot_src_mesh=0, tot_fm_mesh=0,src_fine=0,src_fine_tot=0

flag=0
tot_src_mesh=0
tot_fm_mesh=0

curfile%name=prbname
curfile%ext='src'
curfile%fullname=''
curfile%dir=input_dir
cur_filename=Getfullname(curfile,1)

inquire(file=cur_filename,exist=ex)

IsFile: if(ex  ) then
 allocate (x_pos_src(xnum_s_mesh+1), y_pos_src(ynum_s_mesh+1),&
          z_pos_src(znum_s_mesh+1) )
  
 open(unit=88,file=cur_filename)
!false warning
! warning_message='Reading Src distribution file: '//cur_filename
 write(warning_message,"('Reading Src distribution file: ',A)") trim(cur_filename)
 call DisplayMsg
 cur_var='x_pos_src in .src file'
 if(NumCmtLine(88,flag) .ge. 0) &
  read(88,*,iostat=i_var) (x_pos_src(i),i=1,xnum_s_mesh+1)
 if(i_var .ne. 0) then
   write(warning_message,"('read error, skiped file: ',A)") trim(cur_var)
   call TrapReadError(0)
   return
 endif
 cur_var='y_pos_src in .src file'
 if(NumCmtLine(88,flag) .ge. 0) &
  read(88,*,iostat=i_var) (y_pos_src(i),i=1,ynum_s_mesh+1)
 if(i_var .ne. 0) then
   write(warning_message,"('read error, skiped file: ',A)") trim(cur_var)
   call TrapReadError(0)
   return
 endif
 cur_var='z_pos_src in .src file'
 if(NumCmtLine(88,flag) .ge. 0) &
   read(88,*,iostat=i_var) (z_pos_src(i),i=1,znum_s_mesh+1)
 if(i_var .ne. 0) then
   write(warning_message,"('read error, skiped file: ',A)") trim(cur_var)
   call TrapReadError(0)
   return
 endif
 cur_var='2 comment lines between grid coordinates and src denstisies '
 if(NumCmtLine(88,flag) .ge. 0) then
  read(88,*)
  read(88,*,iostat=i_var)
 endif
 if(i_var .ne. 0) then
   write(warning_message,"('read error, skiped file: ',A)") trim(cur_var)
   call TrapReadError(0)
   return
 endif

 allocate(src_mesh(xnum_s_mesh,ynum_s_mesh,znum_s_mesh))
 
 do k=1, znum_s_mesh
   read(88,*,iostat=i_var) src_mesh(:,:,k)
   if(i_var .ne. 0) then
    write(warning_message,"('read error, skiped file: ',A)") trim(cur_filename)
    call TrapReadError(0)
    return
   endif
  if(k .ne. znum_s_mesh) then
  cur_var='2 comment lines between z levs of src denstisies '
  if(NumCmtLine(88,flag) .ge. 0) then
  read(88,*)
  read(88,*,iostat=i_var)
  endif
  if(i_var .ne. 0) then
   write(warning_message,"('read error, skiped file: ',A)") trim(cur_var)
   call TrapReadError(0)
   return
  endif
  endif
 enddo !endk
else
  return
endif IsFile

close(88)

warning_message='src file done: '// cur_filename
call DisplayMsg

x=x_pos_src(1)
y=y_pos_src(1)
z=z_pos_src(1)

call MapPos2Mesh(x,y,z,mesh,0)
if(mesh(1) .gt. num_zlev .or. mesh(2) .gt. zlevel(1)%ncx &
 .or. mesh(3) .gt. zlevel(1)%ncy) then
write(warning_message,"('src grid is outside the model, skiped src projection: ')") 
call DisplayMsg
return
endif
IsPart=0
cm_start(1)=max(1,mesh(2))
cm_start(2)=max(1,mesh(3))
cm_start(3)=max(1,mesh(1) )
if(cm_start(1) .ne. mesh(2) .or. cm_start(2) .ne. mesh(3) &
  .or. cm_start(3) .ne. mesh(1) )  IsPart=1

x=x_pos_src(xnum_s_mesh+1)
y=y_pos_src(ynum_s_mesh+1)
z=z_pos_src(znum_s_mesh+1)
call MapPos2Mesh(x,y,z,mesh,0)
if(mesh(1) .eq. 0 .or. mesh(2) .eq. 0 &
 .or. mesh(3) .eq. 0) then
write(warning_message,"('src grid is outside the model, skiped src projection: ')") 
call DisplayMsg
return
endif

cm_end(1)=min(mesh(2),zlevel(1)%ncx)
cm_end(2)=min(mesh(3),zlevel(1)%ncy)
cm_end(3)=min(mesh(1),num_zlev)
if(cm_end(1) .ne. mesh(2) .or. cm_end(2) .ne. mesh(3) &
  .or. cm_end(3) .ne. mesh(1) ) IsPart=1

if(IsPart .eq. 1) then
write(warning_message,"('part of src grid is outside the model')") 
call DisplayMsg
endif


allocate(src_mesh_vol(xnum_s_mesh,ynum_s_mesh,znum_s_mesh))
!calculate src grid volume and total
do sk=1, znum_s_mesh

fk=z_pos_src(sk+1)-z_pos_src(sk)
if(fk .le. 0 ) then 
write(warning_message,&
 "('src grid z is not in increasing order: z(',I0,')=',f7.2, '  z(',I0,')=',f7.2)") &
 sk,z_pos_src(sk),sk+1,z_pos_src(sk+1) 
call DisplayMsg
return
endif

do sj=1, ynum_s_mesh
fj=y_pos_src(sj+1)-y_pos_src(sj)
if(fj .le. 0 ) then 
write(warning_message,"('src grid y is not in increasing order: y(',I0,')=',f7.2, '  y(',I0,')=',f7.2)")&
   sj,y_pos_src(sj),sj+1,y_pos_src(sj+1) 
call DisplayMsg
return
endif


do si=1, xnum_s_mesh
fi=x_pos_src(si+1)-x_pos_src(si)
if(fi .le. 0 ) then 
write(warning_message,"('src grid x is not in increasing order: x(',I0,')=',f7.2, '  x(',I0,')=',f7.2)")&
  si,x_pos_src(si),si+1,x_pos_src(si+1) 
call DisplayMsg
return
endif


if(src_mesh(si,sj,sk) .lt. 0 ) then 
write(warning_message,"('src grid value is less than 0, skiped src projection: ','s(',I0,',',I0,',',I0,')=',ES12.5 )") &
 si,sj,sk,src_mesh(si,sj,sk)
call DisplayMsg
return
endif

src_mesh_vol(si,sj,sk)=fk*fj*fi
tot_src_mesh=tot_src_mesh+src_mesh_vol(si,sj,sk)*src_mesh(si,sj,sk)

enddo
enddo
enddo

do k=cm_start(3), cm_end(3)
do j=cm_start(2), cm_end(2)
do i=cm_start(1), cm_end(1)
cen_lft=0
cen_rgt=0
finvol=zlevel(k)%cm_zlev(i,j)%fine_vol

ink_loop: do ink=1,zlevel(k)%cm_zlev(i,j)%fine(3)

call GetPos_k(i,j,k,ink,cen_lft(3),cen_rgt(3),fra_z)
if (cen_rgt(3) .eq. 0 ) cycle ink_loop  

inj_loop: do inj=1, zlevel(k)%cm_zlev(i,j)%fine(2)

call GetPos_j(i,j,k,inj,cen_lft(2),cen_rgt(2),fra_y)
if (cen_rgt(2) .eq. 0 ) cycle inj_loop  

ini_loop: do ini=1,zlevel(k)%cm_zlev(i,j)%fine(1)

call GetPos_i(i,j,k,ini,cen_lft(1),cen_rgt(1),fra_x)
if (cen_rgt(1) .eq. 0 ) cycle ini_loop  

src_fine_tot=0
fi=1

loop_sk: do sk=max(1,cen_lft(3)), cen_rgt(3)
 fk=1
 if(sk .eq. cen_lft(3)) fk=fra_z(1)
 if(sk .eq. cen_rgt(3)) fk=fra_z(2)
 if (fk .eq. 0) cycle loop_sk
loop_sj: do sj=max(1,cen_lft(2)),cen_rgt(2) 
 fj=1
 if(sj .eq. cen_lft(2)) fj=fra_y(1)
 if(sj .eq. cen_rgt(2)) fj=fra_y(2)
 if(fj .eq. 0) cycle loop_sj 
loop_si: do si=max(1,cen_lft(1)),cen_rgt(1)
 fi=1
if(si .eq. cen_lft(1)) fi=fra_x(1)
if(si .eq. cen_rgt(1)) fi=fra_x(2)
if(fi .eq. 0) cycle loop_si
src_fine=fk*fj*fi*src_mesh(si,sj,sk)*src_mesh_vol(si,sj,sk)
src_fine_tot=src_fine_tot+src_fine/finvol
tot_fm_mesh=tot_fm_mesh+src_fine


enddo loop_si !si
enddo loop_sj !sj
enddo loop_sk !sk

if(zlevel(k)%cm_zlev(i,j)%s_flag .eq. 0 .and. &
    src_fine_tot .ne. 0 ) then
zlevel(k)%cm_zlev(i,j)%s_flag=1
num_source=num_source+1
endif
zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)=&
 zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)+src_fine_tot
zlevel(k)%cm_zlev(i,j)%sum_act=zlevel(k)%cm_zlev(i,j)%sum_act+&
	 src_fine_tot
enddo ini_loop
enddo inj_loop
enddo ink_loop
   

enddo !i
enddo !j
enddo !k 

!if(tot_fm_mesh .ne. tot_src_mesh) then
write(*,"('total source magnitude in src grid: ', ES12.5)") &
       tot_src_mesh
write(*,"('total src projected into fine mesh: ', ES12.5)") &
    tot_fm_mesh
write(READLOG,"('total source magnitude in src grid: ', ES12.5)") &
       tot_src_mesh
write(READLOG,"('total src projected into fine mesh: ', ES12.5)") &
    tot_fm_mesh

!endif

deallocate(src_mesh, src_mesh_vol)
deallocate(x_pos_src,y_pos_src,z_pos_src)
return
end subroutine


!read in penmsh format prbname.src, only take each point src to
!one fine mesh
subroutine ReadInFis
use paraset1
use paraset2
use paraset3
use paraset4
use files
use lineread
use ErrControl
use funs

implicit  none
real :: fis_src_mag=0, fis_src_wgt=0
integer cmi,cmj,cmk
real, allocatable :: fis_src_dis(:,:,:)
logical ex

curfile%name='cargosrc'
curfile%ext='dat'
curfile%fullname=''
curfile%dir=input_dir
cur_filename=Getfullname(curfile,1)

inquire(file=cur_filename,exist=ex)

if(.NOT. ex  ) return
  
open(unit=88,file=cur_filename)
 write(warning_message,"('Reading Fission Src distribution file: ',A)") trim(cur_filename)
 call DisplayMsg
 cur_var='srcmag'
read(88,*, end=1002, err=1002)

read(88,*, end=1002, err=1002) fis_src_mag
read(88,*, end=1002, err=1002)
!fixed fission src is the frist src
call TransCmIndex(nscmsh(1),cmi,cmj,cmk)
allocate (fis_src_dis( zlevel(cmk)%cm_zlev(cmi,cmj)%fine(1), &
   zlevel(cmk)%cm_zlev(cmi,cmj)%fine(2), zlevel(cmk)%cm_zlev(cmi,cmj)%fine(3) ))
read(88,*, end=1002, err=1002) fis_src_dis
           
close(88) 

zlevel(cmk)%cm_zlev(cmi,cmj)%src_matrix=&
 zlevel(cmk)%cm_zlev(cmi,cmj)%src_matrix+fis_src_dis
smag(1)=smag(1)+fis_src_mag
zlevel(cmk)%cm_zlev(cmi,cmj)%sum_act=smag(1)
fis_src_wgt=fis_src_mag/smag(1)
serg(1:num_group)=serg(1:num_group)*(1-fis_src_wgt)+chig(1:num_group)*fis_src_wgt

 write(warning_message,"('Induced Fission Src is combined with spontaneous fission src (in src# 1)')") 
 call DisplayMsg

 return
1002  call TrapReadError(1)
end subroutine

!on output
!k1 starting src mesh for fine mesh fm_k 
!k2 endding src mesh for fine mesh fm_k  
!on input
!k1, k2, for last fine mesh
!fk(1): first src mesh fraction in fine mesh fm_k
!fk(2): last srt mesh fraction in fine mesh fm_k  
subroutine GetPos_k(cm_i,cm_j,cm_k,fm_k,k1,k2,fk)
use paraset4
use paraset2

integer cm_i,cm_j,cm_k, fm_k,k1,k2
real fk(2)

real top,bot
! integer start
integer i

!top=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%cm_z(fm_k)+0.5*&
!    zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(3)

!bot=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%cm_z(fm_k)-0.5*&
!    zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(3)

bot=z_lev_pos(cm_k)+(fm_k-1)*zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(3)
top=bot+zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(3)

!fine mesh is outside src grid
if(top .le. z_pos_src(1) .or.  bot .ge. z_pos_src(znum_s_mesh+1)) then
k1=0
k2=0
fk=0
return
endif  
!calculate the starting src mesh k1
if(fm_k .eq. 1) then
 do i=1,znum_s_mesh+1
  if (bot .le. z_pos_src(i)) then
   k1=i-1
   exit
  endif
 enddo
 
else
 k1=k2
endif

!find the fraction fk(1) inside the fine mesh for the first src mesh
if(k1 .gt. 0) then
fk(1)=(min(z_pos_src(k1+1),top)-bot)/(z_pos_src(k1+1)-z_pos_src(k1))
else
!first src mesh is inside current fine mesh
!if not, fk(2) will be used
fk(1)=1
endif

!setup k2 and fk(2)
top_loop: do i=k1+1, znum_s_mesh+1

 if(z_pos_src(i) .gt. top) then
  k2=i-1
  fk(2)=(top-max(bot,z_pos_src(k2)) )/(z_pos_src(i)-z_pos_src(k2))
  exit top_loop
 endif
enddo top_loop
!last src mesh top is under fine mesh top
if(i .eq. znum_s_mesh+2) then
 k2=znum_s_mesh
 fk(2)=(z_pos_src(k2+1)-max(bot,z_pos_src(k2)) )/(z_pos_src(k2+1)-z_pos_src(k2))
endif

return

end subroutine


subroutine GetPos_j(cm_i,cm_j,cm_k,fm_k,k1,k2,fk)
use paraset4
use paraset2

integer cm_i,cm_j,cm_k, fm_k,k1,k2
real fk(2)

real top,bot
! integer start
integer i
!top=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%cm_y(fm_k)+0.5*&
!    zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(2)

!bot=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%cm_y(fm_k)-0.5*&
!    zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(2)
bot=zlevel(1)%y_cm_pos(cm_j)+(fm_k-1)*zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(2)
top=bot+zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(2)
!fine mesh is outside src grid
if(top .le. y_pos_src(1) .or.  bot .ge. y_pos_src(ynum_s_mesh+1)) then
k1=0
k2=0
fk=0
return
endif  


if(fm_k .eq. 1) then	
 do i=1,ynum_s_mesh+1
  if (bot .le. y_pos_src(i)) then
   k1=i-1
   exit
  endif
 enddo
 
else
 k1=k2
endif

if(k1 .gt. 0) then
fk(1)=(min(y_pos_src(k1+1),top)-bot)/(y_pos_src(k1+1)-y_pos_src(k1))
else
fk(1)=1
endif


top_loop: do i=k1+1, ynum_s_mesh+1

 if(y_pos_src(i) .gt. top) then
  k2=i-1
  fk(2)=(top-max(bot,y_pos_src(k2)))/(y_pos_src(i)-y_pos_src(k2))
  exit top_loop
 endif
enddo top_loop

if(i .eq. ynum_s_mesh+2) then
 k2=ynum_s_mesh
 fk(2)=(y_pos_src(k2+1)-max(bot,y_pos_src(k2)) )/(y_pos_src(k2+1)-y_pos_src(k2))
endif


return

end subroutine

subroutine GetPos_i(cm_i,cm_j,cm_k,fm_k,k1,k2,fk)
use paraset4
use paraset2

integer cm_i,cm_j,cm_k, fm_k,k1,k2
real fk(2)

real top,bot
!integer start
integer i
!top=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%cm_x(fm_k)+0.5*&
!    zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(1)

!bot=zlevel(cm_k)%cm_zlev(cm_i,cm_j)%cm_x(fm_k)-0.5*&
!    zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(1)

bot=zlevel(1)%x_cm_pos(cm_i)+(fm_k-1)*zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(1)
top=bot+zlevel(cm_k)%cm_zlev(cm_i,cm_j)%finsize(1)

!fine mesh is outside src grid
if(top .le. x_pos_src(1) .or.  bot .ge. x_pos_src(xnum_s_mesh+1)) then
k1=0
k2=0
fk=0
return
endif  

if(fm_k .eq. 1) then
 do i=1,xnum_s_mesh+1
  if (bot .le. x_pos_src(i)) then 
   k1=i-1
   exit
  endif
 enddo
  
else
 k1=k2
endif

if(k1 .gt. 0) then
fk(1)=(min(x_pos_src(k1+1),top)-bot)/(x_pos_src(k1+1)-x_pos_src(k1))
else
fk(1)=1
endif


top_loop: do i=k1+1, xnum_s_mesh+1

 if(x_pos_src(i) .gt. top) then
  k2=i-1
  fk(2)=(top-max(bot,x_pos_src(k2)))/(x_pos_src(i)-x_pos_src(k2))
  exit top_loop
 endif
enddo top_loop

if(i .eq. xnum_s_mesh+2) then
 k2=xnum_s_mesh
 fk(2)=(x_pos_src(k2+1)-max(bot,x_pos_src(k2)) )/(x_pos_src(k2+1)-x_pos_src(k2))
endif


return

end subroutine
