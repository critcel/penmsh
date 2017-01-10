!subcode1.f90 read input files
!subroutine read penmsh input file
!Subroutines:
!            ReadInput:  Read inputfiles
!            TrapReadError: Handle read errors
subroutine ReadInput
use files
use funs
use ErrControl
use lineread
use paraset1
use paraset2
use paraset4
use mIn4deck
integer i,j,k,m, IsAllCmat
integer bit_num_zlev,cm_num_1z
integer cmx,cmy
integer overlay,overlay_typ,overlay_num_abs, over_num
integer name_len,dir_len,ext_len,fullname_len
! integer err_open
! integer err_read

!overlay block
integer block, num_block

!for cm type3
integer cm_typ, num_reg

!for repeat overlay block
integer :: rep_ijk(3), rep_i, rep_j, rep_k
integer :: temp_int
integer :: mycm_type, mycm_x, mycm_y
real    :: cmtyp
logical :: ex

penmsh_inp%dir=input_dir
penmsh_inp%fullname='penmsh.inp'

write(READLOG,*) 'opening penmsh.inp'

dir_len=len_trim(penmsh_inp%dir)
name_len=len_trim(penmsh_inp%fullname)
write(form,"( '(','A',I0,',A1,','A',I0,')' )" ) dir_len, name_len
write(cur_filename,fmt=form) penmsh_inp%dir ,'/',penmsh_inp%fullname
cur_fileunit=31
open(UNIT=cur_fileunit, file=trim(cur_filename),STATUS='OLD',err=1001)
write(READLOG,*) 'penmsh.inp open success'

cur_var='problem name'
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit, '(A)',err=1002,end=1002) prbname
endif

write(READLOG,*) 'problem name : ', trim(prbname)
write(READLOG,*)

prbname_len=len_trim(prbname)

cur_var='num_zlev,num_material,flag_mathematica'
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit, *,err=1002,end=1002) num_zlev,num_material,flag_mathematica
endif

write(READLOG,*) 'num of z level:  ', num_zlev
write(READLOG,*) 'num of material: ', num_material
write(READLOG,*) 'mathematica:     ', flag_mathematica
write(READLOG,*) 

if (num_zlev .le. 0 .or. num_material .le. 0) then

 write(READLOG,*) 'num of z level or num of material 0 detected' 
  write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1000)'
  call TrapInputError(1)

endif

allocate(z_lev_pos(num_zlev+1),max_zfine(num_zlev))
allocate(ratio_x(num_zlev),ratio_y(num_zlev),ratio_z(num_zlev))

cur_var='z lev position'
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit, *,err=1002,end=1002) (z_lev_pos(j),j=1,num_zlev+1)
endif

call check_real(z_lev_pos, num_zlev+1)

write(READLOG, *) 'z level boundary: ', z_lev_pos
write(READLOG,*) 

cur_var='max zfine for each z lev'
max_zfine=0
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit,*,err=1002,end=1002) max_zfine
endif

write(READLOG,*) 'max z fine mesh for each z level:'
write(READLOG,*) max_zfine	
write(READLOG,*)

m=0
do k=1, num_zlev
if(max_zfine(k) .le. 0) then
  m=-1
  exit
endif

enddo
if(m .eq. -1) then
write(warning_message,"('Warning 4001:', A, ' <=0 in file: ',A)")&
          trim(cur_var), trim(cur_filename)
call DisplayMsg
endif

cur_var='fine-to-med mesh ratio along x '
ratio_x=0
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit, *,err=1002,end=1002) ratio_x
endif

write(READLOG,*) 'fine to med ratio along x for each z level:'
write(READLOG,*) ratio_x

m=0
do k=1, num_zlev
if(ratio_x(k) .le. 0) then
 ratio_x(k)=1
 m=-1
endif 
enddo
if(m .eq. -1) then
 write(warning_message,"('Warning 4002:', A, ' <=0 in file: ' )")&
          trim(cur_var), trim(cur_filename)
  call DisplayMsg
endif

cur_var='fine-to-med mesh ratio along y '
ratio_y=0
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  read(cur_fileunit, *,err=1002,end=1002) ratio_y
endif

write(READLOG,*) 'fine to med ratio along y for each z level:'
write(READLOG,*) ratio_y

m=0
do k=1, num_zlev
if(ratio_y(k) .le. 0) then
 ratio_y(k)=1
 m=-1
endif 
enddo
if(m .eq. -1) then
 write(warning_message,"('Warning 4003:', A, ' <=0 in file: ' )")&
          trim(cur_var), trim(cur_filename)
  call DisplayMsg
endif

cur_var='fine-to-med mesh ratio along z '
ratio_z=0
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  read(cur_fileunit, *,err=1002,end=1002) ratio_z
endif

write(READLOG,*) 'fine to med ratio along z for each z level:'
write(READLOG,*) ratio_z

m=0
do k=1, num_zlev
if(ratio_y(k) .le. 0) then
 ratio_y(k)=1
 m=-1
endif 
enddo
if(m .eq. -1) then
 write(warning_message,"('Warning 4004:', A, ' <=0 in file: ' )")&
          trim(cur_var), trim(cur_filename)
  call DisplayMsg
endif

write(READLOG,*) 

cur_var='source_format, xnum_s_mesh ... '
num_group=-1
sn_order=-1
pn_order=-1
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit, *,err=1002,end=1002) s_format, xnum_s_mesh, ynum_s_mesh, znum_s_mesh,num_group, sn_order,pn_order
endif

write(READLOG,*) 'source format: ', s_format
write(READLOG,*) 'x mesh in src grid(.src): ', xnum_s_mesh
write(READLOG,*) 'y mesh in src grid(.src): ', ynum_s_mesh
write(READLOG,*) 'z mesh in src grid(.src): ', znum_s_mesh
write(READLOG,*) 'num of group: '         , num_group
write(READLOG,*) 'Sn order: ', sn_order
write(READLOG,*) 'Pn order: ', pn_order

if(num_group .lt. 0 .or. sn_order .lt. 0 .or. pn_order .lt. 0) then
  write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1001)'
	call TrapInputError(1)
endif

write(READLOG,*)

num_comment=-1
leg_order=-1
len_table=-1
cur_var='xs_format,num_comment...'
num_comment=-1
leg_order=-1
len_table=-1
if(NumCmtLine(cur_fileunit,flag) .ge. 0)  then
 read(cur_fileunit, *,err=1002,end=1002) xs_format,num_comment,leg_order,len_table
endif 

write(READLOG, *) 'cross section information:'
write(READLOG, *) 'xs format:     ', xs_format
write(READLOG, *) 'num_comment:   ', num_comment
write(READLOG, *) 'legdre order:  ', leg_order
write(READLOG, *) 'table length:  ', len_table
write(READLOG, *)

if(num_comment .lt. 0 .or. leg_order .lt. 0 .or. len_table .lt. 0) then
  write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1002)'
  call TrapInputError(1)
endif

cur_var='boundary cond. '
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
! ibback, ibfrnt, jbright,jbleft,kbottom,ktop
 read(cur_fileunit, *,err=1002,end=1002) bcs
endif

write(READLOG, *) 'boundary conditions:'
write(READLOG,*) 'ibback=   ', bcs(1)
write(READLOG,*) 'ibfrnt=   ', bcs(2)
write(READLOG,*) 'jbright=  ', bcs(3)
write(READLOG,*) 'jbleft=   ', bcs(4)
write(READLOG,*) 'kbottom=  ', bcs(5)
write(READLOG,*) 'ktop=     ', bcs(6)
write(READLOG,*) 

select case(s_format)
case(:-1)
 num_src_mat=abs(s_format)
 allocate(id_s_mat(num_src_mat),s_intensity(num_src_mat))
 cur_var='source mat ID'
 id_s_mat=0
 if(NumCmtLine(cur_fileunit,flag) .ge. 0)  read(cur_fileunit,*,err=1002,end=1002) id_s_mat
 write(READLOG,*) 'source material: ', id_s_mat
 
 do k=1, num_src_mat
 if(id_s_mat(k) .le. 0 .or. id_s_mat(k).gt. num_material) then
  write(err_message,*) trim(cur_var)//' : input value error  (ERR1003)'
  call TrapInputError(1)
 endif
 enddo

 cur_var='source intensity'
 if(NumCmtLine(cur_fileunit,flag) .ge. 0)  read(cur_fileunit,*,err=1002,end=1002) s_intensity
 write(READLOG,*) 'source intensity(neg value as total magnitude): ', s_intensity
 
 allocate(s_mag_src(num_material))
 s_mag_src=0
 do i=1,num_src_mat
  s_mag_src(id_s_mat(i))=s_intensity(i)
 enddo


case(1)
 write(warning_message,*) 's_format=1: specified for venus source only in PENMSH, ignored'
 s_format=0
 call TrapError(-1)
case(2)
 write(warning_message,*) 's_format=2: no src in the model, .src file ignored even if it exists'
 call TrapError(-1)
case default
 s_format=0
end select  

close(cur_fileunit)

write(READLOG, *) 'penmsh.inp read complete'
write(READLOG, *) '********************************************************'
write(READLOG,*)


allocate (zlevel(num_zlev))

prb_inp%name=prbname
prb_inp%dir=input_dir
prb_inp%ext='inp'
prb_inp%fullname=''

prb_inp%fullname=Getfullname(prb_inp, 1)
cur_filename=prb_inp%fullname
inquire(file=cur_filename,exist=ex)
if (ex) then

write(warning_message,"('Warning:  all in one (all z-lev combined) input file found:', A )") trim(cur_filename)
call DisplayMsg

call ReadPrbInp

write(warning_message,"(A, ' read successfully, individual z-lev inp files ignored' )") trim(cur_filename)
call DisplayMsg

z_start=1
z_end=num_zlev
allocate(cgdim(1)%boundary(0:zlevel(1)%ncx),cgdim(2)%boundary(0:zlevel(1)%ncy),&
         cgdim(3)%boundary(0:num_zlev) )


cgdim(1)%boundary=zlevel(1)%x_cm_pos
cgdim(2)%boundary=zlevel(1)%y_cm_pos
cgdim(3)%boundary=z_lev_pos


if(IsZSingle .eq. 0) call ReadIn4
! construct cross section file name prbname.xs
prb_xs%name=prbname
prb_xs%dir=''
prb_xs%ext='xs'
prb_xs%fullname=''
prb_xs%fullname=Getfullname(prb_xs, 0)
return
endif


if( IsZSingle .eq. 0) then
z_start=1
z_end=num_zlev
endif

allocate (inputfile(num_zlev))

do i=z_start, z_end
 call intlen(i,bit_num_zlev)
 write(form,"('(A',I0,',I0,A4)' )") prbname_len
 write(inputfile(i)%fullname, form ) trim(prbname),i,".inp"
 write(READLOG,*)
 write(READLOG,*) 'opening inputfile: ', trim(inputfile(i)%fullname)

 dir_len=len_trim(input_dir)
 name_len=len_trim(inputfile(i)%fullname)
 write(form,"('(A',I0,',A1,A',I0,')' )" ) dir_len, name_len
 write(cur_filename,form) input_dir,'/',inputfile(i)%fullname
 cur_fileunit=i+50
 open(cur_fileunit, file=cur_filename,STATUS='OLD',ERR=1001)
 write(READLOG,*) trim(inputfile(i)%fullname), ' :  open success'
 write(READLOG,*)
 
 write(READLOG,*) 'information for z level: ', i
 
 cur_var='ncx,ncy... '
 
 zlevel(i)%ncx=0
 zlevel(i)%ncy=0
zlevel(i)%num_zfine=0
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
    read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%ncx,zlevel(i)%ncy,zlevel(i)%num_zfine  
 endif 
 
 write(READLOG,*) 'number of coarse mesh along x(ncx):   ', zlevel(i)%ncx
 write(READLOG,*) 'number of coarse mesh along x(ncy):   ', zlevel(i)%ncy
 write(READLOG,*) 'max z fine mesh for this z level:     ', zlevel(i)%num_zfine
 write(READLOG,*)

if(zlevel(i)%ncx .le. 0 .or.  zlevel(i)%ncy .le. 0 ) then
  write(err_message,*) trim(cur_var)//'less than zero:  input value error or read failure (ERR1004)'
  call TrapInputError(1)
endif

allocate(zlevel(i)%cm_zlev(zlevel(i)%ncx,zlevel(i)%ncy)) 

num_need=zlevel(i)%ncy*zlevel(i)%ncx
allocate(fidobank(num_need))

cur_var='x fine mesh for each cm'

num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo
m=0
do k=1,zlevel(i)%ncy
do j=1,zlevel(i)%ncx
m=m+1
zlevel(i)%cm_zlev(j,k)%fine(1)=nint(fidobank(m))
enddo
enddo

 
 write(READLOG,*) 'num of x fine mesh for each coarse mesh in this z level: '
  
 do k=1, zlevel(i)%ncy
   write(READLOG,*) (zlevel(i)%cm_zlev(j,k)%fine(1), j=1,zlevel(i)%ncx)
 enddo
 write(READLOG,*)

do k=1, zlevel(i)%ncy
 do j=1, zlevel(i)%ncx
   if(zlevel(i)%cm_zlev(j,k)%fine(1) .le. 0) then
    write(err_message,*) trim(cur_var)//': input value error or read failure (ERR1005)'
    call TrapInputError(1)
   endif
 enddo
enddo

 cur_var='y fine mesh for each cm'
 
 num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo
m=0
do k=1,zlevel(i)%ncy
do j=1,zlevel(i)%ncx
m=m+1
zlevel(i)%cm_zlev(j,k)%fine(2)=nint(fidobank(m))
enddo
enddo
 
 write(READLOG,*) 'num of y fine mesh for each coarse mesh in this z level: '
 do k=1, zlevel(i)%ncy
   write(READLOG,*) (zlevel(i)%cm_zlev(j,k)%fine(2), j=1,zlevel(i)%ncx)
 enddo
 write(READLOG,*)

do k=1, zlevel(i)%ncy
 do j=1, zlevel(i)%ncx
   if(zlevel(i)%cm_zlev(j,k)%fine(2) .le. 0) then
    write(err_message,*) trim(cur_var)//': input value error or read failure (ERR1006)'
    call TrapInputError(1)
   endif
 enddo
enddo

if(zlevel(i)%num_zfine .lt. 0) then
  cur_var='z fine mesh for each cm'
  
num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo
m=0
do k=1,zlevel(i)%ncy
do j=1,zlevel(i)%ncx
m=m+1
zlevel(i)%cm_zlev(j,k)%fine(3)=nint(fidobank(m))
enddo
enddo

elseif(zlevel(i)%num_zfine .gt. 0) then
 
  do k=1,zlevel(i)%ncy
   do j=1,zlevel(i)%ncx

    zlevel(i)%cm_zlev(j,k)%fine(3)=zlevel(i)%num_zfine
   
   enddo
  enddo
 
 else
    err_message='max z fine number =0' 
    call TrapInputError(1)
 endif
 
 write(READLOG,*) 'num of z fine mesh for each coarse mesh in this z level: '
 do k=1, zlevel(i)%ncy
   write(READLOG,*) (zlevel(i)%cm_zlev(j,k)%fine(3), j=1,zlevel(i)%ncx)
 enddo
 write(READLOG,*)

do k=1, zlevel(i)%ncy
 do j=1, zlevel(i)%ncx
   if(zlevel(i)%cm_zlev(j,k)%fine(3) .le. 0) then
    write(err_message,*) trim(cur_var)//': input value error or read failure (ERR1007)'
    call TrapInputError(1)
   endif
 enddo
enddo
 
 allocate(zlevel(i)%x_cm_pos(1+zlevel(i)%ncx))
 allocate(zlevel(i)%y_cm_pos(1+zlevel(i)%ncy))
 
 cur_var='x_cm_pos '
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%x_cm_pos
 endif
 
 call check_real(zlevel(i)%x_cm_pos,1+zlevel(i)%ncx)

 write(READLOG,*) 'coarse mesh position along x: '
 write(READLOG,*) zlevel(i)%x_cm_pos
 write(READLOG,*)

 cur_var='y_cm_pos'
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%y_cm_pos
 endif

call check_real(zlevel(i)%y_cm_pos,1+zlevel(i)%ncy)

 write(READLOG,*) 'coarse mesh position along y: '
 write(READLOG,*) zlevel(i)%y_cm_pos
 write(READLOG,*)

! re-size mesh
if (max_fm_size.gt.0) then
 do k=1,zlevel(i)%ncy
  do j=1,zlevel(i)%ncx
   zlevel(i)%cm_zlev(j,k)%fine(1)=ceiling((zlevel(i)%x_cm_pos(j+1)-zlevel(i)%x_cm_pos(j))/max_fm_size)*zlevel(i)%cm_zlev(j,k)%fine(1)
   zlevel(i)%cm_zlev(j,k)%fine(2)=ceiling((zlevel(i)%y_cm_pos(k+1)-zlevel(i)%y_cm_pos(k))/max_fm_size)*zlevel(i)%cm_zlev(j,k)%fine(2)
   zlevel(i)%cm_zlev(j,k)%fine(3)=ceiling((z_lev_pos(i+1)-z_lev_pos(i))/max_fm_size)*zlevel(i)%cm_zlev(j,k)%fine(3)
   !write(*,'(6I4)') j,k,i ,zlevel(i)%cm_zlev(j,k)%fine(1:3)
  enddo
 enddo
endif

IsAllCmat=1
cur_var='cm_type '
num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

m=0
do k=1,zlevel(i)%ncy
do j=1,zlevel(i)%ncx
m=m+1
cmtyp=fidobank(m)/10.0-int(fidobank(m)/10.0)
cmtyp=10*cmtyp
zlevel(i)%cm_zlev(j,k)%cm_type=nint(fidobank(m))
if(abs(nint(cmtyp)) .gt. 1 .and. IsAllCmat .eq. 1) IsAllCmat=0

enddo
enddo

 
 write(READLOG,*) 'coarse mesh type in this z level: '
 do k=1, zlevel(i)%ncy
   write(READLOG,*) (zlevel(i)%cm_zlev(j,k)%cm_type, j=1,zlevel(i)%ncx)
 enddo
 write(READLOG,*)

cur_var='num_subregion '

num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

m=0
do k=1,zlevel(i)%ncy
do j=1,zlevel(i)%ncx
m=m+1
zlevel(i)%cm_zlev(j,k)%num_subregion=nint(fidobank(m))
if(abs(zlevel(i)%cm_zlev(j,k)%num_subregion) .ne. 1 .and. IsAllCmat .eq. 1) IsAllCmat=0
enddo
enddo
deallocate(fidobank)

do k=1, zlevel(i)%ncy
do j=1, zlevel(i)%ncx
   if(zlevel(i)%cm_zlev(j,k)%num_subregion .eq. 0) then
    write(err_message,*) trim(cur_var)//': input value error or read failure (ERR1008)'
    call TrapInputError(1)
   endif
 enddo
enddo
 

 write(READLOG,*) 'number of subregions in corase mesh  in this z level: '
 do k=1, zlevel(i)%ncy
   write(READLOG,*) (zlevel(i)%cm_zlev(j,k)%num_subregion, j=1,zlevel(i)%ncx)
 enddo
 write(READLOG,*)

!read cm mat all in once if cm_type=1 and num_region=1 uniformly
if(IsAllCmat .eq. 1) then

cur_var='cm_mat '
num_need=zlevel(i)%ncy*zlevel(i)%ncx
allocate(fidobank(num_need))
num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

m=0
do cmy=1,zlevel(i)%ncy
do cmx=1,zlevel(i)%ncx
m=m+1
cm_num_1z=cmx+(cmy-1)*zlevel(i)%ncx
zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num=nint(fidobank(m))
if(zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num .gt. num_material .or. zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num .le. 0) then
   write(err_message, "('Material # > num of materials at cm ',I0,' zlev ',I0)") cm_num_1z, i
   call TrapInputError(1)
endif 
write(READLOG, "(' Mat num for coarse mesh ', I0, ' :  ', I0)" )  &
   cm_num_1z, zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num
enddo
enddo

deallocate(fidobank)

!tranditional read
else

!read in cm boundary
do cmy=1, zlevel(i)%ncy
 do cmx=1, zlevel(i)%ncx
   cm_num_1z=cmx+(cmy-1)*zlevel(i)%ncx
write(cur_var,"('cm',I0,' mat_num/boundary')")  cm_num_1z

cm_typ=mod(abs(zlevel(i)%cm_zlev(cmx,cmy)%cm_type),10)
num_reg=abs(zlevel(i)%cm_zlev(cmx,cmy)%num_subregion)

write(READLOG,"('Coarse mesh ',I0,' : type ',I0, ' with ',I0, ' region(s)')") &
  cm_num_1z, cm_typ, num_reg   

if(num_reg .ne. 1) then
     
allocate(zlevel(i)%cm_zlev(cmx,cmy)%r_lft(num_reg),&
   zlevel(i)%cm_zlev(cmx,cmy)%r_rgt(num_reg),&
   zlevel(i)%cm_zlev(cmx,cmy)%theta_bot(num_reg),&
   zlevel(i)%cm_zlev(cmx,cmy)%theta_top(num_reg),&
   zlevel(i)%cm_zlev(cmx,cmy)%cm_mat(num_reg) )

num_need=num_reg
allocate(fidobank(num_need))
num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(i)%cm_zlev(cmx,cmy)%r_lft(:)=fidobank(:)

write(form, "('(A26,A3,',I0,'(ES12.5,2x))'  )") num_reg
write(READLOG,form ) ' subregion left boundary  ', ' : ', &
    zlevel(i)%cm_zlev(cmx,cmy)%r_lft


num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(i)%cm_zlev(cmx,cmy)%r_rgt(:)=fidobank(:)

write(READLOG,form ) ' subregion right boundary ',' : ', &
    zlevel(i)%cm_zlev(cmx,cmy)%r_rgt

num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(i)%cm_zlev(cmx,cmy)%theta_bot(:)=fidobank(:)

write(READLOG,form ) ' subregion bottom boundary',' : ', &
    zlevel(i)%cm_zlev(cmx,cmy)%theta_bot

num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(i)%cm_zlev(cmx,cmy)%theta_top(:)=fidobank(:)


write(READLOG,form ) ' subregion top boundary   ',' : ', &
   zlevel(i)%cm_zlev(cmx,cmy)%theta_top

num_read=0
lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(i)%cm_zlev(cmx,cmy)%cm_mat(:)=fidobank(:)

write(form,"('(A25,I0,A3,',I0,'(I0,2x))' )" ) num_reg
write(READLOG, form)' Mat num for coarse mesh ', cm_num_1z,':  ', &
   zlevel(i)%cm_zlev(cmx,cmy)%cm_mat

if(maxval(zlevel(i)%cm_zlev(cmx,cmy)%cm_mat) .gt. num_material .or. minval(zlevel(i)%cm_zlev(cmx,cmy)%cm_mat) .lt. 0) then
  write(err_message, "('Material # > num of materials at cm ',I0,' zlev ',I0)") cm_num_1z, i
   call TrapInputError(1)
endif   
 zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num=zlevel(i)%cm_zlev(cmx,cmy)%cm_mat(1)

deallocate(fidobank)

else  !only one region
 
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num
endif
   
if(cmy .eq. zlevel(i)%ncy .and. cmx .eq. zlevel(i)%ncx) flag(2)=1
  

if(NumCmtLine(cur_fileunit,flag) .eq. 0) then
 flag(2)=0
 if(NumCmtLine(cur_fileunit,flag) .ge. 0)  read(cur_fileunit,*,err=1002,end=1002)
 if(NumCmtLine(cur_fileunit,flag) .ge. 0)  read(cur_fileunit,*,err=1002,end=1002)
 if(NumCmtLine(cur_fileunit,flag) .ge. 0)  read(cur_fileunit,*,err=1002,end=1002)
	
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
   read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num
 endif

endif
  
flag(2)=0
if(zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num .gt. num_material .or. zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num .le. 0) then
   write(err_message, "('Material # > num of materials at cm ',I0,' zlev ',I0)") cm_num_1z, i
   call TrapInputError(1)
endif 

write(form,"('(A25,I0,A3,',I0,'(I0,2x))'  )") num_reg
write(READLOG,form ) ' Mat num for coarse mesh ' &
   , cm_num_1z,':  ', zlevel(i)%cm_zlev(cmx,cmy)%cm_mat_num
   
   
endif !num_reg
    
enddo
enddo

endif !IsAllCmat


write(READLOG,*)
write(READLOG,*)

!read in overlay 
do cmy=1, zlevel(i)%ncy
do cmx=1, zlevel(i)%ncx
mycm_type=zlevel(i)%cm_zlev(cmx,cmy)%cm_type
if(mycm_type .gt. 0 .and. mycm_type .lt. 10) cycle
cm_num_1z=cmx+(cmy-1)*zlevel(i)%ncx

!Coarse meshes share same overlay structure
! if(zlevel(i)%cm_zlev(cmx,cmy)%cm_type .gt. 10) then
 !last digit is the cm type, the others is the number of a previous CM, as which current CM has the same overlay stucture 
 !e.g. cm_type=43, means cm type 3, and it has the same overlay structure as CM 4 
if(mycm_type .gt. 10) then
mycm_num=int(zlevel(i)%cm_zlev(cmx,cmy)%cm_type/10)
mycm_x=mod(mycm_num-1,zlevel(i)%ncx)+1
mycm_y=int((mycm_num-1)/zlevel(i)%ncx)+1

if(mycm_num .gt. cm_num_1z) then
 write(err_message, "('shared overlay structure has to appear in previous CM, error in cm_type for CM:',I0 )") cm_num_1z
 call TrapInputError(1)
elseif(zlevel(i)%cm_zlev(mycm_x,mycm_y)%cm_type .ge. 0 .and. zlevel(i)%cm_zlev(mycm_x,mycm_y)%cm_type .le. 10) then
 write(err_message, "('shared overlay structure previous CM is not overlayed, error in cm_type for CM:',I0 )") cm_num_1z
 call TrapInputError(1)
else
  write(READLOG,"('Overlay information for coarse mesh: ',I0,1x)") &
       cm_num_1z
  write(READLOG,"(' overlay structure same as: ',I0,1x)") &
       mycm_num	   	
endif

num_block=zlevel(i)%cm_zlev(mycm_x,mycm_y)%num_block
allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(num_block))


do block=1,num_block
!handle repeat structure
  zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num=&
   zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%overlay_num
!  write(READLOG,*)
!  write(cur_var,"('num of overlay for cm:  ',I0)")  cm_num_1z
!  if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
!   read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num
!  endif
  over_num=zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%overlay_num
  overlay_num_abs=abs(over_num)
  
   
allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    overlay_typ(overlay_num_abs))
allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
       overlay_mat(overlay_num_abs))

allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   num_rep(3,overlay_num_abs), &
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   dist_rep(3,overlay_num_abs) )
     
zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep=1
zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep=0.0
	

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(:)= &
   zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%overlay_typ(:)


!zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j)=&
!	zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%overlay_typ(j)
zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(:)= &
   zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%overlay_mat(:)


allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos, &
   overlay_num_abs))

allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay_num_abs))

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(:,:)=&
     zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%overlay_bon(:,:)
	  
do overlay=1,overlay_num_abs

if(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay) .lt. 0) then


zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(:,overlay)=&
  zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%num_rep(:,overlay)

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep(:,overlay)=&
 zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%dist_rep(:,overlay)
endif !repeat structure

!different mat number for lattice
if(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .lt. 0) then

rep_ijk(:)=zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(:,overlay)
allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay)%d3_int(rep_ijk(1), rep_ijk(2),rep_ijk(3) ) )

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay)%d3_int(:,:,:)= &
   zlevel(i)%cm_zlev(mycm_x,mycm_y)%overlay_block(block)%mat_rep(overlay)%d3_int(:,:,:)


endif  !diff. mat number

enddo !overlay number
enddo !enddo block

cycle
endif 
 

! if(zlevel(i)%cm_zlev(cmx,cmy)%cm_type .le. -1) then
 !last digit is the cm type, the others is the number of overlay block
 !e.g. cm_type=-43, means cm type 3, and 4 overlay blocks 
 num_block=int(-zlevel(i)%cm_zlev(cmx,cmy)%cm_type/10)
 if(num_block .eq. 0) num_block=1
 zlevel(i)%cm_zlev(cmx,cmy)%num_block=num_block
 allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(num_block))
 write(READLOG,*)
 write(READLOG,"('Overlay information for coarse mesh: ',I0,1x)") &
       cm_num_1z	
 write(READLOG,"(1x,I0,' block(s) ')") num_block	     
 do block=1,num_block
!handle repeat structure
  zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num=0
  write(READLOG,*)
  write(cur_var,"('num of overlay for cm:  ',I0)")  cm_num_1z
  if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
   read(cur_fileunit,*,err=1002,end=1002) zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num
  endif
  over_num=zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num
  overlay_num_abs=abs(over_num)
  
  write(READLOG,"('overlay block:  ',I0)") block
  write(READLOG,"(' num of overlay   : ', I0)")  &
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num
   
  if(overlay_num_abs .eq. 0)then
    write(err_message, "('number of overlay should not be zero on z-level ',I0, ' CM ', I0)") i, cm_num_1z
    call TrapInputError(1)
  endif
   
  allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    overlay_typ(overlay_num_abs))
  allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
       overlay_mat(overlay_num_abs))

  allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    num_rep(3,overlay_num_abs), &
    zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    dist_rep(3,overlay_num_abs) )
     
  zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep=1
  zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep=0.0
	
  write(cur_var,"('type of overlay for cm:  ',I0)")  cm_num_1z	

num_need=overlay_num_abs
allocate(fidobank(num_need))

num_read=0

lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   overlay_typ(:)=fidobank(:)

deallocate(fidobank)
!  flag(1)=-1
!  if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
!   if(flag(1) .lt. overlay_num_abs) then
!    write(err_message, &
!     "(' ERR2000: only ',I0,' values found, ',I0,' expected for overlay type at cm ', I0,' on zlev ',I0)") &
!     flag(1),overlay_num_abs, cm_num_1z,i
!    call TrapInputError(1)
!   endif
!  read(cur_fileunit,*,err=1002,end=1002) (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j), &
!     j=1, overlay_num_abs)
!  endif
  if(over_num .lt. 0) then
  do j=1, overlay_num_abs
  temp_int=zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j)
  zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j)=-abs(temp_int)
  enddo
  endif
  
  write(form,"('(A20,',I0,'(I0,2x))' )") overlay_num_abs 	  

  write(READLOG,form) ' overlay type     : ', &
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ
  
  write(cur_var,"('mat num of overlay for cm:  ',I0)")  cm_num_1z
num_need=overlay_num_abs
allocate(fidobank(num_need))

num_read=0

lineseg=''
do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido
enddo
zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   overlay_mat(:)=fidobank(:)

deallocate(fidobank)				
!if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
!if(flag(1) .lt. overlay_num_abs) then
!  write(err_message, &
!   "('ERR2001: only ',I0,' values found, ',I0,' expected for overlay mat at cm ',I0,' on zlev ',I0)")&
!    flag(1),overlay_num_abs, cm_num_1z,i
!  call TrapInputError(1)
!endif
!read(cur_fileunit,*,err=1002,end=1002) (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(j), &
!    j=1, overlay_num_abs)
!endif

do j=1, overlay_num_abs
  temp_int=zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(j)
  if(abs(temp_int) .gt. num_material .or. temp_int .eq. 0 ) then
  write(err_message, "('Overlay Material # > num of materials at cm=',I0,' zlev=',I0, 'block=', I0, 'overlay=', I0)") & 
     cm_num_1z, i, block, j
   call TrapInputError(1)
  endif
  if(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j) .gt. 0) &
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(j)=abs(temp_int)
enddo

 	
write(form,"('(A20,',I0,'(I0,2x))' )") overlay_num_abs 	  
write(READLOG,form) ' overlay mat      : ', &
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat

allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos, &
   overlay_num_abs))

allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay_num_abs))

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos,&
      overlay_num_abs)=0
	  
do overlay=1,overlay_num_abs
 if (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .gt. num_material) then  
  write(err_message, "('ERR2002: Overlay Mat # > num of materials at cm=',I0,' zlev=',I0, 'block=', I0, 'overlay=', I0)") & 
     cm_num_1z, i, block, overlay
  call TrapInputError(1)
endif

if(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .lt. 0 .and. &
 zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay) .gt. 0) then
write(err_message, "('ERR2003: Overlay Mat # < 0 ',I0,' zlev ',I0,' overlay ',I0)") &
  cm_num_1z, i, overlay
 call TrapInputError(1)
endif
 overlay_typ=zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay)
 flag(1)=-1
 write(cur_var,"('overlay boundaries of overlay ',I0,' for cm: ',I0)")  overlay,cm_num_1z

 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos,overlay)=flag(1)
   if(flag(1) .lt. map_overlay_typ2bon(abs(overlay_typ),cm_num_1z)) then
    write(err_message,"('ERR2003: not enough entries for boundaries of overlay ',I0, ' in cm ',I3)")&
        overlay,cm_num_1z
    call TrapInputError(1)
   endif
 read(cur_fileunit,*,err=1002,end=1002)& 
  (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(j,overlay),j=1, &
     min(flag(1),maxos-1) )  
 endif

 write(READLOG,"(' overlay boundary for overlay: ',I0)" ) overlay
 write(form,"('(1X,',I0,'(f8.3,2X))' )") min(flag(1),maxos-1) 
 write(READLOG,form) (zlevel(i)%cm_zlev(cmx,cmy)%&
  overlay_block(block)%overlay_bon(j, overlay),j=1, min(flag(1),maxos-1) )
	  
! enddo
	  
!handel different block	  
!overlay_num < 0 repeat structure
 if(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay) .lt. 0) then
! do overlay=1, overlay_num_abs
  flag(1)=-1
  write(cur_var,"('rep struc data in Block ',I0, ', Overlay ', I0, ' for cm: ',I0)")  block, overlay, cm_num_1z
  if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  if(flag(1) .lt. 6) then
  write(err_message,"('ERR2003: not enough entries for rep stuct data of block ',I0, ' in cm ',I3)")&
    block,cm_num_1z
  call TrapInputError(1)
endif
read(cur_fileunit,*,err=1002,end=1002) (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(m,overlay),m=1,3),& 
  (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep(m,overlay),m=1,3) 

write(READLOG,"(' Rep. Structure info for overlay: ',I0)") overlay
write(READLOG,"('  repeat times along +x,+y,+z : ', 3(I0,2x))" ) &
 (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(m,overlay),m=1,3)
write(READLOG,"('  dist. between elements along +x,+y,+z : ', 3(f9.3,2x))" ) &
 (zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep(m,overlay),m=1,3)
endif
endif !repeat structure
!different mat number for lattice
if(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .lt. 0) then

rep_ijk(:)=zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(:,overlay)
allocate(zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay)%d3_int(rep_ijk(1), rep_ijk(2),rep_ijk(3) ) )
num_need=rep_ijk(1)*rep_ijk(2)*rep_ijk(3)

allocate(fidobank(num_need))

num_read=0

lineseg=''
do while(num_read .lt. num_need) 
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido
enddo

k=0
do rep_k=1, rep_ijk(3)
do rep_j=1, rep_ijk(2)
do rep_i=1, rep_ijk(1)

k=k+1
if(fidobank(k) .le. 0) then
write(err_message,"('ERR3003: lattice neg./zero mat num, z_lev=',I0, ' cm_ij=', 2(I0,1x), 'block=',I0, 'overlay=',I0)")&
    i, cmx,cmy, block, overlay
call TrapInputError(1)
endif

zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   mat_rep(overlay)%d3_int(rep_i, rep_j, rep_k)=fidobank(k)

end do !rep_i
end do !rep_j
end do !rep_k

write(READLOG,"(' Rep. mat num for overlay: ',I0)") overlay
write(READLOG,"( 20(f4.0,1x) )") fidobank 

deallocate(fidobank)				

endif  !diff. mat number

enddo !overlay number

! endif
flag(1)=0
	 
write(READLOG,"('done overlay block:  ',I2)") block

	  
enddo !enddo block
! endif  !endif block

write(READLOG,*) '------------------------------------'
enddo !enddo cmi
enddo !enddo cmj
 
close(cur_fileunit)
write(READLOG,*)	
write(READLOG,*)  trim(inputfile(i)%fullname), '  : read complete'
write(READLOG, *) '********************************************************'
write(READLOG,*)

enddo  !i zlevel

allocate(cgdim(1)%boundary(0:zlevel(z_start)%ncx),cgdim(2)%boundary(0:zlevel(z_start)%ncy),&
         cgdim(3)%boundary(0:num_zlev) )


cgdim(1)%boundary=zlevel(z_start)%x_cm_pos
cgdim(2)%boundary=zlevel(z_start)%y_cm_pos
cgdim(3)%boundary=z_lev_pos


if(IsZSingle .eq. 0) call ReadIn4
! construct cross section file name prbname.xs
prb_xs%name=prbname
prb_xs%dir=''
prb_xs%ext='xs'
prb_xs%fullname=''


prb_xs%fullname=Getfullname(prb_xs, 0)

1000 return

1001  write(err_message,*) 'open error: file does not exist'
      call TrapOpenError(1)

1002  call TrapReadError(1)

end subroutine

! trap open file errors
subroutine TrapOpenError(stopsign)
integer stopsign
 
  if(stopsign .eq. 0) then
   call TrapError(-1)
  else
   call TrapError(1)
  endif

return
end subroutine


! trap read errors
subroutine TrapReadError(stopsign)
use ErrControl
integer stopsign

if(stopsign .eq. 0) then
 call TrapError(-2)
else
 call TrapError(2)
endif

return
end subroutine

subroutine DisplayMsg
  
 call TrapError(-30)
end subroutine
! trap input errors
subroutine TrapInputError(stopsign)
integer stopsign

if(stopsign .eq. 0) then
 call TrapError(-3)
else
 call TrapError(3)
endif


return
end subroutine

! trap memory allocation error
subroutine TrapMemError(stopsign)
integer stopsign

if(stopsign .eq. 0) then
 call TrapError(-4)
else
 call TrapError(4)
endif

return
end subroutine


! trap errs&warning 
subroutine TrapError(ErrClass)
use ErrControl
use files
integer ErrClass

integer:: i=0 

!Error traps
if (ErrClass .gt. 0) then
 select case (ErrClass)

  !Open file error
  case (1)
    write(*,*) 'open file error'
    write(READLOG,*)  'open file error'
    write(*,*) trim(cur_filename),': ',trim(err_message)
    write(*,*) 'read log file: ',trim(outputfile(1)%fullname)
  
    write(READLOG,*) trim(cur_filename),': ',trim(err_message)
    write(READLOG,*) 'program terminated' 

  !Read file error
  case (2)
    write(*,*) 'Read file error'
    write(READLOG,*)  'Read file error'
    write(*,*) trim(cur_filename),': ',trim(err_message)
    write(*,*) 'read log file: ',trim(outputfile(1)%fullname)
  
    write(READLOG,*) trim(cur_filename),': ',trim(err_message)
    write(READLOG,*) 'program terminated' 


  !Input value error
  case (3) 
    write(*,*) "Input value error in file: ", trim(cur_filename)
    write(READLOG,*)  "Input value error in file: ", trim(cur_filename)
    
	write(*,*)  trim(err_message)
    write(READLOG,*) trim(cur_filename),': ',trim(err_message)
    if(len_trim(err_message_2) .ne. 0) then
     write(READLOG,*) trim(err_message_2)
     write(*,*)  trim(err_message_2)
	 err_message_2=""
	endif

	write(READLOG,*) 'program terminated' 
    write(*,*) 'read log file: ',trim(outputfile(1)%fullname)
 
   !Memory allocation error
   case(4) 
     write(*,*) 'memory allocation error'
	 write(READLOG,*) 'memory allocation error'
    
	 write(*,*) trim(cur_filename),': ',trim(err_message)
     write(*,*) 'read log file: read.log'
  
     write(READLOG,*) trim(cur_filename),': ',trim(err_message)
     write(READLOG,*) 'program terminated' 

  case default
    write(*,"('ErrClass=',I0, ' unexpected')") ErrClass
    write(READLOG,"('ErrClass=',I0, ' unexpected')") ErrClass
  end select
 
 !to die nicely 
  close(READLOG)
!  if (cur_fileunit .ne. 0) close(cur_fileunit)
  stop 'program terminated'

!Warning traps
elseif (ErrClass .lt. 0) then
 
 select case(ErrClass)
  case(-1)
    write(*,*) 'Warning 1004: '
	write(*,*) trim(warning_message)
    write(READLOG,*) 'Warning 1004: ', warning_message
  case(-2)
    write(*,*) 'Warning 1003: ',trim(warning_message)
    write(READLOG,*) 'Warning 1003: ', warning_message
  case(-30)  !false warning 
   write(*,*) trim(warning_message)
   write(READLOG,*) warning_message
  case(-3)
  warning_times(2)=warning_times(2)+1
  if(warning_times(2) .le. max_warning) then
    write(*,*) "Input value warning in file: ", trim(cur_filename)
    write(READLOG,*)  "Input value warning in file: ", trim(cur_filename)
    write(*,*) 'Warning 1002: ',trim(warning_message)
    write(READLOG,*) 'Warning 1002: ', warning_message
	 write(READLOG,*) trim(warning_message_2)
     write(*,*)  trim(warning_message_2)
  endif
  case(-4)
  warning_times(5)=warning_times(5)+1
  if(warning_times(5) .le. max_warning) then
    write(*,*) 'Warning 1005: ',trim(warning_message)
    write(READLOG,*) 'Warning 1005: ', warning_message
  endif
  case (-5)
  warning_times(1)=warning_times(1)+1
  if(warning_times(1) .le. max_warning) then
    write(*,*) 'Warning 1001: ', trim(warning_message)
    write(READLOG,*) 'Warning 1001: ', warning_message
  endif 
  
 case (-6)
  warning_times(3)=warning_times(3)+1
  if(warning_times(3) .le. max_warning) then
    write(*,*) 'Warning 5001: ', trim(warning_message)
    write(READLOG,*) 'Warning 5001: ', warning_message
  endif

  case default
    write(*,"('ErrClass=',I0, ' unexpected')") ErrClass
	write(READLOG,"('ErrClass=',I0, ' unexpected')") ErrClass
	close(READLOG)
    if (cur_fileunit .ne. 0) close(cur_fileunit)
    stop 'program terminated'
 end select

 
!ErrClass=0 not used
else
    write(*,"('ErrClass=',I0, ' unexpected')") ErrClass
	write(READLOG,"('ErrClass=',I0, ' unexpected')") ErrClass
    close(READLOG)
!    if (cur_fileunit .ne. 0) close(cur_fileunit)
    stop 'program terminated'
endif

return
end subroutine

!check the scanned real number whether or not reasonable
subroutine check_real(object,len_object)
use ErrControl
integer len_object
real   object(len_object)

integer i

do i=2, len_object
  if(object(i) .le. object(i-1)) then
    write(err_message,*) trim(cur_var)//' is not in an increasing order'
	call TrapInputError(1)
  endif
enddo

return
end subroutine

!subside warning message

subroutine SubWarning
use ErrControl

integer i
do i=1,10
  if(warning_times(i) .gt. max_warning) then
    write(warning_message,"('Warning ',I4,' happened ',I0, ' times, only ',I0,' displayed, use -w 10 for more' )")&
          i+1000, warning_times(i),max_warning
	call DisplayMsg
  endif
 enddo

end subroutine


!decode lineseg to fidobank
!lineseg is disintegrated into sections
subroutine FeedFido
use mIn4deck
use ErrControl

integer i,k,ierr
integer len_seg
integer :: sec_start=0, sec_end=0, pos_fido=0

character(LEN=50)  secfido
integer :: len_sec=0

integer :: num_has=0
integer :: operee1, operee2, oper_start, oper_end
real    :: operee3, operee4 , step_i
integer pos_i

lineseg=adjustl(lineseg)
len_seg=len_trim(lineseg)

pos_fido=scan(lineseg, '/!')
if(pos_fido .ne. 0) then 
 len_seg=pos_fido-1
endif

!preprocssing the line
do i=1, len_seg
 if(lineseg(i:i) .eq. ',') lineseg(i:i)=' '
enddo

sec_start=1
sec_end=len_seg

do while(sec_start .le. len_seg )

if( num_read .eq. num_need) then
num_read=num_read+1
exit
else if( num_read .gt. num_need) then
exit
endif

!positiong the end of the current sec
sec_end=len_seg
do i=sec_start, len_seg
if(lineseg(i:i) .eq. ' ') then
  sec_end=i
  exit
endif
enddo

len_sec=sec_end-sec_start+1

if(len_sec .gt. 50) then
write(err_message,"('individual fidoset length > 50 in the following line: ')") 
err_message_2=lineseg
call TrapInputError(1)
endif

secfido=lineseg(sec_start : sec_end)

pos_fido=scan(secfido, 'RQZI')


if(pos_fido .eq. 0) then  !No fido character is found
num_has=1
read(secfido, *, iostat=ierr) fidobank(num_read+1)

if(ierr .ne. 0) then
write(err_message,"('error in read non-fido input in the following line: ')") 
 err_message_2=lineseg  
call TrapInputError(1)
endif

num_read=num_read+num_has

else !FIDO character is found

!handle FIDO Character R and Q
select case (secfido(pos_fido:pos_fido))

case('R')
secfido(pos_fido:pos_fido)=' '
read(secfido, *, iostat=ierr) operee1, operee3
secfido(pos_fido:pos_fido)='R'

if(ierr .ne. 0) then
write(err_message,"('error in read fido input in the following line: ')") 
err_message_2=trim(lineseg)  
call TrapInputError(1)
endif

num_has=operee1
oper_start=num_read+1
oper_end=min(num_read+num_has, num_need)
fidobank(oper_start:oper_end)=operee3
num_read=num_read+num_has

case('Q')
secfido(pos_fido:pos_fido)=' '
read(secfido, *, iostat=ierr) operee1, operee2
secfido(pos_fido:pos_fido)='Q'

if(ierr .ne. 0) then
write(err_message,"('error in read fido input in the following line: ')")  
err_message_2=lineseg 
call TrapInputError(1)
endif

num_has=operee1*operee2

if(num_read .lt. operee2) then 
write(err_message,"('FIDO Q error in the following line: (not enough previous numbers to repeat)')")  
err_message_2=trim(lineseg)
call TrapInputError(1)
endif

do i=1, operee1
oper_start=num_read+(i-1)*operee2+1
oper_end=min(oper_start+operee2-1, num_need)
fidobank(oper_start: oper_end)=fidobank(num_read-operee2+1:num_read)
enddo
num_read=num_read+num_has

case('Z')
  secfido(pos_fido:pos_fido)=' '
  read(secfido, *, iostat=ierr) operee1
  if(ierr .ne. 0) then
    secfido(pos_fido:pos_fido)='Z'
    write(err_message,"('error in read fido input in the following line:  ')")  
    err_message_2=lineseg 
    call TrapInputError(1)
   endif
   num_has=operee1
   if (num_has+num_read .le. num_need) then
     fidobank(num_read+1:num_read+num_has)=0
   else 
     fidobank(num_read+1:num_need)=0
   endif
   num_read=num_read+num_has
case('I')
   secfido(pos_fido:pos_fido)=' '
   pos_i=scan(secfido, ':')
   if(pos_i .eq. 0) then 
     secfido(pos_fido:pos_fido)='I'
     write(err_message,"('wrong format for FIDO character I: ', A)")  trim(secfido)
     err_message_2=lineseg 
     call TrapInputError(1)
   endif
   secfido(pos_i:pos_i)=' '
   read(secfido, *, iostat=ierr) operee1, operee3, operee4

   if(ierr .ne. 0) then
     secfido(pos_fido:pos_fido)='I'
     write(err_message,"('error in read fido I:  ', A)")  trim(secfido)
     err_message_2=lineseg 
     call TrapInputError(1)
   endif
   num_has=operee1+2
   step_i=(operee4-operee3)/(operee1+1)
   k=1
   do while (k .le. min(num_has, num_need-num_read) )
     fidobank(num_read+k)=operee3+(k-1)*step_i
     k=k+1
   enddo
 
   if (num_has+num_read .le. num_need ) fidobank(num_read+num_has)=operee4
   num_read=num_read+num_has

end select

endif !FIDO/no FIDO

!for next section   
sec_start=0
do i=sec_end+1, len_seg
  if(lineseg(i:i) .ne. ' ') then
    sec_start=i
    exit
  endif
enddo
if (sec_start .eq. 0) sec_start=len_seg+1

enddo !end do while

!if( num_need .gt. num_read) then
!write(err_message,"('ERR3004: not enough entries (FIDO) , Found=',I0, ' Required=', I0)")&
!    num_read, num_need
!err_message_2=trim(lineseg)
!call TrapInputError(1)
!endif

if (num_need .lt. num_read) then

write(warning_message,"('too many entres (FIDO) , Found>=',I0, ' Required=', I0)")&
    num_read, num_need
warning_message_2=trim(lineseg)
call TrapInputError(0)

endif


end subroutine

!check the entries in lineseg
!return 
! lineflag=1  : lineseg is a string
! lineflag=0  : lineseg is a series of fido entries or numbers
subroutine CheckLineseg
use mIn4deck

integer i
integer len_seg
integer ::  pos_space=0, pos_digit=0, pos_fidoc=0

lineseg=adjustl(lineseg)
len_seg=len_trim(lineseg)
! firstsec=''

lineflag=0
do i=1, len_seg
if(lineseg(i:i) .eq. ' ') then
 pos_space=i
 exit
endif
enddo

if(pos_space .eq. 0) pos_space=len_seg

do i=1, pos_space

 pos_digit=index(digchar, lineseg(i:i))
 
 if (pos_digit .eq. 0 ) then
   if( pos_fidoc .ne. 0) then
     lineflag=1
     exit
   endif
   pos_fidoc=index(fidochar, lineseg(i:i))
   if (pos_fidoc .eq. 0) then
    lineflag=1
    exit
   endif
 else
   pos_fidoc=0
 endif

enddo
end subroutine

!read all-in-one input file (combined all z levels) prb.inp
subroutine ReadPrbInp
use files
use funs
use ErrControl
use lineread
use paraset1
use paraset2
use paraset4
use mIn4deck

integer i
integer :: ierr=0, len_eff=0

integer pos_mark
integer cmx,cmy, cmz, cm_idx, matn
integer :: cijk(3)=0, fijk(3)=0

integer :: block, num_block, temp_int

integer overlay,overlay_typ,overlay_num_abs, over_num
integer :: rep_ijk(3), rep_i, rep_j, rep_k

character(len=50) :: cm_var=''
integer, allocatable :: cm_read_flag(:,:,:)

 write(READLOG,*) 'opening inputfile: ', trim(cur_filename)

 cur_fileunit=88
 open(cur_fileunit, file=cur_filename,STATUS='OLD',ERR=1001)
 
  
 cur_var='ncx,ncy... '
 
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
    read(cur_fileunit,*,err=1002,end=1002) cmx, cmy  
 endif 

if(cmx .le. 0 .or.  cmy .le. 0 ) then
  write(err_message,*) trim(cur_var)//'less than zero:  input value error or read failure (ERR1004)'
  call TrapInputError(1)
endif

do i=1 , num_zlev
 zlevel(i)%ncx=cmx
 zlevel(i)%ncy=cmy
 zlevel(i)%num_zfine=0
 allocate(zlevel(i)%cm_zlev(zlevel(i)%ncx,zlevel(i)%ncy)) 
 allocate(zlevel(i)%x_cm_pos(1+zlevel(i)%ncx))
 allocate(zlevel(i)%y_cm_pos(1+zlevel(i)%ncy))
 
enddo

 write(READLOG,*) 'number of coarse mesh along x(ncx):   ', cmx
 write(READLOG,*) 'number of coarse mesh along x(ncy):   ', cmy
 write(READLOG,*)

 num_cmesh(1)=cmx
 num_cmesh(2)=cmy
 num_cmesh(3)=num_zlev
 
 allocate(cm_read_flag(cmx,cmy,num_zlev))
 cm_read_flag=0
 
 cur_var='x_cm_pos '
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  read(cur_fileunit,*,err=1002,end=1002) zlevel(1)%x_cm_pos
 endif
 
 call check_real(zlevel(1)%x_cm_pos,1+zlevel(1)%ncx)

 write(READLOG,*) 'coarse mesh position along x: '
 write(READLOG,*) zlevel(1)%x_cm_pos
 write(READLOG,*)

 cur_var='y_cm_pos'
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  read(cur_fileunit,*,err=1002,end=1002) zlevel(1)%y_cm_pos
 endif

call check_real(zlevel(1)%y_cm_pos,1+zlevel(1)%ncy)

 write(READLOG,*) 'coarse mesh position along y: '
 write(READLOG,*) zlevel(1)%y_cm_pos
 write(READLOG,*)

do i=2 , num_zlev
 zlevel(i)%x_cm_pos=zlevel(1)%x_cm_pos
 zlevel(i)%y_cm_pos=zlevel(1)%y_cm_pos
enddo
 cur_var='corase mesh cards'

ierr=0
loop_line : do while (ierr .eq. 0)

in4line=''
read(cur_fileunit,"(A)",iostat=ierr) in4line
if(ierr .ne. 0 ) exit loop_line

in4line=adjustl(trim(in4line))
len_eff=len_trim(in4line)

!skip blank line
if(len_eff .eq. 0) cycle loop_line 
!check comment line
if( index(cmtchar, in4line(1:1)) .ne. 0) cycle loop_line

!check coarsh mesh marker
pos_mark=index(in4line(1:len_eff), cmmarker)
if( pos_mark .eq. 0) then
  write(err_message,"('unexpected line: ', A )")  in4line(1:len_eff)
  call TrapInputError(1)
endif

!read cm index (i,j,k) after the marker
read(in4line(pos_mark+3:len_eff), *, iostat=ierr) cijk
if(ierr .ne. 0) then
  write(err_message,"('coarse mesh reading index error: ', A )")  in4line(1:len_eff)
  call TrapInputError(1)
endif

do i=1, 3
if(cijk(i) .le. 0 .or. cijk(i) .gt. num_cmesh(i) ) then
  write(err_message,"('coarse mesh index error: ', A )")  in4line(1:len_eff)
  call TrapInputError(1)
endif
enddo

cmx=cijk(1)
cmy=cijk(2)
cmz=cijk(3)
cm_idx=num_cmesh(1)*num_cmesh(2)*(cmz-1)+num_cmesh(1)*(cmy-1) + cmx

write(cm_var, "('cm(', I0, ',' , I0 , ',' , I0, ')' )" ) cijk

write(cur_var, "(A, ' mat number')" ) trim(cm_var)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit,*,err=1002,end=1002)  matn
endif
   
zlevel(cmz)%cm_zlev(cmx,cmy)%cm_mat_num=abs(matn)

if(zlevel(cmz)%cm_zlev(cmx,cmy)%cm_mat_num .gt. num_material .or. zlevel(cmz)%cm_zlev(cmx,cmy)%cm_mat_num .le. 0) then
   write(err_message, "('Material # > num of materials at : ',A )") trim(cm_var)
   call TrapInputError(1)
endif 

write(READLOG, "(A, 2x, '  :  ', I0)") trim(cur_var), zlevel(cmz)%cm_zlev(cmx,cmy)%cm_mat_num

write(cur_var, "('cm(', I0, ',' , I0 , ',' , I0, ') fine mesh number along x, y, and z' )" ) cijk
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
 read(cur_fileunit,*,err=1002,end=1002) zlevel(cmz)%cm_zlev(cmx,cmy)%fine
endif

! re-size mesh
if (max_fm_size.gt.0) then

   zlevel(cmz)%cm_zlev(cmx,cmy)%fine(1)=ceiling((zlevel(cmz)%x_cm_pos(cmx+1)-zlevel(cmz)%x_cm_pos(cmx))/max_fm_size)*zlevel(cmz)%cm_zlev(cmx,cmy)%fine(1)
   zlevel(cmz)%cm_zlev(cmx,cmy)%fine(2)=ceiling((zlevel(cmz)%y_cm_pos(cmy+1)-zlevel(cmz)%y_cm_pos(cmy))/max_fm_size)*zlevel(cmz)%cm_zlev(cmx,cmy)%fine(2)
   zlevel(cmz)%cm_zlev(cmx,cmy)%fine(3)=ceiling((z_lev_pos(cmz+1)-z_lev_pos(cmz))/max_fm_size)*zlevel(cmz)%cm_zlev(cmx,cmy)%fine(3)

endif

do i=1, 3
if(zlevel(cmz)%cm_zlev(cmx,cmy)%fine(i) .le. 0) then
   write(err_message, "('fine mesh number less than zero at : ',A )") trim(cm_var)
   call TrapInputError(1)
endif 
enddo

write(READLOG, "(A, 2x, '  :  ', 3(I0,1x) )") trim(cur_var), zlevel(cmz)%cm_zlev(cmx,cmy)%fine

zlevel(cmz)%cm_zlev(cmx,cmy)%num_subregion=1
zlevel(cmz)%cm_zlev(cmx,cmy)%cm_type=1

!read overlay
if(matn .gt. 0) then
cm_read_flag(cmx,cmy,cmz)=cm_read_flag(cmx,cmy,cmz)+1
cycle loop_line
endif

i=cmz
zlevel(cmz)%cm_zlev(cmx,cmy)%cm_type=-1

 num_block=1
 allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(num_block))
 write(READLOG,*)
 write(READLOG,"('Overlay information for coarse mesh: ',A)") trim(cur_var)

 do block=1,num_block
!handle repeat structure
  write(READLOG,"('overlay block:  ',I2)") block
  zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num=0
  write(READLOG,*)
  write(cur_var,"('num of overlay for cm:  ',A)")  trim(cm_var)

  if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
   read(cur_fileunit,*,err=1002,end=1002) zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num
  endif
  over_num=zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num
  overlay_num_abs=abs(over_num)
  
  write(READLOG,"(' num of overlay   : ', I0)")  &
   zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_num

  allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    overlay_typ(overlay_num_abs))
  allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
       overlay_mat(overlay_num_abs))

  allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    num_rep(3,overlay_num_abs), &
    zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
    dist_rep(3,overlay_num_abs) )
     
  zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep=1
  zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep=0.0

  write(cur_var,"('type of overlay for cm:  ',A)")  trim(cm_var)

num_need=overlay_num_abs
allocate(fidobank(num_need))

num_read=0

lineseg=''

do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido

enddo

zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   overlay_typ(:)=fidobank(:)

deallocate(fidobank)

  if(over_num .lt. 0) then
  do j=1, overlay_num_abs
  temp_int=zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j)
  zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j)=-abs(temp_int)
  enddo
  endif
  
  write(form,"('(A20,',I0,'(I0,2x))' )") overlay_num_abs   

  write(READLOG,form) ' overlay type     : ', &
   zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ
  
  write(cur_var,"('mat num of overlay for cm:  ',A)")  trim(cm_var)
num_need=overlay_num_abs
allocate(fidobank(num_need))

num_read=0

lineseg=''
do while (num_read .lt. num_need)
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido
enddo
zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   overlay_mat(:)=fidobank(:)

deallocate(fidobank)				

do j=1, overlay_num_abs
  temp_int=zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(j)
  if(abs(temp_int) .gt. num_material .or. temp_int .eq. 0 ) then
  write(err_message, "('Overlay Material # > num of materials at ', A , 'overlay=', I0)") & 
     cm_var, j
   call TrapInputError(1)
  endif
  if(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(j) .gt. 0) &
   zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(j)=abs(temp_int)
enddo

 
write(form,"('(A20,',I0,'(I0,2x))' )") overlay_num_abs 	  
write(READLOG,form) ' overlay mat      : ', &
   zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat

allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos, &
   overlay_num_abs))

allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay_num_abs))

zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos,&
      overlay_num_abs)=0
	  
do overlay=1,overlay_num_abs
 if (zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .gt. num_material) then  
  write(err_message, "('ERR2002: Overlay Mat # > num of materials at ',A, 'overlay=', I0)") & 
     trim(cm_var), overlay
  call TrapInputError(1)
endif

if(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .lt. 0 .and. &
 zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay) .gt. 0) then
write(err_message, "('ERR2003: Overlay Mat # < 0 at ', A , ' overlay ',I0)") &
  trim(cm_var), overlay
 call TrapInputError(1)
endif
 overlay_typ=zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay)
 flag(1)=-1
 write(cur_var,"('overlay boundaries of overlay ',I0,' for : ',A)")  overlay,trim(cm_var)

 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
   zlevel(i)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(maxos,overlay)=flag(1)
   if(flag(1) .lt. map_overlay_typ2bon(abs(overlay_typ),cm_idx)) then
    write(err_message,"('ERR2003: not enough entries for boundaries of overlay ',I0, ' in cm ',A)")&
        overlay,trim(cm_var)
    call TrapInputError(1)
   endif
 read(cur_fileunit,*,err=1002,end=1002)& 
  (zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_bon(j,overlay),j=1, &
     min(flag(1),maxos-1) )  
 endif

 write(READLOG,"(' overlay boundary for overlay: ',I0)" ) overlay
 write(form,"('(1X,',I0,'(f8.3,2X))' )") min(flag(1),maxos-1) 
 write(READLOG,form) (zlevel(cmz)%cm_zlev(cmx,cmy)%&
  overlay_block(block)%overlay_bon(j, overlay),j=1, min(flag(1),maxos-1) )

!overlay_typ < 0 repeat structure
 if(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_typ(overlay) .lt. 0) then
! do overlay=1, overlay_num_abs
  flag(1)=-1
  write(cur_var,"('rep struc data in Block ',I0, ', Overlay ', I0, ' for: ', A)")  block, overlay, trim(cm_var)
  if(NumCmtLine(cur_fileunit,flag) .ge. 0) then
  if(flag(1) .lt. 6) then
  write(err_message,"('ERR2003: not enough entries for rep stuct data of block ',I0, ' in  ',A)")&
    block,trim(cm_var)
  call TrapInputError(1)
endif
read(cur_fileunit,*,err=1002,end=1002) (zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(m,overlay),m=1,3),& 
  (zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep(m,overlay),m=1,3) 

write(READLOG,"(' Rep. Structure info for overlay: ',I0)") overlay
write(READLOG,"('  repeat times along +x,+y,+z : ', 3(I0,2x))" ) &
 (zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(m,overlay),m=1,3)
write(READLOG,"('  dist. between elements along +x,+y,+z : ', 3(f9.3,2x))" ) &
 (zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%dist_rep(m,overlay),m=1,3)
endif
endif !repeat structure
!different mat number for lattice
if(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%overlay_mat(overlay) .lt. 0) then

rep_ijk(:)=zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%num_rep(:,overlay)
allocate(zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)%mat_rep(overlay)%d3_int(rep_ijk(1), rep_ijk(2),rep_ijk(3) ) )
num_need=rep_ijk(1)*rep_ijk(2)*rep_ijk(3)

allocate(fidobank(num_need))

num_read=0

lineseg=''
do while(num_read .lt. num_need) 
if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

call FeedFido
enddo

k=0
do rep_k=1, rep_ijk(3)
do rep_j=1, rep_ijk(2)
do rep_i=1, rep_ijk(1)

k=k+1
if(fidobank(k) .le. 0) then
write(err_message,"('ERR3003: lattice neg./zero mat num, z_lev=',I0, ' cm_ij=', 2(I0,1x), 'block=',I0, 'overlay=',I0)")&
    cmz, cmx,cmy, block, overlay
call TrapInputError(1)
endif

zlevel(cmz)%cm_zlev(cmx,cmy)%overlay_block(block)% &
   mat_rep(overlay)%d3_int(rep_i, rep_j, rep_k)=fidobank(k)

end do !rep_i
end do !rep_j
end do !rep_k

write(READLOG,"(' Rep. mat num for overlay: ',I0)") overlay
write(READLOG,"( 20(f4.0,1x) )") fidobank 

deallocate(fidobank)				

endif  !diff. mat number

enddo !overlay number

! endif
flag(1)=0
 
write(READLOG,"('done overlay block:  ',I2)") block
 
enddo !enddo block

cm_read_flag(cmx,cmy,cmz)=cm_read_flag(cmx,cmy,cmz)+1

enddo loop_line

close(cur_fileunit)

do cmz=1, num_cmesh(3)
do cmy=1, num_cmesh(2)
do cmx=1, num_cmesh(1)

if(cm_read_flag(cmx,cmy,cmz) .eq. 0) then
write(err_message,"('cm is not defined:  cm(',I0, ',',I0,',',I0, ')' )") cmx, cmy, cmz
call TrapInputError(1)
endif

if(cm_read_flag(cmx,cmy,cmz) .gt. 1) then
write(err_message,"('cm is not defined more than once:  cm(',I0, ',',I0,',',I0, ')' )") cmx, cmy, cmz
call TrapInputError(1)
endif


enddo
enddo
enddo

write(READLOG,*) '------------------------------------'

 
write(READLOG,*)	
write(READLOG,*)  trim(cur_filename), '  : read complete'
write(READLOG, *) '********************************************************'
write(READLOG,*)


return

1001  write(err_message,*) 'open error'
      call TrapOpenError(1)

1002  call TrapReadError(1)
end subroutine
