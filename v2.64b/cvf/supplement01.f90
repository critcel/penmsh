!supplement01: for heart model
! heart prb in supplement01
Module heart
use files

!@dist(1,2,3): x,y,z size of the model
!@num_pix(1,2,3): x, y,z pixel size 
 real dist(3)
 integer num_pix(3)
!dist_pix : distance per pixel
 real  dist_pix(3)
 integer num_point_hrt, num_binary
 

type tHrtfile
   
 integer num_data
 integer:: imax
!data_type=0  : 4 byte integer
!data_type=1  : 4 byte real
 integer :: data_type=0
 integer, dimension(:,:,:), allocatable  :: image
 real, dimension(:,:,:), allocatable :: rmage 
 integer*2,dimension(:,:,:), allocatable :: matrix
 real, dimension(:), allocatable::  list
 
end type

type (tHrtfile), pointer :: hrt_data(:)
type(fileinfo), pointer ::  hrt_file(:) 

type(fileinfo),save :: hrt_deck
 
end module


!read heart model binary activity and attenuation input files

subroutine ReadHrtInput
use heart
use files
use ErrControl
use paraset1
use paraset2
use paraset3
use paraset4
use lineread

integer i,j,k,m
integer dir_len,name_len

character*100 nameline
integer len_nameline


hrt_deck%dir=input_dir

if (hrt_filename .ne. '') then
 hrt_deck%fullname=hrt_filename
else
 hrt_deck%fullname='heart.inp'
endif

write(*,"('Opening input deck:', A)")   hrt_deck%fullname
write(READLOG,"('Opening input deck:', A)")   hrt_deck%fullname


dir_len=len_trim(hrt_deck%dir)
name_len=len_trim(hrt_deck%fullname)
write(form,"( '(','A',I0,',A1,','A',I0,')' )" ) dir_len, name_len
write(cur_filename,fmt=form) hrt_deck%dir ,'/',hrt_deck%fullname
cur_fileunit=1
open(UNIT=1, file=trim(cur_filename),STATUS='OLD',err=1001)
write(READLOG,*) 'heart input deck open success'

cur_var='problem name'
if(NumCmtLine(1,flag) .ge. 0) then
 read(1,'(A)',err=1002,end=1002) prbname
endif

write(READLOG,*) 'problem name : ', trim(prbname)
write(READLOG,*)

prbname_len=len_trim(prbname)

cur_var='number of pixels along x, y, and z'
num_pix=0
if(NumCmtLine(1,flag) .ge. 0) then
 read(1,*,err=1002,end=1002) num_pix
endif

write(READLOG,*) 'num of pixels along x:  ', num_pix(1)
write(READLOG,*) 'num of pixels along y:  ', num_pix(2)
write(READLOG,*) 'num of pixels along z:  ', num_pix(3)
write(READLOG,*) 

m=0
do k=1, 3
if(num_pix(k) .le. 0) then
  m=-1
  exit
endif
enddo
if(m .eq. -1) then
 write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1501)'
 call TrapInputError(1)
endif

cur_var='pixel size along x, y, and z'
dist_pix=0
if(NumCmtLine(1,flag) .ge. 0) then
 read(1,*,err=1002,end=1002) dist_pix
endif

write(READLOG,*) 'pixel size along x:  ', dist_pix(1)
write(READLOG,*) 'pixel size along y:  ', dist_pix(2)
write(READLOG,*) 'pixel size along z:  ', dist_pix(3)
write(READLOG,*) 

m=0
do k=1, 3
if(dist_pix(k) .le. 0) then
  m=-1
  exit
endif
enddo
if(m .eq. -1) then
 write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1502)'
 call TrapInputError(1)
endif

num_point_hrt=0

cur_var='binary file name(s)'
if(NumCmtLine(1,flag) .ge. 0) then
 read(1,"(A)",err=1002,end=1002) nameline
endif

nameline=ADJUSTL(nameline)
len_nameline=len_trim(nameline)
num_binary=index(nameline(1:len_nameline),' ')
if (num_binary .eq. 0 .or. num_binary .ge. len_nameline) then
num_binary=1
else
num_binary=2
endif

allocate(hrt_file(num_binary), hrt_data(num_binary) )

cur_var='binary data file type'
flag(1)=-1

if(NumCmtLine(1,flag) .ge. 0) then
 if(flag(1) .lt. 2 .or. num_binary .eq. 1) then
  read(1,*,err=1002,end=1002) hrt_data(1)%data_type
 else if (flag(1) .ge. 2 .and. num_binary .eq. 2) then
  read(1,*,err=1002,end=1002) hrt_data(1)%data_type,hrt_data(2)%data_type
 endif
endif

flag(1)=0

if(num_binary .eq. 1) then
 read(nameline, *) hrt_file(1)%fullname
else
 read(nameline, *) hrt_file(1)%fullname, hrt_file(2)%fullname
endif

write(READLOG,*) 'mat binary filename:  ', hrt_file(1)%fullname
if (num_binary .gt. 1) &
write(READLOG,*) 'src binary filename:  ', hrt_file(2)%fullname
write(READLOG,*) 

cur_var='number of coarse mesh along x, y, and z'
num_cmesh=0
if(NumCmtLine(1,flag) .ge. 0) then
 read(1,*,err=1002,end=1002) num_cmesh
endif

write(READLOG,*) 'num of coarse mesh along x:  ', num_cmesh(1)
write(READLOG,*) 'num of coarse mesh along y:  ', num_cmesh(2)
write(READLOG,*) 'num of coarse mesh along z:  ', num_cmesh(3)
write(READLOG,*) 

m=0
do k=1, 3
if(num_cmesh(k) .le. 0) then
  m=-1
  exit
endif
enddo
if(m .eq. -1) then
 write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1503)'
 call TrapInputError(1)
endif

!setup meshing
num_cm=1
do i=1,3
 num_cm=num_cm*num_cmesh(i)
enddo
num_cm_zlev=num_cmesh(1)*num_cmesh(2)
allocate(cgdim(1)%pixb(num_cmesh(1)),cgdim(2)%pixb(num_cmesh(2)),cgdim(3)%pixb(num_cmesh(3)))
allocate(cgdim(1)%boundary(0:num_cmesh(1)),cgdim(2)%boundary(0:num_cmesh(2)),&
         cgdim(3)%boundary(0:num_cmesh(3)) )

do i=1, 3
 cgdim(i)%pixb=0
enddo

cur_var='number of pixels per coarse mesh along x, y, and z'
do i=1,3
 cgdim(i)%pixb=0
enddo

if(NumCmtLine(1,flag) .ge. 0) then
 do i=1,3
 read(1,*,err=1002,end=1002) cgdim(i)%pixb
 enddo
endif

write(READLOG,*) 'num of pixel per coarse mesh along x:  ', cgdim(1)%pixb
write(READLOG,*) 'num of pixel per coarse mesh along y:  ', cgdim(2)%pixb
write(READLOG,*) 'num of pixel per coarse mesh along z:  ', cgdim(3)%pixb
write(READLOG,*) 

m=0
do k=1, 3
do i=1, num_cmesh(k)
if(cgdim(k)%pixb(i) .le. 0) then
  m=-1
  exit
endif
enddo
enddo

if(m .eq. -1) then
 write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1504)'
 call TrapInputError(1)
endif

cur_var='source_format, xnum_s_mesh ... '

num_group=-1
sn_order=-1
pn_order=-1
if(NumCmtLine(1,flag) .ge. 0) then
 read(1,*,err=1002,end=1002) s_format, xnum_s_mesh, ynum_s_mesh, znum_s_mesh,num_group, sn_order,pn_order
endif

write(READLOG,*) 'source format: ', s_format
write(READLOG,*) 'x mesh in src grid(.src): ', xnum_s_mesh
write(READLOG,*) 'y mesh in src grid(.src): ', ynum_s_mesh
write(READLOG,*) 'z mesh in src grid(.src): ', znum_s_mesh
write(READLOG,*) 'num of group: '         , num_group
write(READLOG,*) 'Sn order: ', sn_order
write(READLOG,*) 'Pn order: ', pn_order

if(num_group .lt. 0 .or. sn_order .lt. 0 .or. pn_order .lt. 0) then
  write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1505)'
 call TrapInputError(1)
endif

write(READLOG,*)


cur_var='xs_format,num_comment...'
num_comment=-1
leg_order=-1
len_table=-1
if(NumCmtLine(1,flag) .ge. 0)  then
 read(1,*,err=1002,end=1002) xs_format,num_comment,leg_order,len_table
endif 

write(READLOG, *) 'cross section information:'
write(READLOG, *) 'xs format:     ', xs_format
write(READLOG, *) 'num_comment:   ', num_comment
write(READLOG, *) 'legdre order:  ', leg_order
write(READLOG, *) 'table length:  ', len_table
write(READLOG, *)

if(num_comment .lt. 0 .or. leg_order .lt. 0 .or. len_table .lt. 0) then
  write(err_message,*) trim(cur_var)//' : input value error or read failure (ERR1506)'
  call TrapInputError(1)
endif

cur_var='boundary cond. '
if(NumCmtLine(1,flag) .ge. 0) then 
! ibback, ibfrnt, jbright,jbleft,kbottom,ktop
 read(1,*,err=1002,end=1002) bcs
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
 if(NumCmtLine(1,flag) .ge. 0)  read(1,*,err=1002,end=1002) id_s_mat
 write(READLOG,*) 'source material: ', id_s_mat
 
 cur_var='source intensity'
 if(NumCmtLine(1,flag) .ge. 0)  read(1,*,err=1002,end=1002) s_intensity
 write(READLOG,*) 'source intensity(neg value as total magnitude): ', s_intensity

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


write(*,"('phantom-like input deck done:', A)")   hrt_deck%fullname
write(READLOG,"('phantom-like input deck done:', A)")   hrt_deck%fullname

call ReadIn4
call InitMeshing


!handle binary files
dir_len=len_trim(hrt_deck%dir)
name_len=len_trim(hrt_file(1)%fullname)
write(form,"( '(','A',I0,',A1,','A',I0,')' )" ) dir_len, name_len
write(cur_filename,fmt=form) hrt_deck%dir ,'/',hrt_file(1)%fullname
cur_fileunit=61
write(*,"('Opening binary file :', A)")   cur_filename
write(READLOG,"('Opening binary:', A)")   cur_filename

open(UNIT=61, file=trim(cur_filename),FORM='BINARY',STATUS='OLD',err=1001)
write(READLOG,"('heart-phantom-like mat binary file open success: ',A)") trim(cur_filename) 

if (hrt_data(1)%data_type .le. 0 ) then
allocate (hrt_data(1)%image(num_pix(1),num_pix(2),num_pix(3)))
read(61,ERR=1002) hrt_data(1)%image
else
allocate (hrt_data(1)%rmage(num_pix(1),num_pix(2),num_pix(3)))
read(61,ERR=1002) hrt_data(1)%rmage
endif


close(61)

 
call ProcessData(1)

if (hrt_data(1)%data_type .le. 0) then
deallocate(hrt_data(1)%image)
else
deallocate(hrt_data(1)%rmage)
endif

call DefineMatMatrix


if (num_binary .eq. 2) then
dir_len=len_trim(hrt_deck%dir)
name_len=len_trim(hrt_file(2)%fullname)
write(form,"( '(','A',I0,',A1,','A',I0,')' )" ) dir_len, name_len
write(cur_filename,fmt=form) hrt_deck%dir ,'/',hrt_file(2)%fullname
cur_fileunit=62
write(*,"('Opening binary file :', A)")   cur_filename
write(READLOG,"('Opening binary:', A)")   cur_filename
open(UNIT=62, file=trim(cur_filename),FORM='BINARY',STATUS='OLD',err=1001)
write(READLOG,"('heart-phantom-like src binary file open success: ',A)") trim(cur_filename) 

!allocate (hrt_data(2)%image(num_pix(1),num_pix(2),num_pix(3)))

if (hrt_data(2)%data_type .eq. 0) then
allocate (hrt_data(2)%image(num_pix(1),num_pix(2),num_pix(3)))
read(62,ERR=1002) hrt_data(2)%image
else
allocate (hrt_data(2)%rmage(num_pix(1),num_pix(2),num_pix(3)))
read(62,ERR=1002) hrt_data(2)%rmage
endif

close(62)


call ProcessData(2)

if (hrt_data(2)%data_type .eq. 0) then
deallocate(hrt_data(2)%image)
else
deallocate(hrt_data(2)%rmage)
endif


endif

call DefineSrcMatrix

close (1)
call AssignPenmshVaribles

open (unit=91, file=trim(prbname)//"_data.out")
do i=1, num_binary
write(91, *)
write(91, '(" different values read in file", A," : ", I0 )') &
    trim(hrt_file(i)%fullname), hrt_data(i)%num_data
write(91, *)
if(i .eq. 1) then
write(91, "('    #          value in phantom      num assigned')") 
else
write(91, "('    #          value in phantom')") 
endif
do j=1, hrt_data(i)%num_data 
 if( i .eq. 1) then
  if( hrt_data(1)%data_type .ge. 0) then 
   write(91,"(I6,8x,ES12.5,11x,I0 )")  j, hrt_data(i)%list(j), j
  else
   write(91,"(I6,8x,ES12.5, 11x, I0)")  j, hrt_data(i)%list(j),int(hrt_data(i)%list(j))
  endif
 else
   write(91,"(I6,8x,ES12.5)")  j, hrt_data(i)%list(j)
 endif
enddo
write(91, *)

write(91,"('   total  ', I4)")  hrt_data(i)%num_data

enddo
close (91)
write(*,*)
write(READLOG,*) 
write(*,"('binary file(s) done ')") 
write(READLOG,"('binary file(s) done ')") 

write(READLOG,"('data list in file : ',A )") trim(prbname)//"_data.out"
write(*,"('data list in file : ',A )") trim(prbname)//"_data.out"

write(*,*)
write(READLOG,*) 


1000 return

1001  write(err_message,*) 'open error: file does not exist'
      call TrapOpenError(1)

1002  call TrapReadError(1)

end subroutine

! processing the attenuation and activity file
subroutine ProcessData(filenum)
use heart

integer filenum

integer i,j,k,m,numdata
integer daflag,ivalue
real,dimension(:), allocatable:: dl_temp
real value 

allocate (hrt_data(filenum)%matrix(num_pix(1),num_pix(2),num_pix(3)))
numdata=1
allocate (hrt_data(filenum)%list(numdata))
! count how many different numbers in the attenuation file(material number)
if (hrt_data(filenum)%data_type .le. 0) then
value= hrt_data(filenum)%image(1,1,1)
else
value=hrt_data(filenum)%rmage(1,1,1)
endif

hrt_data(filenum)%list(1)= value
hrt_data(filenum)%imax=0

do k=1,num_pix(3)
  do j=1,num_pix(2)
    do i=1,num_pix(1)
       daflag=0
       if (hrt_data(filenum)%data_type .le. 0 ) then
         value= hrt_data(filenum)%image(i,j,k)
       else
         value= hrt_data(filenum)%rmage(i,j,k)
    endif
 
    do m=1,numdata
      if(value .eq. hrt_data(filenum)%list(m)) daflag=m
    enddo
    if (daflag .eq. 0) then
      allocate (dl_temp(numdata))
   dl_temp=hrt_data(filenum)%list
   deallocate(hrt_data(filenum)%list)
   allocate(hrt_data(filenum)%list(numdata+1))
         do m=1,numdata
     hrt_data(filenum)%list(m)=dl_temp(m)
   enddo
   hrt_data(filenum)%list(numdata+1)=value 
   deallocate(dl_temp)
   numdata=numdata+1
!      hrt_data(filenum)%matrix(i,j,k)=numdata
      if(hrt_data(filenum)%data_type .eq. -1) then
        ivalue=hrt_data(filenum)%image(i,j,k)
               hrt_data(filenum)%matrix(i,j,k)=ivalue
               if(hrt_data(filenum)%imax .lt. ivalue) hrt_data(filenum)%imax=ivalue
      else
              hrt_data(filenum)%matrix(i,j,k)=numdata
      endif
    else 
!      hrt_data(filenum)%matrix(i,j,k)=daflag  
      if(hrt_data(filenum)%data_type .eq. -1 ) then
     hrt_data(filenum)%matrix(i,j,k)=hrt_data(filenum)%image(i,j,k)
   else
        hrt_data(filenum)%matrix(i,j,k)=daflag
   endif

    endif

    enddo
  enddo
enddo    
hrt_data(filenum)%num_data=numdata

return

end subroutine


subroutine DefineMatMatrix
use heart
use paraset1
use paraset2
use paraset3
use paraset4
use files
use ErrControl

integer i,j,k,m, matone
integer x0,y0,z0,x,y,z,ini,inj,ink

if(hrt_data(1)%data_type .eq. -1) then
num_material=hrt_data(1)%imax
else
num_material=hrt_data(1)%num_data
endif

allocate(mat_vol(num_material),mat_fm(num_material))
mat_vol=0
mat_fm=0
total_vol=(z_lev_pos(num_zlev+1)-z_lev_pos(1))*&
          (zlevel(1)%y_cm_pos(1+zlevel(1)%ncy)-zlevel(1)%y_cm_pos(1))*&
    (zlevel(1)%x_cm_pos(1+zlevel(1)%ncx)-zlevel(1)%x_cm_pos(1))

do k=1, num_cmesh(3)
 do j=1,num_cmesh(2)
  do i=1, num_cmesh(1) 
  
   x0=0
   y0=0
   z0=0
  
  do m=1,k-1
    z0=z0+cgdim(3)%pixb(m)
  enddo
  z0=z0+1
 
  do m=1,j-1
    y0=y0+cgdim(2)%pixb(m)
  enddo
  y0=y0+1
 
  do m=1,i-1
    x0=x0+cgdim(1)%pixb(m)
  enddo
  x0=x0+1
 
 allocate(zlevel(k)%cm_zlev(i,j)%fm_num_mat(num_material))
 zlevel(k)%cm_zlev(i,j)%fm_num_mat=0

 do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
   z=z0+ink-1
    
 do inj=1, zlevel(k)%cm_zlev(i,j)%fine(2)
      y=y0+inj-1
   do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
    x=x0+ini-1
    matone=hrt_data(1)%matrix(x,y,z) 
    if(matone .le. 0) then
  write(err_message, &
"('negtive or zero mat. number found (ERR1601) at: voxel(', 2(I0,1x),I0 ')=', I0 )" ) x, y, z, matone
        call TrapInputError(1)
    endif
    zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=matone    
       zlevel(k)%cm_zlev(i,j)%fm_num_mat(matone)=&
       zlevel(k)%cm_zlev(i,j)%fm_num_mat(matone)+1
       
    mat_fm(matone)=mat_fm(matone)+1
       mat_vol(matone)=mat_vol(matone)+zlevel(k)%cm_zlev(i,j)%fine_vol
       
   enddo  !end ini
   enddo     !end inj
 enddo       !end ink

zlevel(k)%cm_zlev(i,j)%cm_num_mat=0
do m=1, num_material
 if(zlevel(k)%cm_zlev(i,j)%fm_num_mat(m) .ne. 0) &
  zlevel(k)%cm_zlev(i,j)%cm_num_mat=zlevel(k)%cm_zlev(i,j)%cm_num_mat+1
enddo
if  (zlevel(k)%cm_zlev(i,j)%cm_num_mat .eq. 1) then
zlevel(k)%cm_zlev(i,j)%uni_zfine=0
else
zlevel(k)%cm_zlev(i,j)%uni_zfine=2
endif

enddo !end i
enddo !end j
enddo !end k

call OutputMatdis

return
end subroutine


subroutine DefineSrcMatrix
use heart
use paraset1
use paraset2
use paraset3
use paraset4
use ErrControl
use files


integer i,j,k,m
integer x0,y0,z0,x,y,z,ini,inj,ink

real srcden

allocate(s_mag_src(num_material))
s_mag_src=0

do i=1,num_src_mat
 s_mag_src(id_s_mat(i))=s_intensity(i)
enddo

do k=1, num_cmesh(3)
 do j=1,num_cmesh(2)
  do i=1, num_cmesh(1) 
   
   x0=0
   y0=0
   z0=0
  
  do m=1,k-1
    z0=z0+cgdim(3)%pixb(m)
  enddo
  z0=z0+1
 
  do m=1,j-1
    y0=y0+cgdim(2)%pixb(m)
  enddo
  y0=y0+1
 
  do m=1,i-1
    x0=x0+cgdim(1)%pixb(m)
  enddo
  x0=x0+1
  
  
 ! zlevel(k)%cm_zlev(i,j)%uni_zfine=2
 do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
   z=z0+ink-1
    
 do inj=1, zlevel(k)%cm_zlev(i,j)%fine(2)
      y=y0+inj-1
   do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
    x=x0+ini-1
    
       srcden=s_mag_src(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)) 
       if(srcden .ne. 0) then
        if(srcden .gt. 0) then
      zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)=srcden
     else        
       zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)=abs(srcden) &
       / mat_vol(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)) !* &
   ! zlevel(k)%cm_zlev(i,j)%fine_vol 
     endif
    endif
    
 if (num_binary .eq. 2) then  
    if(hrt_data(2)%list(hrt_data(2)%matrix(x,y,z)) .ge. 0) then
     zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)=zlevel(k)%cm_zlev(i,j)%&
     src_matrix(ini,inj,ink)+hrt_data(2)%list(hrt_data(2)%matrix(x,y,z))
 else
      write(warning_message, &
 "('Negtive source intensity found in file ', A, ' at voxel(', 2(I0,1x),I0 ')=',ES12.5, ' negtive value ignored' )") &
 trim(cur_filename), x,y,z, hrt_data(2)%list(hrt_data(2)%matrix(x,y,z))
      call TrapError(-6)
    endif
    endif
                                                    
    zlevel(k)%cm_zlev(i,j)%sum_act=zlevel(k)%cm_zlev(i,j)%sum_act+zlevel(k)%cm_zlev(i,j)&
         %src_matrix(ini,inj,ink)
       
   enddo  !end ini
   enddo     !end inj
 enddo       !end ink
 
 if(zlevel(k)%cm_zlev(i,j)%sum_act .gt. 0 ) &
   zlevel(k)%cm_zlev(i,j)%s_flag=1
    
 
enddo !end i
enddo !end j
enddo !end k

call ReadInSrc

return
end subroutine



subroutine InitMeshing
use heart
use funs
use paraset1
use paraset2
use paraset3
use paraset4

integer i,j,k,in
integer err

do i=1,3
  dist(i)=num_pix(i)*dist_pix(i)
enddo


do j=1,3
cgdim(j)%boundary(0)=0.0
do i=1,num_cmesh(j)
cgdim(j)%boundary(i)=cgdim(j)%boundary(i-1)+cgdim(j)%pixb(i)*dist_pix(j)
enddo
enddo


num_zlev=num_cmesh(3)
allocate(zlevel(num_zlev), z_lev_pos(num_zlev+1))
z_lev_pos=cgdim(3)%boundary


do k=1,num_zlev
 zlevel(k)%ncx=num_cmesh(1)
 zlevel(k)%ncy=num_cmesh(2)
 allocate(zlevel(k)%x_cm_pos(1+zlevel(k)%ncx))
 allocate(zlevel(k)%y_cm_pos(1+zlevel(k)%ncy))
 zlevel(k)%x_cm_pos=cgdim(1)%boundary
 zlevel(k)%y_cm_pos=cgdim(2)%boundary
enddo

!making xs file name
prb_xs%name=prbname
prb_xs%dir='.'
prb_xs%ext='xs'
prb_xs%fullname=''
prb_xs%fullname=Getfullname(prb_xs,0)

tot_num_fm=0
num_cm=num_cmesh(1)*num_cmesh(2)*num_cmesh(3)


do k=1, num_cmesh(3)
 allocate(zlevel(k)%cm_zlev(num_cmesh(1),num_cmesh(2))) 
 do j=1,num_cmesh(2)
   do i=1, num_cmesh(1) 
    zlevel(k)%cm_zlev(i,j)%fine(1)=cgdim(1)%pixb(i)
    zlevel(k)%cm_zlev(i,j)%fine(2)=cgdim(2)%pixb(j)
 zlevel(k)%cm_zlev(i,j)%fine(3)=cgdim(3)%pixb(k)
 zlevel(k)%cm_zlev(i,j)%finsize=dist_pix
 zlevel(k)%cm_zlev(i,j)%med=zlevel(k)%cm_zlev(i,j)%fine

   zlevel(k)%cm_zlev(i,j)%s_flag=0
   zlevel(k)%cm_zlev(i,j)%sum_act=0.0
   zlevel(k)%cm_zlev(i,j)%tot_fm=zlevel(k)%cm_zlev(i,j)%fine(1)*&
       zlevel(k)%cm_zlev(i,j)%fine(2)*zlevel(k)%cm_zlev(i,j)%fine(3)
   zlevel(k)%cm_zlev(i,j)%fine_vol=zlevel(k)%cm_zlev(i,j)%finsize(1)*zlevel(k)%cm_zlev(i,j)%&
                            finsize(2)*zlevel(k)%cm_zlev(i,j)%finsize(3)
   tot_num_fm=tot_num_fm+zlevel(k)%cm_zlev(i,j)%fine(1)*zlevel(k)%cm_zlev(i,j)%fine(2)&
            *zlevel(k)%cm_zlev(i,j)%fine(3)
    allocate(zlevel(k)%cm_zlev(i,j)%mat_matrix(zlevel(k)%cm_zlev(i,j)%fine(1),&
            zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)),STAT=err)
    if (err .ne. 0)  then 
  stop 'mat_matrix allocation error'
    else
     zlevel(k)%cm_zlev(i,j)%mat_matrix=0
 endif
    allocate(zlevel(k)%cm_zlev(i,j)%src_matrix(zlevel(k)%cm_zlev(i,j)%fine(1),&
            zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)),STAT=err)
    if (err .ne. 0)  then
  stop 'src_matrix allocation error'
    else 
  zlevel(k)%cm_zlev(i,j)%src_matrix=0
 endif

   allocate(zlevel(k)%cm_zlev(i,j)%cm_x( zlevel(k)%cm_zlev(i,j)%fine(1) ),&
             zlevel(k)%cm_zlev(i,j)%cm_y( zlevel(k)%cm_zlev(i,j)%fine(2) ),&
             zlevel(k)%cm_zlev(i,j)%cm_z( zlevel(k)%cm_zlev(i,j)%fine(3) ) )  
    
 do in=1,zlevel(k)%cm_zlev(i,j)%fine(1)
     zlevel(k)%cm_zlev(i,j)%cm_x(in)=zlevel(1)%x_cm_pos(i)+&
           (in-0.5)*zlevel(k)%cm_zlev(i,j)%finsize(1)
    enddo

    do in=1,zlevel(k)%cm_zlev(i,j)%fine(2)
     zlevel(k)%cm_zlev(i,j)%cm_y(in)=zlevel(1)%y_cm_pos(j)+&
           (in-0.5)*zlevel(k)%cm_zlev(i,j)%finsize(2)
    enddo

    do in=1,zlevel(k)%cm_zlev(i,j)%fine(3)
      zlevel(k)%cm_zlev(i,j)%cm_z(in)=z_lev_pos(k)+&
         (in-0.5)*zlevel(k)%cm_zlev(i,j)%finsize(3)
 enddo

   enddo
 enddo  
enddo

end subroutine

! assign other parameters for the heart model
subroutine AssignPenmshVaribles
use heart
use paraset1
use paraset2
use paraset3
use paraset4
use files
  
!integer i,j,k

ndmeth_global=2

!zlevel(1)%cm_zlev(1,1)%s_flag=1
!zlevel(1)%cm_zlev(1,1)%src_matrix=0
!zlevel(1)%cm_zlev(1,1)%sum_act=6.0E+10

!do k=1, zlevel(1)%cm_zlev(1,1)%fine(3)
!do j=1, zlevel(1)%cm_zlev(1,1)%fine(2)
!do i=1, zlevel(1)%cm_zlev(1,1)%fine(1)

!zlevel(1)%cm_zlev(1,1)%src_matrix(i,j,k)=1.923060E+07

!enddo
!enddo
!enddo

return
end subroutine 
