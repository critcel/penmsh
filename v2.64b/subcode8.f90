!output mat balance file prbname_out.mba file
 
subroutine OutputMatdis
use paraset1
use paraset2
use paraset4
use ErrControl
use files
use funs


integer,dimension(:),allocatable :: g
logical ex
character*120 a_buffer
integer*8 :: mat_fm_total=0

integer i,j,k,m,n
integer ini,inj,ink
real :: xyz(3)=0, vol=0

total_mat_vol=0

do i=1,num_material
 total_mat_vol=total_mat_vol+mat_vol(i)
enddo

allocate(mat_tgt(num_material),g(num_material))
g=num_material+1

allocate (cen_mat(3,num_material) )
cen_mat=0
!calculate mass center for each material
do k=1, num_cmesh(3)
do j=1, num_cmesh(2)
do i=1, num_cmesh(1)

do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
  xyz(3)=zlevel(k)%cm_zlev(i,j)%cm_z(ink) 
do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
  xyz(2)=zlevel(k)%cm_zlev(i,j)%cm_y(inj)
do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
  xyz(1)=zlevel(k)%cm_zlev(i,j)%cm_x(ini)

m=zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)
vol=zlevel(k)%cm_zlev(i,j)%fine_vol
do n=1, 3
cen_mat(n,m)=cen_mat(n,m)+xyz(n)*vol
enddo

enddo !ini
enddo !inj
enddo !ink


enddo !i
enddo !j
enddo !k

do m=1, num_material
if(mat_vol(m) .ne. 0) then
cen_mat(:,m)=cen_mat(:,m)/mat_vol(m)
do i=1,3
if(cen_mat(i,m) .lt. 1e-7 ) cen_mat(i,m)=0.0
enddo
endif
enddo

cur_filename=trim(prbname)//'_out.mba'
 
open(unit=51, file=cur_filename)


write(51,"('Material Balance Summary for Model: ',A)") trim(prbname)
write(51,*) 
write(51,"('Section 1 : Coarse Mesh infomation')")
write(51,"(A,A)") ' CM #   x_size   y_size   z_size        fm num      tot_fm          ', &
                  'fm size             fm vol  num_mat    # of fm per mat'
!"(' CM #',3x, 'x_size', 3x,'y_size', 3x,'z_size',8x,'fm num',5x,'tot_fm', 4x, 'fm vol',2x, 'num_mat',4x,'# of fm per mat')")
write(form,"('(37x,A1,5x,A1,5x,A1,15x,A1,7x,A1,7x,A1,24x,A3,',I0,'(I3,5x))' )" ) num_material
write(51,form) 'x','y','z', 'x', 'y','z', 'mat',  &
     (m,m=1,num_material)
do k=1, num_zlev
do j=1,zlevel(k)%ncy
do i=1, zlevel(k)%ncx 
write(form,"('(I5,1x,3(f8.3,1x),3(I5,1x),I7,1x,3(f7.2,1x),2x, ES9.3,3x,I2,6x,',I0,'(I7,1x))' )") num_material
write(51,form) & 
 (k-1)*num_cm_zlev+(j-1)*num_cmesh(1)+i,zlevel(k)%x_cm_pos(i+1)-zlevel(k)%x_cm_pos(i),&
 zlevel(k)%y_cm_pos(j+1)-zlevel(k)%y_cm_pos(j),z_lev_pos(k+1)-z_lev_pos(k),&
 zlevel(k)%cm_zlev(i,j)%fine(1),zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3),&
 zlevel(k)%cm_zlev(i,j)%tot_fm, zlevel(k)%cm_zlev(i,j)%finsize, zlevel(k)%cm_zlev(i,j)%fine_vol,&
 zlevel(k)%cm_zlev(i,j)%cm_num_mat,(zlevel(k)%cm_zlev(i,j)%fm_num_mat(m),m=1,num_material) 
enddo
enddo
enddo

k=0

write(51,*)
write(51,"('Section 2 : Model mat. infomation')")
write(51,*)
 
curfile%name=prbname
curfile%ext='mba'
curfile%fullname=''
curfile%dir=input_dir
cur_filename=Getfullname(curfile,1)

inquire(file=cur_filename,exist=ex)
if(ex ) then

write(51,"('Contents of file: ',A)") trim(cur_filename)
write(51,*)

open(unit=44,file=cur_filename)
write(warning_message,"('Reading Mat vol file:' ,A)") trim(cur_filename)
call DisplayMsg

do i=1,6
 read(44,'(A)',err=1001,end=1001) a_buffer
 write(51,*) trim(a_buffer)
enddo

do i=1,num_material
 read(44,*,err=1001,end=1001) mat_tgt(i)%nam,mat_tgt(i)%num,mat_tgt(i)%vol_tgt,mat_tgt(i)%den
 backspace(44)

 read(44,'(A)') a_buffer
 write(51,*) trim(a_buffer)
 mat_tgt(i)%vol_mod=mat_vol(mat_tgt(i)%num)  
 g(mat_tgt(i)%num)=i
 k=k+1
 tgt_tot_vol=tgt_tot_vol+mat_tgt(i)%vol_tgt
 tgt_tot_mass=tgt_tot_mass+mat_tgt(i)%vol_tgt*mat_tgt(i)%den
 mod_tot_mass=mod_tot_mass+mat_tgt(i)%vol_mod*mat_tgt(i)%den

enddo
goto 1020 
1001 write(51,"('read error(first 6 lines are comments): ',A)") trim(cur_filename)
 
1020 close(44)

if( k .eq. num_material) then !read success
write(warning_message,"('Reading Mat vol file done, file:' ,A)") trim(cur_filename)
call DisplayMsg

else !read error, skip the file
write(warning_message,"('read error, skiped file: ',A)") trim(cur_filename)
call TrapReadError(0)
endif

else

write(51,"('input file does not exist: ',A)") trim(cur_filename)
write(51,"('comparison skipped ')")

do i=1,num_material
 mat_tgt(i)%vol_mod=mat_vol(i)  
 tgt_tot_vol=tgt_tot_vol+mat_tgt(i)%vol_tgt
 tgt_tot_mass=tgt_tot_mass+mat_tgt(i)%vol_tgt*mat_tgt(i)%den
enddo
endif !read in .mba done 

write(51,*)
write(51,"('Model Material Inventory/Volume Table')")
write(51,*)
if(k .eq. num_material) then
write(51,"('Targeted Mean Density: ',ES14.6,' g/cm3 from .mba file')") &
    tgt_tot_mass/tgt_tot_vol 
write(51,"('     Targeted Volume : ',ES14.6,' cm3    Mass:  ', ES14.6, ' g')") &
   tgt_tot_vol,tgt_tot_mass
write(51,"('        Model Volume : ',ES14.6,' cm3    Mass:  ', ES14.6, ' g')") &
   total_mat_vol,mod_tot_mass
    
else
 write(51,"('        Model Volume : ',ES14.6,' cm3  in ',I0,1x  ' total fine meshes')") &
  total_mat_vol,tot_num_fm
 write(51,*)
endif
write(51,*)
write(51,"(' Material    Mat      Fine        Model (cm3)      % of      Mdl/Trg     Model Mass    Model    ' ) ")
write(51,"('   name      No.     Meshes        volume          Total     V-Ratio     Excess (g)    Mass (g) ' )  ")
write(51,"('----------  ------  --------     -------------    -------    ---------  ------------  ---------- ')")

do i=1,num_material
mat_fm_total=mat_fm(i)+mat_fm_total
if(g(i) .ne. num_material+1) then
write(51,"(A10, 1x I5,4x, I8, 3x, ES14.6,2x, f9.3, 5x, ES9.3,3x,ES10.3, 3x, ES10.3)") &
 mat_tgt(g(i))%nam, mat_tgt(g(i))%num, mat_fm(i),mat_tgt(g(i))%vol_mod, &
 mat_tgt(g(i))%vol_mod*100/total_mat_vol, mat_tgt(g(i))%vol_mod/mat_tgt(g(i))%vol_tgt, &
 (mat_tgt(g(i))%vol_mod-mat_tgt(g(i))%vol_tgt)*mat_tgt(g(i))%den,&
 mat_tgt(g(i))%vol_mod*mat_tgt(g(i))%den
else
write(51,"(A10, 1x I5,4x, I8,3x, ES14.6,2x, f9.3)") &
'unknown',i, mat_fm(i),mat_vol(i),mat_vol(i)*100/total_mat_vol
endif

enddo

write(51,*)
write(51,"(A10, 1x, 9x, I8,3x,ES14.6,3x)") &
   'TOTAL',  mat_fm_total,total_mat_vol

write(51,*)
write(51, "('Material mass center ')")
write(51,"(' Material    Mat           mass center coodinates    ' ) ")
write(51,"('   name      No.           x              y             z     ' )  ")
write(51,"('----------  -----    ------------   ------------   ------------ ')")

do i=1, num_material
if(g(i) .ne. num_material+1) then
write(51,"(A10, 2x I5, 4x, ES12.5, 3x, ES12.5, 3x ES12.5)") &
 mat_tgt(g(i))%nam, i, cen_mat(:,i)
else
write(51,"(A10, 2x I5, 4x, ES12.5, 3x, ES12.5, 3x ES12.5)") &
 'unknown', i, cen_mat(:,i)
endif
enddo
   


write(51,*)
write(51, "('Material Vol in Cm^3')")
write(51,*)
write(51,"('Total Vol    :',2x, ES14.6)") total_vol
write(51,"('Total Mat Vol:',2x, ES14.6)") total_mat_vol

close(51)

write(warning_message,"('Mat vol comparison output file done, file: ' ,A)") &
    trim(prbname)//'_out.mba'
call DisplayMsg

return
end subroutine


!output file : flux summary
!format
! x, y, z, vol, mat#, fixsrc, flx-g1, flx-g2, .....
subroutine MxOutputFlx
use paraset1
use paraset2
use paraset3
use paraset4
use files
use funs
use ErrControl

implicit none
integer, parameter ::  num_pre=6
real :: point(num_pre)=0.0

integer g,m,n
integer cmi, cmj,cmk, fmi, fmj, fmk
integer i,j
integer :: num_vol=0 , fm=0 , fm_bon=0

character(LEN=100) filename_flux
!integer, allocatable :: fm_bon(:)

if(num_flx_out .eq. 0) then
!write one binary file
   curfile%name=prbname
   curfile%dir='.'
   if(IsFluxout .eq. 2) then
    write(curfile%ext, "('adj.bin')" ) 
   else 
    write(curfile%ext, "('fwd.bin')") 
   endif
   curfile%fullname=''

   cur_filename=Getfullname(curfile, 0)
   cur_fileunit=123

   open (unit=cur_fileunit, file=cur_filename, FORM='UNFORMATTED', ACCESS='STREAM', err=1001)

do cmk=1, num_cmesh(3)
do cmj=1, num_cmesh(2)
do cmi=1, num_cmesh(1)
 
 point(4)=zlevel(cmk)%cm_zlev(cmi,cmj)%fine_vol

 do fmk=1, zlevel(cmk)%cm_zlev(cmi,cmj)%fine(3)
   point(3)=zlevel(cmk)%cm_zlev(cmi,cmj)%cm_z(fmk)
 
 do fmj=1, zlevel(cmk)%cm_zlev(cmi,cmj)%fine(2)
   point(2)=zlevel(cmk)%cm_zlev(cmi,cmj)%cm_y(fmj)

 do fmi=1, zlevel(cmk)%cm_zlev(cmi,cmj)%fine(1)
   point(1)=zlevel(cmk)%cm_zlev(cmi,cmj)%cm_x(fmi)
   
   point(5)=zlevel(cmk)%cm_zlev(cmi,cmj)%mat_matrix(fmi,fmj,fmk)
  if(s_format .eq. 2 .and. IsFluxout .eq. 1) then
   point(6)=0
   do g=1, num_group
     point(6)=point(6) + zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk)      
   enddo
  else
   point(6)=zlevel(cmk)%cm_zlev(cmi,cmj)%src_matrix(fmi,fmj,fmk)
  endif

   if (IsFluxout .eq. 2) then   !adjoint flux reverse group order
   write(cur_fileunit)  point, &
    (zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk), g=num_group, 1, -1) 
   else
    write(cur_fileunit)  point, &
    (zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk), g=1, num_group)
   endif
  enddo
 enddo
 enddo
enddo
enddo
enddo
if(IsFluxout .eq. 1) then
   write(*, "('flux output file done (binary) ',  A)")   trim(cur_filename)
   write(READLOG,"('flux output file done (binary) ', A)")  trim(cur_filename)
else
  write(*,"('adjoint flux output file done (binary ', A, ' GROUP ORDER FLIPPED')")  m,trim(cur_filename)
  write(READLOG,"('adjoint flux output file done (volumn ', I0, ') :', A, ' GROUP ORDER FLIPPED')")  m,trim(cur_filename)
endif
close(cur_fileunit)
return
endif

!write multi-vol ASCII file, max vol=99
if(num_flx_out .gt. 99) num_flx_out=99
if(num_flx_out .lt. 0) num_flx_out=abs(num_flx_out)
if(num_flx_out .gt. tot_num_fm) num_flx_out=tot_num_fm

i=int(tot_num_fm/num_flx_out)
j=tot_num_fm - num_flx_out*i

fm=0
m=0
fm_bon=0
do cmk=1, num_cmesh(3)
do cmj=1, num_cmesh(2)
do cmi=1, num_cmesh(1)
 
 point(4)=zlevel(cmk)%cm_zlev(cmi,cmj)%fine_vol

 do fmk=1, zlevel(cmk)%cm_zlev(cmi,cmj)%fine(3)
   point(3)=zlevel(cmk)%cm_zlev(cmi,cmj)%cm_z(fmk)
 
 do fmj=1, zlevel(cmk)%cm_zlev(cmi,cmj)%fine(2)
   point(2)=zlevel(cmk)%cm_zlev(cmi,cmj)%cm_y(fmj)

 do fmi=1, zlevel(cmk)%cm_zlev(cmi,cmj)%fine(1)
   point(1)=zlevel(cmk)%cm_zlev(cmi,cmj)%cm_x(fmi)
   
   point(5)=zlevel(cmk)%cm_zlev(cmi,cmj)%mat_matrix(fmi,fmj,fmk)
  if(s_format .eq. 2 .and. IsFluxout .eq. 1) then
   point(6)=0
   do g=1, num_group
     point(6)=point(6) + zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk)      
   enddo
  else
   point(6)=zlevel(cmk)%cm_zlev(cmi,cmj)%src_matrix(fmi,fmj,fmk)
  endif
   fm=fm+1
!multi-vol files
   if(fm .gt. fm_bon) then
   if(m .ne. 0) then
    if(IsFluxout .eq. 1) then
      write(*, "('flux output file done (volumn ', I0, '): ', A)")  m, trim(cur_filename)
      write(READLOG,"('flux output file done (volumn ', I0, '): ', A)") m, trim(cur_filename)
    else
      write(*,"('adjoint flux output file done (volumn ', I0, '): ', A, ' GROUP ORDER FLIPPED')")  m,trim(cur_filename)
      write(READLOG,"('adjoint flux output file done (volumn ', I0, '): ', A, ' GROUP ORDER FLIPPED')")  m,trim(cur_filename)
    endif
    close(cur_fileunit)
   endif
   m=m+1
   fm_bon=fm_bon+i
   if(m .le. j) fm_bon=fm_bon+1
   curfile%name=prbname
   curfile%dir='.'
   if(IsFluxout .eq. 2) then
    write(curfile%ext, "('adj.m', I1, I1)" ) int(m/10), mod(m,10)
   else 
    write(curfile%ext, "('fwd.m', I1, I1)") int(m/10), mod(m,10)
   endif
   curfile%fullname=''

   cur_filename=Getfullname(curfile, 0)
   cur_fileunit=123

   open (unit=cur_fileunit, file=cur_filename, err=1001)
   write(form, "('Flux distribution file:', A, ' generated at ' )" ) trim(cur_filename)
   call TimeStamp(cur_fileunit)
   backspace(cur_fileunit)  !reverse to the blank line in TimeStamp
   write(cur_fileunit,"('Volumn ', I0 ' of ',I0, ': fine mesh ', I0, ' to ',I0, ' of total ', I0)") &
     m, num_flx_out, fm, fm_bon, tot_num_fm
!call intlen(num_group, len)
   write(form, "( '(1x, A1,14x, A1,14x, A1,14x, A3,12x, A4,11x, A6,9x,' I0,'(I4, 9x) )' )") &
     num_group
   write(cur_fileunit, form) 'x', 'y', 'z' , 'vol', 'mat#', 'fixsrc',  (g, g=1, num_group)

   write(form, "('(',I0,'(ES12.5,3x)', I0, '(ES12.5,1x) )' )" ) num_pre, num_group
   endif  !fm_bon
   
   if (IsFluxout .eq. 2) then   !adjoint flux reverse group order
   write(cur_fileunit, form)  point, &
    (zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk), g=num_group, 1, -1) 
   else
    write(cur_fileunit, form)  point, &
    (zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk), g=1, num_group) 
   endif
 enddo
 enddo
 enddo
enddo
enddo
enddo

    if(IsFluxout .eq. 1) then
      write(*, "('flux output file done (volumn ', I0, '): ', A)")  m, trim(cur_filename)
      write(READLOG,"('flux output file done (volumn ', I0, '): ', A)") m, trim(cur_filename)
    else
      write(*,"('adjoint flux output file done (volumn ', I0, '): ', A, ' GROUP ORDER FLIPPED')")  m,trim(cur_filename)
      write(READLOG,"('adjoint flux output file done (volumn ', I0, '): ', A, ' GROUP ORDER FLIPPED')")  m,trim(cur_filename)
    endif

close(cur_fileunit)


return 

1001  write(err_message,*) 'open file error: no permission for creating file'
      call TrapOpenError(1)
end subroutine
