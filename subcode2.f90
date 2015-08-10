!subcode2.f90 : meshing and defining mat number on each fine mesh
!Subroutines included in this file:
!          Meshing:  Check the Gemetry and assign some gemetry varibles
!     AssignMatNum:  Assign material number for each fine mesh.     
!     WriteDotOut:   output mat dist. file .out, penmsh output file
 
subroutine Meshing
use paraset1
use paraset2
use paraset4
use ErrControl
use funs

integer i,j,k,in
integer ncm

tot_num_fm=0

num_cm_zlev=zlevel(1)%ncx*zlevel(1)%ncy
num_cm=num_cm_zlev*num_zlev
num_cmesh(1)=zlevel(z_start)%ncx
num_cmesh(2)=zlevel(z_start)%ncy
num_cmesh(3)=num_zlev

!Alocate space to calculate model mat vols
allocate(mat_vol(num_material),mat_fm(num_material))
mat_vol=0
mat_fm=0
total_vol=(z_lev_pos(num_zlev+1)-z_lev_pos(1))*&
   (zlevel(z_start)%y_cm_pos(1+zlevel(z_start)%ncy)-zlevel(z_start)%y_cm_pos(1))*&
   (zlevel(z_start)%x_cm_pos(1+zlevel(z_start)%ncx)-zlevel(z_start)%x_cm_pos(1))

if(IsZSingle .eq. 0) then
do k=2, num_zlev
 if(zlevel(k)%ncx .ne. zlevel(1)%ncx) then
  write(cur_filename,*) 'different ncx(# of cm along x) detected on:'
  write(err_message,"('zlev=',I0,' and ','zlev=',I0)") 1,k
  call TrapInputError(1)
 endif
 if(zlevel(k)%ncy .ne. zlevel(1)%ncy) then
  write(cur_filename,*) 'different ncy(# of cm along x) detected on:'
  write(err_message,"('zlev=',I0,' and ','zlev=',I0) ") 1,k
  call TrapInputError(1)
 endif
 if(real_compare(zlevel(1)%x_cm_pos,zlevel(k)%x_cm_pos,zlevel(1)%ncx) .eq. 0) then
   write(warning_message,&
"('different x coarse mesh position detected on zlev=',I0,' and ','zlev= ',I0,'zlev 1 data applied to all zlevs' )") 1,k
   call TrapInputError(0)
 endif
 if(real_compare(zlevel(1)%y_cm_pos,zlevel(k)%y_cm_pos,zlevel(1)%ncy) .eq. 0) then
   write(warning_message,& 
 "('different y coarse mesh position detected on zlev=',I0,' and ','zlev= ',I0,'zlev 1 data applied to all zlevs' )") 1,k
   call TrapInputError(0)
 endif
enddo
endif

do k=z_start, z_end
 do j=1,zlevel(k)%ncy
   do i=1, zlevel(k)%ncx   
    zlevel(k)%cm_zlev(i,j)%med(1)=zlevel(k)%cm_zlev(i,j)%fine(1)/ratio_x(k)
    zlevel(k)%cm_zlev(i,j)%med(2)=zlevel(k)%cm_zlev(i,j)%fine(2)/ratio_y(k)
    zlevel(k)%cm_zlev(i,j)%med(3)=zlevel(k)%cm_zlev(i,j)%fine(3)/ratio_z(k)
    tot_num_fm=tot_num_fm+zlevel(k)%cm_zlev(i,j)%fine(1)*zlevel(k)%cm_zlev(i,j)%fine(2)&
       *zlevel(k)%cm_zlev(i,j)%fine(3)
    ncm=(k-1)*num_cm_zlev+(j-1)*num_cmesh(1)+i

    zlevel(k)%cm_zlev(i,j)%finsize(1)=(zlevel(k)%x_cm_pos(i+1)-zlevel(k)%x_cm_pos(i))/zlevel(k)%cm_zlev(i,j)%fine(1)
    zlevel(k)%cm_zlev(i,j)%finsize(2)=(zlevel(k)%y_cm_pos(j+1)-zlevel(k)%y_cm_pos(j))/zlevel(k)%cm_zlev(i,j)%fine(2)
    zlevel(k)%cm_zlev(i,j)%finsize(3)=(z_lev_pos(k+1)-z_lev_pos(k))/zlevel(k)%cm_zlev(i,j)%fine(3)
    zlevel(k)%cm_zlev(i,j)%fine_vol=zlevel(k)%cm_zlev(i,j)%finsize(1)*zlevel(k)%cm_zlev(i,j)%&
                                 finsize(2)*zlevel(k)%cm_zlev(i,j)%finsize(3)
    zlevel(k)%cm_zlev(i,j)%cm_vol=(zlevel(k)%x_cm_pos(i+1)-zlevel(k)%x_cm_pos(i))*&
           (zlevel(k)%y_cm_pos(j+1)-zlevel(k)%y_cm_pos(j))*&
                (z_lev_pos(k+1)-z_lev_pos(k))
    
    zlevel(k)%cm_zlev(i,j)%tot_fm=zlevel(k)%cm_zlev(i,j)%fine(1)*&
       zlevel(k)%cm_zlev(i,j)%fine(2)*zlevel(k)%cm_zlev(i,j)%fine(3)
    
    allocate(zlevel(k)%cm_zlev(i,j)%cm_x( zlevel(k)%cm_zlev(i,j)%fine(1) ),&
         zlevel(k)%cm_zlev(i,j)%cm_y( zlevel(k)%cm_zlev(i,j)%fine(2) ),&
         zlevel(k)%cm_zlev(i,j)%cm_z( zlevel(k)%cm_zlev(i,j)%fine(3) ) )  
    
 do in=1,zlevel(k)%cm_zlev(i,j)%fine(1)
    zlevel(k)%cm_zlev(i,j)%cm_x(in)=zlevel(k)%x_cm_pos(i)+&
        (in-0.5)*zlevel(k)%cm_zlev(i,j)%finsize(1)
 enddo

 do in=1,zlevel(k)%cm_zlev(i,j)%fine(2)
   zlevel(k)%cm_zlev(i,j)%cm_y(in)=zlevel(k)%y_cm_pos(j)+&
       (in-0.5)*zlevel(k)%cm_zlev(i,j)%finsize(2)
 enddo

 do in=1,zlevel(k)%cm_zlev(i,j)%fine(3)
   zlevel(k)%cm_zlev(i,j)%cm_z(in)=z_lev_pos(k)+&
   (in-0.5)*zlevel(k)%cm_zlev(i,j)%finsize(3)
 enddo
 if(mod(zlevel(k)%cm_zlev(i,j)%fine(1),ratio_x(k)) .ne. 0) then 
!   write(*,*) 'fine to med ratio along x is not dividable on:' 
!      write(*,"('zlev=',I3,'CM= ',I5)") k, ncm
      zlevel(k)%cm_zlev(i,j)%med(1)=int(zlevel(k)%cm_zlev(i,j)%fine(1)/ratio_x(k))+1
 endif
 if(mod(zlevel(k)%cm_zlev(i,j)%fine(2),ratio_y(k)) .ne. 0) then 
!   write(*,*) 'fine to med ratio along y is not dividable on:' 
!      write(*,"('zlev=',I3,'CM=',I5)") k, ncm
   zlevel(k)%cm_zlev(i,j)%med(2)=int(zlevel(k)%cm_zlev(i,j)%fine(2)/ratio_y(k))+1
 endif
 if(mod(zlevel(k)%cm_zlev(i,j)%fine(2),ratio_z(k)) .ne. 0) then 
!      write(*,*) 'fine to med ratio along z is not dividable on:' 
!      write(*,"('zlev=',I3,'CM=',I5)") k, ncm
   zlevel(k)%cm_zlev(i,j)%med(3)=int(zlevel(k)%cm_zlev(i,j)%fine(3)/ratio_z(k))+1
   endif 
        
   enddo
 enddo  
enddo


end subroutine

! assign array in each coarse mesh
subroutine AssignMatNum
use paraset1
use paraset2
use paraset4
use ErrControl
use files
use funs
! integer filenum

integer i,j,k
! real*8 finsize(3), cm_vol, fine_vol 
!for cm3
integer cm_typ, num_reg
! integer i_var

allocate(mat_list(num_material))

do k=z_start, z_end
do j=1,zlevel(k)%ncy
do i=1, zlevel(k)%ncx

zlevel(k)%cm_zlev(i,j)%s_flag=0
cm_typ=mod(abs(zlevel(k)%cm_zlev(i,j)%cm_type),10)
num_reg=abs(zlevel(k)%cm_zlev(i,j)%num_subregion)

call InitCMmat(i,j,k,cm_typ,num_reg)


if(zlevel(k)%cm_zlev(i,j)%cm_type .le. -1 .or. zlevel(k)%cm_zlev(i,j)%cm_type .gt. 10) then
 call SetupCMmat(i,j,k) 
else
 call FinalizeCMmat(i,j,k,num_reg)
endif

enddo
enddo
enddo

deallocate(mat_list)

if(IsZSingle .eq. 0) call OutputMatdis

end subroutine

!assign source magnitude for each fine mesh
subroutine AssignSrcMag
use paraset2
use paraset4
use ErrControl
use paraset1, only : max_src, z_start, z_end !, IsZsingle
use files

integer i,j,k,ini,inj,ink
!integer err
real srcden

! num_source=0

!uniform src distributed in one or some material
do k=z_start, z_end
 
 do j=1,zlevel(k)%ncy
     
  do i=1, zlevel(k)%ncx

!allocate(zlevel(k)%cm_zlev(i,j)%src_matrix(zlevel(k)%cm_zlev(i,j)%fine(1),&
!    zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)),STAT=err)
!if (err .gt. 0)  stop 'src_matrix allocation error'
!zlevel(k)%cm_zlev(i,j)%src_matrix=0    

zlevel(k)%cm_zlev(i,j)%sum_act=0
if(zlevel(k)%cm_zlev(i,j)%s_flag .ge. 1) then
! num_source=num_source+1
! allocate(zlevel(k)%cm_zlev(i,j)%src_matrix(zlevel(k)%cm_zlev(i,j)%fine(1),&
!          zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)),STAT=err)
! if (err .gt. 0)  stop 'mat_matrix allocation error'
! zlevel(k)%cm_zlev(i,j)%src_matrix=0
 if(zlevel(k)%cm_zlev(i,j)%cm_num_mat .ne. 1 ) then
  do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
   do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
    do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
      srcden=s_mag_src(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)) 
     if(srcden .ne. 0) then
      if(srcden .gt. 0) then
   zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)=srcden
   else        
      zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)=abs(srcden) &
      / mat_vol(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)) !* &
      ! zlevel(k)%cm_zlev(i,j)%fine_vol 
      endif
      zlevel(k)%cm_zlev(i,j)%sum_act=zlevel(k)%cm_zlev(i,j)%sum_act+&
         zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)
     endif
      
      enddo !ini
    enddo !inj
  enddo !ink
 else !only one mat in this cm
 srcden=s_mag_src(zlevel(k)%cm_zlev(i,j)%mat_matrix(1,1,1))
 if(srcden .ge. 0) then
 zlevel(k)%cm_zlev(i,j)%src_matrix=srcden
 else
 zlevel(k)%cm_zlev(i,j)%src_matrix=  &
      -srcden      &
      /mat_vol(zlevel(k)%cm_zlev(i,j)%cm_mat_num) !* &
 ! zlevel(k)%cm_zlev(i,j)%fine_vol
 endif
  zlevel(k)%cm_zlev(i,j)%sum_act= &   
 zlevel(k)%cm_zlev(i,j)%src_matrix(1,1,1)*zlevel(k)%cm_zlev(i,j)%tot_fm
 !   /mat_vol(zlevel(k)%cm_zlev(i,j)%cm_mat_num)* &
    
 endif !cm_type
 endif !s_flag=1     
enddo  !i
enddo  !j 
enddo  !k

!read in src input
call ReadInSrc


!get the maxium src in all fine meshes for plotting src
max_src=0
do k=1, num_zlev
 do j=1,zlevel(k)%ncy
  do i=1, zlevel(k)%ncx
  if(max_src .lt. maxval(zlevel(k)%cm_zlev(i,j)%src_matrix) )&
    max_src=maxval(zlevel(k)%cm_zlev(i,j)%src_matrix)
enddo
enddo
enddo


end subroutine


!output .out files (penmsh output material layout file)
subroutine WriteDotOut
use paraset1
use paraset2
use paraset4
use funs
use files

integer i,j,k
integer index_cm,k_fine,j_fine,i_fine

allocate(matfile(num_cmesh(3)))
do k=1,num_cmesh(3)
  write(form,"( '(','A',I0, ',I0,A4', ')' )" ) prbname_len
  write(matfile(k)%fullname, form ) prbname, k, '.out'  
  matfile(k)%dir=penmsh_inp%dir
  open(31,file=matfile(k)%dir//trim(matfile(k)%fullname))
  do j=1,num_cmesh(2)
    do i=1, num_cmesh(1)
      index_cm=num_cm_zlev*(k-1)+(j-1)*num_cmesh(1)+i
      write(1,"('**nmattp=  ',I4)") index_cm
      do k_fine=1,zlevel(k)%cm_zlev(i,j)%fine(3)
    
       write(1,"('***fine zmesh= ',I3)") k_fine
       do j_fine=1,zlevel(k)%cm_zlev(i,j)%fine(2)
         write(form,"( '(', I0, '(I0,1X)', ')' )" ) zlevel(k)%cm_zlev(i,j)%fine(1)
         write(1, form) & 
            (zlevel(k)%cm_zlev(i,j)%mat_matrix(i_fine,j_fine,k_fine),i_fine=&
            1,zlevel(k)%cm_zlev(i,j)%fine(1) )
       enddo
      enddo
    
     enddo
  enddo
  close(31)
enddo

end subroutine

!finalize mat_vol, mat_fm, s_flag without .src for non-overlay cm
subroutine FinalizeCMmat(i,j,k,num_reg)
use paraset4
use paraset2

integer i,j,k,num_reg

integer ini,inj
integer matn


if(num_reg .ne. 1) then
zlevel(k)%cm_zlev(i,j)%uni_zfine=1
do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
 matn=zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)
 mat_fm(matn)=&
   mat_fm(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1))+zlevel(k)%cm_zlev(i,j)%fine(3)
 mat_vol(matn)=&
  mat_vol(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1))+zlevel(k)%cm_zlev(i,j)%fine(3)*&
  zlevel(k)%cm_zlev(i,j)%fine_vol
 zlevel(k)%cm_zlev(i,j)%fm_num_mat(matn)=&
  zlevel(k)%cm_zlev(i,j)%fm_num_mat(matn)+zlevel(k)%cm_zlev(i,j)%fine(3)

enddo
enddo


else
!handle only one region
 matn=zlevel(k)%cm_zlev(i,j)%mat_matrix(1,1,1)
 mat_fm(matn)=&
   mat_fm(matn)+zlevel(k)%cm_zlev(i,j)%tot_fm
 mat_vol(matn)=&
   mat_vol(matn)+zlevel(k)%cm_zlev(i,j)%cm_vol
 zlevel(k)%cm_zlev(i,j)%fm_num_mat(matn)=&
  zlevel(k)%cm_zlev(i,j)%tot_fm
endif

!setup zlevel(k)%cm_zlev(i,j)%cm_num_mat
zlevel(k)%cm_zlev(i,j)%cm_num_mat=0
do matn=1, num_material
 if(zlevel(k)%cm_zlev(i,j)%fm_num_mat(matn) .ne. 0)  &
  zlevel(k)%cm_zlev(i,j)%cm_num_mat=zlevel(k)%cm_zlev(i,j)%cm_num_mat+1 
enddo


!setup s_flag before .src
if(s_format .lt. 0 ) then

do matn=1, num_material
 if(zlevel(k)%cm_zlev(i,j)%fm_num_mat(matn)*s_mag_src(matn) .ne. 0)  &
  zlevel(k)%cm_zlev(i,j)%s_flag=zlevel(k)%cm_zlev(i,j)%s_flag+1 
enddo

endif  !s_format


return
end subroutine


subroutine InitCMmat(i,j,k,cm_typ,num_reg)
use paraset4
use paraset2
use ErrControl
use files
use funs
use constants
integer i,j,k,cm_typ,num_reg

integer err
integer reg
integer ini,inj,ink,cm_num_1z
real x1,y1, r1, y_bot,y_top

allocate(zlevel(k)%cm_zlev(i,j)%mat_matrix(zlevel(k)%cm_zlev(i,j)%fine(1),&
           zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)),STAT=err)
if (err .gt. 0)  then
 write(err_message,"( 'mat_matrix allocation error' )" )
 call TrapMemError(1)
endif
zlevel(k)%cm_zlev(i,j)%mat_matrix=0


allocate(zlevel(k)%cm_zlev(i,j)%src_matrix(zlevel(k)%cm_zlev(i,j)%fine(1),&
    zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)),STAT=err)
if (err .gt. 0)  then
 write(err_message,"( 'src_matrix allocation error' )" )
 call TrapMemError(1)
endif
zlevel(k)%cm_zlev(i,j)%src_matrix=0    


allocate( zlevel(k)%cm_zlev(i,j)%fm_num_mat(num_material) )
zlevel(k)%cm_zlev(i,j)%fm_num_mat=0
cm_num_1z=i+(j-1)*zlevel(k)%ncx
select case(cm_typ)

case(1)
if(num_reg .ne. 1) then
do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
 if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then 
   y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)
 else
   y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)-zlevel(k)%y_cm_pos(j) 
 endif
 do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
 if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then
   x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)
 else !relative coordinates
   x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)-zlevel(k)%x_cm_pos(i)
 endif
 subloop1 : do reg=1, num_reg
 if( (x1-zlevel(k)%cm_zlev(i,j)%r_lft(reg) )*(x1-zlevel(k)%cm_zlev(i,j)%r_rgt(reg) )&
  .le. 0 .and. (y1-zlevel(k)%cm_zlev(i,j)%theta_bot(reg))*(y1-zlevel(k)%cm_zlev(i,j)%theta_top(reg)&
   ) .le. 0 ) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)=zlevel(k)%cm_zlev(i,j)%cm_mat(reg)
  exit subloop1
  endif
 enddo subloop1
  
if(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1) .eq. 0) then
write(warning_message,&
"('CM ',I0,1x,'on zlev ',I0,1x,':fm(i,j,:): ', 2(I0,1x),' not covered for CM type 1')") cm_num_1z,k,ini,inj
   call TrapError(-4)
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)=zlevel(k)%cm_zlev(i,j)%cm_mat(num_reg)   
endif 

enddo !ini
enddo !inj 

do ink=2,zlevel(k)%cm_zlev(i,j)%fine(3)
 zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,ink)=&
 zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,1)
enddo

else  !only one region
zlevel(k)%cm_zlev(i,j)%mat_matrix=zlevel(k)%cm_zlev(i,j)%cm_mat_num

endif 

case(2)
if(num_reg .ne. 1) then
  do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
  if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then 
    y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)
  else
   y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)-zlevel(k)%y_cm_pos(j) 
  endif
   do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
   if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then
     x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)
   else !relative coordinates
    x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)-zlevel(k)%x_cm_pos(i)
   endif
    r1=sqrt(x1**2+y1**2)
    subloop2: do reg=1, num_reg
    if( (r1-zlevel(k)%cm_zlev(i,j)%r_lft(reg) )*(r1-zlevel(k)%cm_zlev(i,j)%r_rgt(reg) )&
     .le. 0 .and. (y1-zlevel(k)%cm_zlev(i,j)%theta_bot(reg))*(y1-zlevel(k)%cm_zlev(i,j)%theta_top(reg))&
    .le. 0 ) then
     zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)=zlevel(k)%cm_zlev(i,j)%cm_mat(reg)
     exit subloop2
     endif
    enddo subloop2
  
  if(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1) .eq. 0) then
   write(warning_message,&
"('CM ',I0,1x,'on zlev ',I0,1x,':fm(i,j,:): ', 2(I0,1x),'  not covered for CM type 2')") cm_num_1z,k,ini,inj
   call TrapError(-4)
   zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)=zlevel(k)%cm_zlev(i,j)%cm_mat(num_reg)   
  endif 

 enddo !ini
enddo !inj 

do ink=2,zlevel(k)%cm_zlev(i,j)%fine(3)
 zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,ink)=&
 zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,1)
enddo

else  !only one region
zlevel(k)%cm_zlev(i,j)%mat_matrix=zlevel(k)%cm_zlev(i,j)%cm_mat_num

endif 


case(3)
if(num_reg .ne. 1) then
 do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
 if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then 
    y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)
  else
   y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)-zlevel(k)%y_cm_pos(j) 
  endif
 
 do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
  if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then
     x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)
   else !relative coordinates
    x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)-zlevel(k)%x_cm_pos(i)
  endif
  
  r1=sqrt(x1**2+y1**2)
  subloop3: do reg=1, num_reg
  if(zlevel(k)%cm_zlev(i,j)%theta_bot(reg) .eq. 0) then
   if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then
    y_bot=zlevel(k)%y_cm_pos(j)
   else
     y_bot=0
   endif
  else
   y_bot=x1*tan(pi/180.0*zlevel(k)%cm_zlev(i,j)%theta_bot(reg))
  endif
  if(zlevel(k)%cm_zlev(i,j)%theta_top(reg) .eq. 0 .or. zlevel(k)%cm_zlev(i,j)%theta_top(reg) .eq. 90 ) then
   if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then
    y_top=zlevel(k)%y_cm_pos(j+1)
   else
    y_top=zlevel(k)%y_cm_pos(j+1)-zlevel(k)%y_cm_pos(j)
   endif
   
  else
   y_top=x1*tan(pi/180.0*zlevel(k)%cm_zlev(i,j)%theta_top(reg))
  endif
  if( (r1-zlevel(k)%cm_zlev(i,j)%r_lft(reg) )*(r1-zlevel(k)%cm_zlev(i,j)%r_rgt(reg) )&
     .le. 0 .and. (y1-y_bot)*(y1-y_top) .le. 0 ) then
   zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)=&
  zlevel(k)%cm_zlev(i,j)%cm_mat(reg)
   exit subloop3
  endif
    
  enddo subloop3
  
 if(zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1) .eq. 0) then
  write(warning_message,&
"('CM ',I0,1x,'on zlev ',I0,1x,':fm(i,j,:): ', 2(I0,1x),'  not covered for CM type 3')") cm_num_1z,k,ini,inj   
  call TrapError(-4)
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,1)=zlevel(k)%cm_zlev(i,j)%cm_mat(num_reg)   
 endif  
 enddo !ini
enddo !inj 

do ink=2,zlevel(k)%cm_zlev(i,j)%fine(3)
   zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,ink)=&
   zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,1)
enddo

else  !one region
zlevel(k)%cm_zlev(i,j)%mat_matrix=zlevel(k)%cm_zlev(i,j)%cm_mat_num
endif

case default
 zlevel(k)%cm_zlev(i,j)%mat_matrix=zlevel(k)%cm_zlev(i,j)%cm_mat_num

end select

!free a little memory
if(num_reg .ne. 1) then
 deallocate(zlevel(k)%cm_zlev(i,j)%cm_mat, zlevel(k)%cm_zlev(i,j)%r_rgt,&
 zlevel(k)%cm_zlev(i,j)%r_lft,zlevel(k)%cm_zlev(i,j)%theta_bot,& 
 zlevel(k)%cm_zlev(i,j)%theta_top )
endif

return
end subroutine


!handle overlay in one CM
subroutine SetupCMmat(i,j,k)
use paraset1
use paraset2
use paraset4
use ErrControl
use files
use funs
use constants
integer i,j,k

type qcircle
 integer mt1, mt2
end type
type tpoint
 real x,y
end type
integer m, matn
integer ini,inj,ink
real x1,y1,z1
real d, sec_r1, sec_r2

integer :: overlay_num_abs=0,zfin_flag=1
integer matloop

integer:: mo_typ=0
type(qcircle) :: qc(4)
type(tpoint):: qt(4)
!for 31 32 33 34
integer tri(3,4),tri_idx
integer block, num_block
integer inii,injj,inkk
integer rep_i,rep_j,rep_k, my_mat
integer cm_num_zlev
logical :: IsZthr=.TRUE.

data qc/qcircle(1,3),qcircle(2,3),qcircle(2,4),qcircle(1,4)/
data tri/3,2,4, &
         4,1,3,&
  1,2,4,&
  2,1,3/
mat_list=0  

cm_num_zlev=num_cmesh(1)*(j-1)+i
!num_block=int(-zlevel(k)%cm_zlev(i,j)%cm_type/10)

num_block=zlevel(k)%cm_zlev(i,j)%num_block
if(num_block .eq. 0) num_block=1
zlevel(k)%cm_zlev(i,j)%uni_zfine=1
zfin_flag=0

loop_ink: do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then 
  z1=zlevel(k)%cm_zlev(i,j)%cm_z(ink) 
else
  z1=zlevel(k)%cm_zlev(i,j)%cm_z(ink)-z_lev_pos(k) 
endif

!z1=zlevel(k)%cm_zlev(i,j)%cm_z(ink) 
do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then 
  y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)
else
  y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)-zlevel(k)%y_cm_pos(j) 
endif
!y1=zlevel(k)%cm_zlev(i,j)%cm_y(inj)
do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)

if(zlevel(k)%cm_zlev(i,j)%num_subregion .gt. 0) then
  x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)
else !relative coordinates
  x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)-zlevel(k)%x_cm_pos(i)
endif
!x1=zlevel(k)%cm_zlev(i,j)%cm_x(ini)

do block=1, num_block
 overlay_num_abs=abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_num)
 do m=1, overlay_num_abs

  repk_loop: do rep_k=1,zlevel(k)%cm_zlev(i,j)%overlay_block(block)%num_rep(3,m)
  do rep_j=1,zlevel(k)%cm_zlev(i,j)%overlay_block(block)%num_rep(2,m)
  do rep_i=1,zlevel(k)%cm_zlev(i,j)%overlay_block(block)%num_rep(1,m)
my_mat=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_mat(m)
if(my_mat .lt. 0 .and. zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_typ(m) .lt. 0) then
 my_mat=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%mat_rep(m)%d3_int(rep_i, rep_j, rep_k)
endif
mo_typ=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_typ(m)
mo_typ=abs(mo_typ)
overlay_type: select case (mo_typ)

!retagular 
case (1)

 Ax=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+&
     (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Ay=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)+&
     (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Bx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+&
      (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 By=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)+&
      (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .ge. 6) then
 zfin_flag=1
 Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 if((z1 .ge. Zx)  .and. (z1 .le. Zy) ) then
 if( (x1 .ge. Ax) .and. (x1 .le. Ay) ) then
 if( (y1 .ge. Bx) .and. (y1 .le. By) ) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
 endif
      
 endif 
 endif
 else
if( (x1 .ge. Ax) .and. (x1 .le. Ay) ) then
 if( (y1 .ge. Bx) .and. (y1 .le. By) ) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
 endif
 endif 
 endif   
                  
! cilinder along z
case (2)
 center_x=(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+zlevel(k)% &
 cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m))/2+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m) 
 center_y=(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+zlevel(k)% &
 cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m))/2+&
 (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
  
 r=abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)-zlevel(k)%cm_zlev(i,j)%&
   overlay_block(block)%overlay_bon(1,m))/2
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .eq. 6) then
 zfin_flag=1
 Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
  (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 if((z1 .ge. Zx)  .and. (z1 .le. Zy) ) then
  d=((x1-center_x)**2+(y1-center_y)**2)**0.5    
 if( d .le. r) then
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 exit repk_loop
 end if
 endif
 else
  d=((x1-center_x)**2+(y1-center_y)**2)**0.5    
  if( d .le. r) then
   zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
   exit repk_loop
  end if
 endif

! cilinder along y
case (25)
 zfin_flag=1
 center_x=(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+zlevel(k)% &
 cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m))/2+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m) 
 center_z=(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+zlevel(k)% &
 cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m))/2+&
 (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
  
 r=abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)-zlevel(k)%cm_zlev(i,j)%&
   overlay_block(block)%overlay_bon(1,m))/2
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .eq. 6) then

 Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
  (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 if((y1 .ge. Zx)  .and. (y1 .le. Zy) ) then
  d=((x1-center_x)**2+(z1-center_z)**2)**0.5    
 if( d .le. r) then
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 exit repk_loop
 end if
 endif
 else
  d=((x1-center_x)**2+(z1-center_z)**2)**0.5    
  if( d .le. r) then
   zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
   exit repk_loop
  end if
 endif
! cilinder along x
case (26)
 zfin_flag=1
 center_y=(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+zlevel(k)% &
 cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m))/2+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m) 
 center_z=(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+zlevel(k)% &
 cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m))/2+&
 (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
  
 r=abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)-zlevel(k)%cm_zlev(i,j)%&
   overlay_block(block)%overlay_bon(1,m))/2
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .eq. 6) then

 Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
  (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 if((x1 .ge. Zx)  .and. (x1 .le. Zy) ) then
  d=((y1-center_y)**2+(z1-center_z)**2)**0.5    
 if( d .le. r) then
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 exit repk_loop
 end if
 endif
 else
  d=((y1-center_y)**2+(z1-center_z)**2)**0.5    
  if( d .le. r) then
   zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
   exit repk_loop
  end if
 endif      
!triangle
case(3)
 Ax=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Ay=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)+&
  (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 Bx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 By=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)+&
  (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 Cx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Cy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
  (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
!Barycentric Coordinates test if a point is inside a triangle
 b0=((Bx-Ax)*(Cy-Ay)-(Cx-Ax)*(By-Ay))
 if(b0 .eq. 0) then 
  write(err_message,"('points superposition for triangle overlay',I0,' in cm:  ',I3)") m,cm_num_zlev 
  cur_filename=inputfile(k)%fullname
  call TrapInputError(1)
 endif
    
 b1=((Bx-x1)*(Cy-y1)-(Cx-x1)*(By-y1))/b0
 b2=((Cx-x1)*(Ay-y1)-(Ax-x1)*(Cy-y1))/b0
 b3=((Ax-x1)*(By-y1)-(Bx-x1)*(Ay-y1))/b0
          
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .ge. 8) then
  zfin_flag=1
  Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(7,m)+&
   (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
  Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(8,m)+&
   (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
  
  if((z1 .ge. Zx) .and. (z1 .le. Zy) ) then
  if( (b1 .ge. 0) .and. (b2 .ge. 0) .and. (b3 .ge. 0)) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
  end if
  endif
  else 
  if( (b1 .ge. 0) .and. (b2 .ge. 0) .and. (b3 .ge. 0)) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
  end if
  endif
     
!sphere
case (4)
 zfin_flag=1
 center_x=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 center_y=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)+&
  (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 center_z=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+&
  (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m)
 r=abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m))
 d=((x1-center_x)**2+(y1-center_y)**2+(z1-center_z)**2)**0.5
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .ge. 6) then
 Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m)  
 if((z1 .ge. Zx)  .and. (z1 .le. Zy) ) then
 if(d .le. r ) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 end if
 endif
else
if( d .le. r) then
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 exit repk_loop
end if 
endif
!cone
!case (5)
!zfin_flag=1
!x1=x1+(rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
!y1=y1+(rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
!z1=z1+(rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m)
!h=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)
!d=abs(z1-zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m))+&
! abs(z1-zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+h)
!if(d .eq. abs(h)) then
! r=( (x1-zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m))**2 + &
!  (y1-zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m))**2 + &
!   (z1-zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m))**2 )**0.5
! theta=acos( -1*h*(z1-zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m))/(abs(h)*r) )     
! if (theta .le. atan(abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)/h)) ) then
!  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
!  exit repk_loop
! endif
! endif

!sector
case (51)
 Ax=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+&
   (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Ay=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)+&
   (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 Bx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+&
   (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 By=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)+&
   (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 Cx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
   (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Cy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
   (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
!Barycentric Coordinates test if a point is inside a triangle
  b0=((Bx-Ax)*(Cy-Ay)-(Cx-Ax)*(By-Ay))
 if(b0 .eq. 0) then 
  write(err_message,"('points superposition for triangle overlay in cm:  ',I3)")  cm_num_zlev 
  cur_filename=inputfile(k)%fullname
  call TrapInputError(1)
 endif
 sec_r1=((Bx-Ax)**2+(By-Ay)**2)**0.5
 sec_r2=((Cx-Ax)**2+(Cy-Ay)**2)**0.5
 r=(sec_r1+sec_r2)/2
 if(abs(sec_r1-sec_r2)/sec_r1 .gt. 0.1) then
  write(err_message,&
"('point 2&3 are not on a circle centered at point 1 for sector overlay ', I0,' in cm ',I0)")  m,cm_num_zlev 
  cur_filename=inputfile(k)%fullname
  call TrapInputError(1)
 endif 

b1=((Bx-x1)*(Cy-y1)-(Cx-x1)*(By-y1))/b0
b2=((Cx-x1)*(Ay-y1)-(Ax-x1)*(Cy-y1))/b0
b3=((Ax-x1)*(By-y1)-(Bx-x1)*(Ay-y1))/b0
if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .eq. 8) then
 zfin_flag=1
 Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(7,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(8,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
 if((z1 .ge. Zx) .and. (z1 .le. Zy) ) then
 d=((x1-Ax)**2+(y1-Ay)**2)**0.5
 if( (d .le. r) .and. (b2 .ge. 0) .and. (b3 .ge. 0)) then
   zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
   exit repk_loop
 end if
 endif
else 
d=((x1-Ax)**2+(y1-Ay)**2)**0.5
if( (d .le. r) .and. (b2 .ge. 0) .and. (b3 .ge. 0)) then
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 exit repk_loop
end if
endif

!hexagon by center point and one vertex
case(6)

center_x=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)
center_y=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)

vertex_x=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)
vertex_y=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)

r_out=(center_x-vertex_x)**2+(center_y-vertex_y)**2
r_out=sqrt(r_out)
if(r_out .lt. 10E-7) then
   write(err_message,&
   "('Hexgon vertex and center points too close for Overlay', I0,' in cm ',I0)")  m,cm_num_zlev 
   cur_filename=inputfile(k)%fullname
  call TrapInputError(1)
endif

r_in=r_out*cos30 !sqrt(3)/2


d=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)*cos30
if(center_y .eq. vertex_y) then
 center_x=center_x+(rep_i-1)*d
 vertex_x=vertex_x+(rep_i-1)*d
 center_y=center_y+(rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 vertex_y=vertex_y+(rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 if(mod(rep_i-1,2) .eq. 1) then
  center_y=center_y-0.5*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
  vertex_y=vertex_y-0.5*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 endif

elseif(center_x .eq. vertex_x) then
 center_y=center_y+(rep_j-1)*d
 vertex_y=vertex_y+(rep_j-1)*d

 center_x=center_x+(rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 vertex_x=vertex_x+(rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 if(mod(rep_j-1,2) .eq. 1) then
  center_x=center_x-0.5*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
  vertex_x=vertex_x-0.5*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 endif
else
   write(err_message,&
   "('Hexgon vertex has to be on X or Y axis for Overlay', I0,' in cm ',I0)")  m,cm_num_zlev 
   cur_filename=inputfile(k)%fullname
  call TrapInputError(1)
endif

IsZthr=.TRUE.
if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .ge. 6) then
   zfin_flag=1
   Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
   Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
   if((z1 .lt. Zx) .or. (z1 .gt. Zy) ) IsZthr=.FALSE.
endif

if(IsZthr) then   
    
d=(x1-center_x)**2+(y1-center_y)**2
d=sqrt(d)

if(d .le. r_in) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
endif

if(d .le. r_out) then
! theta=atan2(vertex_y-center_y, vertex_x-center_x)
! gamma=atan2(y1-center_y, x1-center_x)
! beta=gamma-theta
! if(beta .lt. 0) beta=beta+2*pi
! beta=beta-(int(3*beta/pi)-1)*pi/3
 
! r=r_out*cos30/sin(beta)
! if(d .le. r) then
!  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
!  exit repk_loop
! endif
cos_theta=(x1-center_x)*(vertex_x-center_x)+(y1-center_y)*(vertex_y-center_y)
cos_theta=cos_theta/(d*r_out)
sin_theta=sqrt(1-cos_theta*cos_theta)
if(cos_theta .gt. 0.5) then
  dk=cos_theta*cos30+sin_theta*0.5
  dk=r_in/dk
else if( cos_theta .gt. -0.5) then
  dk=r_in/sin_theta
else
  dk=0.5*sin_theta-cos30*cos_theta
  dk=r_in/dk
endif
if(d .le. dk) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
endif
endif

endif


!quarter circle
case (21,22,23,24)
 Ax=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+&
     (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Ay=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)+&
     (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 Bx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+&
      (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 By=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)+&
      (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 
 r=abs(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)- &
  zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m) )

 center_x=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(qc(mo_typ-20)%mt1,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 center_y=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(qc(mo_typ-20)%mt2,m)+&
 (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .ge. 6) then
  zfin_flag=1
  Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(5,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
  Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(6,m)+&
    (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
  if((z1 .ge. Zx) .and. (z1 .le. Zy) ) then
  if( (x1 .ge. Ax) .and. (x1 .le. Ay) ) then
  if( (y1 .ge. Bx) .and. (y1 .le. By) ) then
  d=((x1-center_x)**2+(y1-center_y)**2)**0.5    
  if( d .le. r) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
  endif
  endif
      
endif 
endif
else
if( (x1 .ge. Ax) .and. (x1 .le. Ay) ) then
if( (y1 .ge. Bx) .and. (y1 .le. By) ) then
 d=((x1-center_x)**2+(y1-center_y)**2)**0.5    
 if( d .le. r) then
 zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
 exit repk_loop
 endif
 endif
    
endif
     
endif   

!quarter triangle
case (31,32,33,34)
 qt(1)%x=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(2,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 qt(2)%x=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(1,m)+&
  (rep_i-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(1,m)
 qt(3)%x=qt(2)%x
 qt(4)%x=qt(1)%x
  
 qt(1)%y=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(4,m)+&
  (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 qt(2)%y=qt(1)%y
 qt(3)%y=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(3,m)+&
  (rep_j-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(2,m)
 qt(4)%y=qt(3)%y
 
 tri_idx=mo_typ-30
 Ax=qt(tri(1,tri_idx))%x
 Ay=qt(tri(1,tri_idx))%y
 Bx=qt(tri(2,tri_idx))%x
 By=qt(tri(2,tri_idx))%y
 Cx=qt(tri(3,tri_idx))%x
 Cy=qt(tri(3,tri_idx))%y
! Barycentric Coordinates test if a point is inside a triangle
 b0=((Bx-Ax)*(Cy-Ay)-(Cx-Ax)*(By-Ay))
 if(b0 .eq. 0) then 
  write(err_message,"('points superposition for triangle overlay',I0,' in cm:  ',I3)") m,cm_num_zlev 
  cur_filename=inputfile(k)%fullname
  call TrapInputError(1)
 endif
    
 b1=((Bx-x1)*(Cy-y1)-(Cx-x1)*(By-y1))/b0
 b2=((Cx-x1)*(Ay-y1)-(Ax-x1)*(Cy-y1))/b0
 b3=((Ax-x1)*(By-y1)-(Bx-x1)*(Ay-y1))/b0
          
 if(zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(maxos,m) .ge. 8) then
  zfin_flag=1
  Zx=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(7,m)+&
   (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
  Zy=zlevel(k)%cm_zlev(i,j)%overlay_block(block)%overlay_bon(8,m)+&
   (rep_k-1)*zlevel(k)%cm_zlev(i,j)%overlay_block(block)%dist_rep(3,m) 
  
  if((z1 .ge. Zx) .and. (z1 .le. Zy) ) then
  if( (b1 .ge. 0) .and. (b2 .ge. 0) .and. (b3 .ge. 0)) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
  end if
  endif
  else 
  if( (b1 .ge. 0) .and. (b2 .ge. 0) .and. (b3 .ge. 0)) then
  zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)=my_mat
  exit repk_loop
  end if
  endif
     
!put your overlay here


 
case default overlay_type
stop 'overlay type not supported'  
end select  overlay_type

enddo !rep_i
enddo !rep_j
enddo repk_loop

enddo  !end m
enddo !block

matn=zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)
mat_vol(matn)=mat_vol(matn)+zlevel(k)%cm_zlev(i,j)%fine_vol
mat_list(matn)=mat_list(matn)+1


enddo   !ini  
enddo     !inj


if(ink .eq. 1) then
zlevel(k)%cm_zlev(i,j)%uni_zfine=zfin_flag+1
if(zfin_flag .eq. 0 ) then
 do inkk=2,zlevel(k)%cm_zlev(i,j)%fine(3)
 zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,inkk)= &
  zlevel(k)%cm_zlev(i,j)%mat_matrix(:,:,1) 
 enddo !inkk
 do injj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
 do inii=1,zlevel(k)%cm_zlev(i,j)%fine(1)
 matn=zlevel(k)%cm_zlev(i,j)%mat_matrix(inii,injj,1)
 mat_vol(matn)=mat_vol(matn)+&
  zlevel(k)%cm_zlev(i,j)%fine_vol*(zlevel(k)%cm_zlev(i,j)%fine(3)-1)
 mat_list(matn)=mat_list(matn)+(zlevel(k)%cm_zlev(i,j)%fine(3)-1)
 enddo !endinii
 enddo !endinjj
 exit loop_ink
endif
endif

enddo loop_ink 

     
!define s_flag and num of mat in cm: number of src mat in the CM
zlevel(k)%cm_zlev(i,j)%cm_num_mat=0
!   zlevel(k)%cm_zlev(i,j)%s_flag=0
do matloop=1, num_material
if (mat_list(matloop) .ge. 1) then
 zlevel(k)%cm_zlev(i,j)%cm_num_mat=zlevel(k)%cm_zlev(i,j)%cm_num_mat+1
 mat_fm(matloop)=mat_fm(matloop)+mat_list(matloop)  
if(s_format .lt. 0 ) then
if(s_mag_src(matloop) .ne. 0) zlevel(k)%cm_zlev(i,j)%s_flag=&
  zlevel(k)%cm_zlev(i,j)%s_flag+1 
endif
endif
enddo !matloop
zlevel(k)%cm_zlev(i,j)%fm_num_mat=mat_list
deallocate(zlevel(k)%cm_zlev(i,j)%overlay_block)

return
end subroutine


