! Read in flux info from pentran output
! subroutines in this file(subcode7.f90)
!          AllocateFluxArray:   allocate memory for flux
!          ReadInFluxset:  read flux info
!          NormalizeFlux:  normalize a fluxset

subroutine AllocateFluxArray
use files
use ErrControl
use paraset1
use paraset4

integer i,j,k,set
integer sv
! integer name_len, dir_len, ext_len,fullname_len

allocate(grp_flag(num_group,1+IsRefFlux))
grp_flag=0

allocate(glb_fluxset(1+IsRefFlux))
do set=1, 1+IsRefFlux
  
  allocate(glb_fluxset(set)%max_glb_flux_grp(num_group), &
           glb_fluxset(set)%avg_glb_flux_grp(num_group))
  glb_fluxset(set)%max_glb_flux_grp=0
  glb_fluxset(set)%avg_glb_flux_grp=0
  
!  glb_fluxset(set)%fluxfile_name=prbname
!  if(IsPrbDotFlx(set) .eq. 1) glb_fluxset(set)%fluxfile_name='prb'
enddo
if(IsProcessCurr .eq. 0) then
  glb_fluxset(1)%fluxfile_dir=trim(flux_dir)
  if (IsRefFlux .ge. 1) then
   glb_fluxset(2)%fluxfile_dir=trim(ref_flux_dir)
  endif
else 
  glb_fluxset(1)%fluxfile_dir=trim(curr_dir)
endif
  
do k=1, num_zlev 
  do j=1,zlevel(k)%ncy
   do i=1, zlevel(k)%ncx
    
      allocate(zlevel(k)%cm_zlev(i,j)%fluxset(1+IsRefFLux))
    
      do set=1, 1+IsRefFlux
      allocate(zlevel(k)%cm_zlev(i,j)%fluxset(set)%max_loc_flux_grp(num_group),&
        zlevel(k)%cm_zlev(i,j)%fluxset(set)%avg_loc_flux_grp(num_group),&
        zlevel(k)%cm_zlev(i,j)%fluxset(set)%locg_set(num_group) ,STAT=sv)
          zlevel(k)%cm_zlev(i,j)%fluxset(set)%max_loc_flux_grp=0
        zlevel(k)%cm_zlev(i,j)%fluxset(set)%avg_loc_flux_grp=0
!        allocate(zlevel(k)%cm_zlev(i,j)%fluxset(set)%loc_flux(zlevel(k)%cm_zlev(i,j)%fine(1),&
!          zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3),&
!          num_group), STAT=sv)
          if(sv .ne.0) then
           write(cur_filename,"('Error when trying allocating memory for fluxset ',I0)") set
           write(err_message, *) 'Maybe not enough memory'
           cur_fileunit=0
           call TrapMemError(1)
          endif
         !for curr
!         allocate(zlevel(k)%cm_zlev(i,j)%fluxset(set)%loc_curr(zlevel(k)%cm_zlev(i,j)%fine(1),&
!                  zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3),&
!                          num_group), STAT=sv)
!          if(sv .ne.0) then
!           write(cur_filename,"('Error when trying allocating memory for flux current set ',I0)") set
!           write(err_message, *) 'Maybe not enough memory'
!           cur_fileunit=0
!           call TrapMemError(1)
!          endif
        
        enddo  !loop set
   
   enddo  !loop i
  enddo  !loop j
enddo  !loop k


end subroutine

subroutine ReadInFluxset(set_number)
use files
use funs
use ErrControl
use paraset1
use paraset2
use paraset4
use lineread
integer set_number

integer grp
integer i,j,k,ini,inj,ink, m
integer sv,lncount
character*80  flxname1, flxname2
logical ex1,ex2
real temp, tempbuf(2)
prbname_len=len_trim(prbname)

curfile%ext='flx'
curfile%fullname=''
curfile%dir=glb_fluxset(set_number)%fluxfile_dir


do grp=1, num_group
 
!looking for prb#.flx or prbname#.flx
write(form,"('(A',I0,',I0)')") prbname_len 
write(curfile%name,form ) trim(prbname),grp
flxname1=Getfullname(curfile,1)

inquire(file=flxname1,exist=ex1)

if(ex1 ) then
 cur_filename=flxname1
else 
 write(curfile%name,"('prb',I0)" ) grp
 flxname2=Getfullname(curfile,1)
 inquire(file=flxname2,exist=ex2)
 if(ex2 ) then
  cur_filename=flxname2
 else
  cur_filename=trim(flxname1)//' or '//trim(flxname2)
  write(warning_message,"('Warning:', A, ': not exist')") &
          trim(cur_filename)
  call DisplayMsg
  write(warning_message,"('   Grp ',I0, ' skipped')")  grp
   call DisplayMsg
   grp_flag(grp,set_number)=0
   grp_flag(grp,1)=0
  cycle
 endif
endif

if(set_number .eq. 2 .and. grp_flag(grp,1) .eq. 0) then
  grp_flag(grp,2)=0
  write(warning_message,"(' Second set (-f2) flux Reading Skipped for group ', I0, ' : ', A)")  grp, trim(cur_filename) 
  call DisplayMsg
  cycle
endif

cur_fileunit=grp+FLUXUNIT
open(unit=cur_fileunit,file=cur_filename , &
      status='OLD'  , err=7101) 
 
write(warning_message,"(' Reading Flux file (first line skipped) : ', A)")   trim(cur_filename) 
call DisplayMsg
lncount=1
read(cur_fileunit,*,err=7102,end=7103)

flag(1)=-1
if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
  if(flag(1) .ne. 1) then
   write(warning_message,"('Warning: ', A)")   trim(cur_filename) 
   call DisplayMsg
       write(warning_message,"(3x, I0, ' values found in each line, last column value is used' )") flag(1)
       call DisplayMsg
 endif
endif
! flag(1)=0

 do k=1, num_zlev 
  do j=1,zlevel(k)%ncy
   do i=1, zlevel(k)%ncx
     !allocate arrays
     allocate(zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_flux(zlevel(k)%cm_zlev(i,j)%fine(1),&
         zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)), STAT=sv)
     if(sv .ne.0) then
       write(cur_filename,"('Error when trying allocating memory for fluxset ',I0)") set_number
       write(err_message, *) 'Maybe not enough memory'
       cur_fileunit=0
       call TrapMemError(1)
     endif
     
     do ink=1,zlevel(k)%cm_zlev(i,j)%fine(3)
       do inj=1 , zlevel(k)%cm_zlev(i,j)%fine(2)
         do ini=1 , zlevel(k)%cm_zlev(i,j)%fine(1)
         lncount=lncount+1   
!        read(grp,*, err=503) tmpx,tmpy,tmpz,zlevel(k)%cm_zlev(i,j)%flux(ini,inj,ink,grp)
         read(cur_fileunit,*,end=7103, err=7102) (temp, m=1, flag(1)-1), tempbuf(1)
	 
	 if(set_number .eq. 2) then
	  tempbuf(2)=zlevel(k)%cm_zlev(i,j)%fluxset(1)%locg_set(grp)%locg_flux(ini,inj,ink)
	  if(tempbuf(1) .ne. 0) then
	    tempbuf(1)=(tempbuf(2)-tempbuf(1))/tempbuf(1)
          elseif(tempbuf(2) .ne. 0) then
	    tempbuf(1)=(tempbuf(2)-tempbuf(1))/tempbuf(2)
	  endif
	 endif
	  zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_flux(ini,inj,ink)=tempbuf(1)
         
          zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)=&
              zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)+tempbuf(1)
                 
          if(tempbuf(1) .gt. zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp) ) then
            zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp)=tempbuf(1)
          endif
	     enddo !end ini
           enddo !end inj
         enddo !ink
     !local avg flux by group
       zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)=&
         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)/(&
         zlevel(k)%cm_zlev(i,j)%fine(3)*zlevel(k)%cm_zlev(i,j)%fine(2)*&
                  zlevel(k)%cm_zlev(i,j)%fine(1))
     !add up global avg flux by grp
       glb_fluxset(set_number)%avg_glb_flux_grp(grp)=&
         glb_fluxset(set_number)%avg_glb_flux_grp(grp)+&
         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)*&
         zlevel(k)%cm_zlev(i,j)%cm_vol 

     if(zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp) .gt. &
         glb_fluxset(set_number)%max_glb_flux_grp(grp)) then
                     
       glb_fluxset(set_number)%max_glb_flux_grp(grp)=&
         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp)
     endif
   enddo !i
  enddo !j
 enddo !k
 close(cur_fileunit)
 
 grp_flag(grp, set_number)=1

 glb_fluxset(set_number)%avg_glb_flux_grp(grp)=&
            glb_fluxset(set_number)%avg_glb_flux_grp(grp)/total_vol
 
 glb_fluxset(set_number)%avg_glb_flux=glb_fluxset(set_number)%avg_glb_flux+&
                                 glb_fluxset(set_number)%avg_glb_flux_grp(grp)
 
 if(glb_fluxset(set_number)%max_glb_flux_grp(grp) .gt. glb_fluxset(set_number)%max_glb_flux) then
   glb_fluxset(set_number)%max_glb_flux=glb_fluxset(set_number)%max_glb_flux_grp(grp)
 endif  

enddo !grp

num_group_eff=sum(grp_flag(:,set_number))
if(num_group_eff .ne. 0) then
  glb_fluxset(set_number)%avg_glb_flux=glb_fluxset(set_number)%avg_glb_flux/num_group_eff
endif

if( IsNomalizeFlux .eq. 1) call NomalizeFlux(set_number)

7100 return 

7101 write(err_message,*) 'flux file open error'
     call TrapOpenError(1)

7102 write(err_message,"('flux file read error on line: ', I0) ") lncount 
     Call TrapReadError(1)
7103 write(err_message,*) 'flux file EOF (end of file)'
     Call TrapReadError(1)

return
end subroutine

!for current file: .fjn
subroutine ReadInCurrset(set_number)
use files
use ErrControl
use paraset1
use paraset2
use paraset4
use funs
use lineread
integer set_number

integer grp, sv
integer i,j,k,ini,inj,ink ,m
character*80  flxname1, flxname2
logical ex1,ex2
real temp, tempbuf(5)

prbname_len=len_trim(prbname)
!len_fluxname=len_trim(prbname)

curfile%ext='fjn'
curfile%fullname=''
curfile%dir=curr_dir

do grp=1, num_group
 
!looking for prb#.fjn or prbname#.fjn
write(form,"('(A',I0,',I0)')") prbname_len 
write(curfile%name,form ) trim(prbname),grp
flxname1=Getfullname(curfile,1)

inquire(file=flxname1,exist=ex1)

if(ex1 ) then
 cur_filename=flxname1
else 
 write(curfile%name,"('prb',I0)" ) grp
 flxname2=Getfullname(curfile,1)
 inquire(file=flxname2,exist=ex2)
 if(ex2 ) then
  cur_filename=flxname2
 else
  cur_filename=trim(flxname1)//' or '//trim(flxname2)
  write(warning_message,"('Warning:', A, ': not exist')") &
          trim(cur_filename)
  call DisplayMsg
  write(warning_message,"('   Grp ',I0, ' skipped')")  grp
   call DisplayMsg
  grp_flag(grp,1)=0
  cycle
 endif
endif

 cur_fileunit=grp+FLUXUNIT
 open(unit=cur_fileunit,file=cur_filename , &
      status='OLD'  , err=7101) 
 write(warning_message,"(' Reading Current file(first line skipped) : ', A)")   trim(cur_filename) 
 call DisplayMsg
 read(cur_fileunit,*,err=7102,end=7102)

 flag(1)=-1
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) then 
   if(flag(1) .lt. 5) then
    write(err_message,&
"(' ERR8000: only ',I0,' values found, ',I0,' expected in each line:', ' phi, J, Jx, Jy, Jz')") flag(1), 5
    call TrapInputError(1)
   elseif( flag(1) .gt. 5) then
     write(warning_message,"(3x, I0, ' values found in each line, last 5 rows value are used' )") flag(1)
     call DisplayMsg
  endif
 endif

 do k=1, num_zlev 
  do j=1,zlevel(k)%ncy
   do i=1, zlevel(k)%ncx
!allocate flux array
    allocate(zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_flux(zlevel(k)%cm_zlev(i,j)%fine(1),&
         zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)), STAT=sv)
     if(sv .ne.0) then
       write(cur_filename,"('Error when trying allocating memory for fluxset ',I0)") set_number
       write(err_message, *) 'Maybe not enough memory'
       cur_fileunit=0
       call TrapMemError(1)
     endif
!allocate current array
    allocate(zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_curr(zlevel(k)%cm_zlev(i,j)%fine(1),&
         zlevel(k)%cm_zlev(i,j)%fine(2),zlevel(k)%cm_zlev(i,j)%fine(3)), STAT=sv)
     if(sv .ne.0) then
       write(cur_filename,"('Error when trying allocating memory for currset ',I0)") set_number
       write(err_message, *) 'Maybe not enough memory'
       cur_fileunit=0
       call TrapMemError(1)
     endif
         do ink=1,zlevel(k)%cm_zlev(i,j)%fine(3)
           do inj=1 , zlevel(k)%cm_zlev(i,j)%fine(2)
            do ini=1 , zlevel(k)%cm_zlev(i,j)%fine(1)
            
! read(grp,*, err=503) tmpx,tmpy,tmpz,zlevel(k)%cm_zlev(i,j)%flux(ini,inj,ink,grp)
         read(cur_fileunit,*,end=7102, err=7102) (temp, m=1, flag(1)-1), tempbuf
!          zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%loc_flux(ini,inj,ink,grp),&
! for curr
!         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%loc_curr(ini,inj,ink,grp)%fcurr
         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_flux(ini,inj,ink)=tempbuf(1)
         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_curr(ini,inj,ink)%fcurr=tempbuf(2:5)

	 zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)=&
           zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)+tempbuf(1)
                 
         if(tempbuf(1) .gt. zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp) ) then
            zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp)=tempbuf(1)
                 
         endif

                enddo !end ini
           enddo !end inj
         enddo !ink
     !local avg flux by group
         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)=&
              zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)/(&
          zlevel(k)%cm_zlev(i,j)%fine(3)*zlevel(k)%cm_zlev(i,j)%fine(2)*&
                  zlevel(k)%cm_zlev(i,j)%fine(1))
     !add up global avg flux by grp
         glb_fluxset(set_number)%avg_glb_flux_grp(grp)=&
            glb_fluxset(set_number)%avg_glb_flux_grp(grp)+&
                zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%avg_loc_flux_grp(grp)*&
        zlevel(k)%cm_zlev(i,j)%cm_vol 

     if(zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp) .gt. &
                      glb_fluxset(set_number)%max_glb_flux_grp(grp)) then
                     
                         glb_fluxset(set_number)%max_glb_flux_grp(grp)=&
                         zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%max_loc_flux_grp(grp)
         endif
   enddo !i
  enddo !j
 enddo !k
 grp_flag(grp, set_number)=1

 close(cur_fileunit)
 
 glb_fluxset(set_number)%avg_glb_flux_grp(grp)=&
            glb_fluxset(set_number)%avg_glb_flux_grp(grp)/total_vol
 
 glb_fluxset(set_number)%avg_glb_flux=glb_fluxset(set_number)%avg_glb_flux+&
                                 glb_fluxset(set_number)%avg_glb_flux_grp(grp)
 
 if(glb_fluxset(set_number)%max_glb_flux_grp(grp) .gt. glb_fluxset(set_number)%max_glb_flux) then
   glb_fluxset(set_number)%max_glb_flux=glb_fluxset(set_number)%max_glb_flux_grp(grp)
 endif  

enddo !grp

num_group_eff=sum(grp_flag(:,1))
if(num_group_eff .ne. 0) then
glb_fluxset(set_number)%avg_glb_flux=glb_fluxset(set_number)%avg_glb_flux/num_group_eff
endif

7100 return 

7101 write(err_message,*) '.fjn file open error'
     call TrapOpenError(1)

7102 write(err_message,*) '.fjn file read error'
     Call TrapReadError(1)


return
end subroutine

! Nomalize a fluxset to max flux=1
Subroutine NomalizeFlux(set_number)
use paraset1
use paraset2
use paraset3
use paraset4
use files
use funs
use ErrControl
use mIn4deck
use lineread
use constants

integer set_number

integer i,j,k
integer ini,inj,ink
integer grp
logical ex


if(set_number .eq. 1) then
  allocate (n_factor_grp(num_group) )
  n_factor_grp=n_factor
endif
if(n_factor .lt. tiny .and. n_factor .gt. -tiny ) then
  n_factor_grp=glb_fluxset(set_number)%max_glb_flux
endif

if(set_number .eq. 1 .and. n_factor .lt. -tiny ) then

cur_filename=trim(nfile%fullname)
cur_fileunit=564
inquire(file=cur_filename, exist=ex)
if( .not. ex) then
  write(warning_message,"('Warning:', A, ': not exist, using default global n_factor=1.0')") &
          trim(cur_filename)
  call DisplayMsg
  n_factor_grp=glb_fluxset(set_number)%max_glb_flux
  n_factor=0
else

open(unit=cur_fileunit, file=cur_filename) 
cur_var='group normalization factor'

num_need=num_group
if(allocated(fidobank)) deallocate(fidobank)
allocate(fidobank(num_need))
num_read=0
lineseg=''

do while (num_read .lt. num_need)
 if(NumCmtLine(cur_fileunit,flag) .ge. 0) &
   read(cur_fileunit,"(A)",err=1002,end=1002) lineseg

 call FeedFido

enddo  !do while

do grp=1, num_group
 n_factor_grp(grp)=fidobank(grp)
 if(n_factor_grp(grp) .le. 0) then
    write(err_message,*) trim(cur_var)//': group normalization factor has to be positive '
    call TrapInputError(1)
 endif
enddo

deallocate(fidobank)
write(READLOG,"('group normalization factor (flux=flux/factor) from file: ', A)") trim(cur_filename)
write(READLOG,"(5(ES12.5,2x) )") n_factor_grp
endif  ! .ex. 

endif  !  set_number .eq. 1 and n_factor .lt. 0



do k=1, num_zlev 
 do j=1,zlevel(k)%ncy
    do i=1, zlevel(k)%ncx
      
      do grp=1, num_group
	     if(grp_flag(grp, set_number) .ne. 1) cycle 
          do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
             do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
           do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
                  
              zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_flux(ini,inj,ink)=&
                  zlevel(k)%cm_zlev(i,j)%fluxset(set_number)%locg_set(grp)%locg_flux(ini,inj,ink)/&
                   n_factor_grp(grp)
 !                   glb_fluxset(set_number)%max_glb_flux

                   enddo !ini
             enddo  !inj
          enddo !ink
      enddo !grp
        enddo
  enddo
enddo
return
1002  call TrapReadError(1)
end subroutine




!x y z -> zlevel(mesh(1))%cm_zlev(mesh(2),mesh(3)))%finemesh(mesh(4),mesh(5),mesh(6))
!z_flag=0 default
!z_flag=z_level number 
!z_flag=-z_level number -> zfine=1
subroutine MapPos2Mesh(x,y,z,mesh,z_flag)
use paraset4
real x,y,z
integer mesh(6),z_flag

integer i, j, k !, ini,inj,ink

mesh=0
!z
if(z_flag .eq. 0 ) then
 if(z .eq. z_lev_pos(num_zlev+1))  then
  mesh(1)=num_zlev
 else if (z .gt. z_lev_pos(num_zlev+1)) then
  mesh(1)=num_zlev+1
 else
 do k=1, num_zlev
  if(z .ge. z_lev_pos(k) .and. z .lt. z_lev_pos(k+1) ) then 
   mesh(1)=k
   exit
  endif

 enddo
 endif

elseif( z_flag .gt. 0) then
 mesh(1)=z_flag
else
 mesh(1)=-z_flag
 mesh(6)=1
endif

!y
if(y .eq. zlevel(1)%y_cm_pos(zlevel(1)%ncy+1))  then
  mesh(3)=zlevel(1)%ncy
 else if (y .gt. zlevel(1)%y_cm_pos(zlevel(1)%ncy+1)) then
  mesh(3)=zlevel(1)%ncy+1
 else
 do j=1, zlevel(1)%ncy
  if(y .ge. zlevel(1)%y_cm_pos(j) .and. y .lt. zlevel(1)%y_cm_pos(j+1) ) then 
   mesh(3)=j
   exit
  endif
 enddo

endif

!x
if(x .eq. zlevel(1)%x_cm_pos(zlevel(1)%ncx+1))  then
  mesh(2)=zlevel(1)%ncx
else if(x .gt. zlevel(1)%x_cm_pos(zlevel(1)%ncx+1))  then
  mesh(2)=zlevel(1)%ncx+1
 else
 do i=1, zlevel(1)%ncx
  if(x .ge. zlevel(1)%x_cm_pos(i) .and. x .lt. zlevel(1)%x_cm_pos(i+1) ) then 
   mesh(2)=i
   exit
  endif
 enddo

endif

if( mesh(1) .gt. 0 .and. mesh(2) .gt. 0 .and. mesh(3) .gt. 0 .and. &
 mesh(1) .le. num_zlev .and. mesh(2) .le. zlevel(1)%ncx .and. &
 mesh(3) .le. zlevel(1)%ncy    ) then
  
mesh(4)=int((x-zlevel(mesh(1))%x_cm_pos(mesh(2)))/zlevel(mesh(1))%&
   cm_zlev(mesh(2),mesh(3))%finsize(1)) + 1
mesh(5)=int((y-zlevel(mesh(1))%y_cm_pos(mesh(3)))/zlevel(mesh(1))%&
    cm_zlev(mesh(2),mesh(3))%finsize(2)) + 1
if(z_flag .lt. 0) then
 mesh(6)=1
else
 mesh(6)=int((z-z_lev_pos(mesh(1)))/zlevel(mesh(1))%&
            cm_zlev(mesh(2),mesh(3))%finsize(3)) + 1
endif

endif
return
end subroutine

!average flux moments for each material
!for xs collapsing
Subroutine AvgMatFlux
use paraset1
use paraset2
use paraset3
use paraset4
use files
use funs
use ErrControl

integer i,j,k, g, m
integer cmi,cmj,cmk, fmi,fmj,fmk
integer num_fmesh(3)

real,allocatable :: flux_gm(:,:)

allocate( flux_gm(num_group, num_material))
flux_gm=0


do cmk=1, num_cmesh(3)
do cmj=1, num_cmesh(2)
do cmi=1, num_cmesh(1)

num_fmesh=zlevel(cmk)%cm_zlev(cmi,cmj)%fine
vol=zlevel(cmk)%cm_zlev(cmi,cmj)%fine_vol

do fmk=1, num_fmesh(3)
do fmj=1, num_fmesh(2)
do fmi=1, num_fmesh(1)

m=zlevel(cmk)%cm_zlev(cmi,cmj)%mat_matrix(fmi,fmj,fmk)

do g=1, num_group
 flux_gm(g,m) = flux_gm(g,m)+ vol* &
   zlevel(cmk)%cm_zlev(cmi,cmj)%fluxset(1)%locg_set(g)%locg_flux(fmi,fmj,fmk)
enddo !g

enddo  !fmi
enddo  !fmj
enddo  !fmk

enddo !cmi
enddo !cmj
enddo !cmk

do m=1, num_material
if(mat_vol(m) .ne. 0) then
flux_gm(:,m)=flux_gm(:,m)/mat_vol(m)
endif
enddo
if( IsFluxGM .eq. 1) then
cur_filename='flux.fgm'
else
cur_filename='flux.agm'
endif 
open(unit=53, file=trim(cur_filename), err=1001 )
do m=1, num_material
if(IsFluxGM .eq. 1) then
write(53,"('/Avg. Flux Spectrum for mat :', I0)", err=1002) m
write(53,"(500(ES12.5,X) )", err=1002) flux_gm(:,m)
else
write(53,"('/Avg. Adj. Flux Spectrum (Grp 1, 2, ...) for mat : ', I0)", err=1002) m
write(53,"(500(ES12.5,X))", err=1002) (flux_gm(g,m) , g=num_group,1, -1)
endif

enddo

close(53)

if(IsFluxGM .eq. 1) then
  write(*,"(' average flux per material zone... file:', A)" ) trim(cur_filename)
  write(READLOG,"(' average flux per material zone... file:', A)" ) trim(cur_filename)
else
  write(*,"('average adjoint flux per material zone...(GROUP ORDER FLIPPED), file: ',A)")  trim(cur_filename)
  write(READLOG,"('average adjoint flux per material zone...(GROUP ORDER FLIPPED), file: ',A)")  trim(cur_filename)
endif

return
1001  write(err_message,*) 'open file error: average flux per material zone'
      call TrapOpenError(1)
1002  write(err_message,*) 'error while writing average flux per material zone'
      call TrapOpenError(1)
end subroutine
