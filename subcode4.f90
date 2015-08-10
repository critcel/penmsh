!subroutine Output tecplot binary datafile
subroutine MxOutputBinGemetry
use paraset1
use paraset2
use paraset3
use paraset4
use files
use mMyTecplt


integer dir_len, ext_len, name_len,fullname_len
integer i,j,k
integer ini,inj,ink
integer ncm, ncm_len
integer, dimension(:), allocatable :: effgrp_num

integer grp
integer m,i_len
!IsHrtDat: generate heart only tecplot ascii file 
!IsNomalizeFlux: numalize flux to the max_flux
!ISCompareMCNP: compare pentran flux with MCNP

! real tmp_flux(2)

real x0,y0,z0
real,dimension(:),allocatable :: point
character(len=40) :: zone_name=''
!character(len=480) :: var_addon=''
character(len=20)  var_first5(5)
!integer:: addvar_len_sum=0

!character(len=530) var_tot

!integer num_addvar,len_first5
!integer,dimension(:),allocatable :: addvar_len
!character(len=20),dimension(:), allocatable :: addvar

data var_first5/'X','Y','Z','MatNum', 'SrcStr'/

!if(IsNomalizeFlux*IsProcessFlux .eq. 1 .or. IsNomalizeFlux*IsProcessCurr .eq. 1) then 
!  write(*,*)  'flux to be normalized to 1 (ie. the maxim flux in all groups equals to 1)' 
!  write(READLOG,*) 'flux to be normalized to 1 (ie. the maxim flux in all groups equals to 1)'
 
! do i=1,1+IsRefFlux
!   call NomalizeFlux(i)
! enddo

! endif

!point=0
IsHrtDat=0
IsHrtDat=IsHrtDat*IsProcessHrt
prbname_len=len_trim(prbname)
write(form, "('(A',I0,',A4)') ") prbname_len
write(outputfile(3)%name,form) prbname, '_mix' 
outputfile(3)%ext='plt'
outputfile(3)%dir='input_dir'
name_len=len_trim(outputfile(3)%name)
dir_len=len_trim(outputfile(3)%dir)
ext_len=len_trim(outputfile(3)%ext)

fullname_len=name_len+1+ext_len
do i=1,fullname_len
  if(i .le. name_len) then
    outputfile(3)%fullname(i:i)=outputfile(3)%name(i:i)
   else if(i .eq. name_len+1) then
    outputfile(3)%fullname(i:i)='.'
   else if(i .ge. name_len+2 .and. i .le. fullname_len) then
    outputfile(3)%fullname(i:i)=outputfile(3)%ext(i-name_len-1:i-name_len-1)
  endif
enddo
allocate(pltfile(1))
 pltfile(1)%unit=301
 pltfile(1)%tempunit=pltfile(1)%unit+111
 pltfile(1)%dir=input_dir
 pltfile(1)%filename=outputfile(3)%fullname
 pltfile(1)%byte_order=1
 pltfile(1)%title=prbname

if(num_group_eff .ne. 0) allocate(effgrp_num(num_group_eff))

if(IsProcessFlux .eq. 1 .or. IsProcessCurr .eq. 1) then
j=1
do i=1, num_group
if(grp_flag(i,1) .eq. 1) then
effgrp_num(j)=i
j=j+1
endif
enddo
endif

!for curr
pltfile(1)%NumVar=5
 if(IsProcessFlux .eq. 1) then 
   pltfile(1)%NumVar=5+num_group_eff
 endif
 if (IsProcessCurr .eq. 1) then
   pltfile(1)%NumVar=5+num_group_eff*5
 endif

 allocate(pltfile(1)%VarName(pltfile(1)%NumVar))
 do i=1, 5
  pltfile(1)%VarName(i)=var_first5(i)
 enddo
 if(IsProcessFlux .eq. 1) then

 k=5
 do i=1, num_group
  if(grp_flag(i,1) .eq. 0) cycle
  k=k+1
  call intlen(i, i_len)
  write(form,"('(A4,I',I0,')')") i_len 
  write(pltfile(1)%VarName(k),form) ' Grp',i
 enddo

 endif

!for curr
if(IsProcessCurr .eq. 1) then
 k=5
 do i=1, num_group
  if(grp_flag(i,1) .eq. 0) cycle
  k=k+1
  call intlen(i, i_len)
  write(form,"('(A4,I',I0,')')") i_len 
  write(pltfile(1)%VarName(k),form) ' Grp',i
 enddo

 i=5+num_group_eff
 do j=1, num_group
  if(grp_flag(j,1) .eq. 0) cycle
  i=i+1
  call intlen(j, i_len)
  write(form,"('(A4,I',I0,',A1)')") i_len
  write(pltfile(1)%VarName(i),form) ' Grp',j,'J'
  
  i=i+1
  call intlen(j, i_len)
  write(form,"('(A4,I',I0,',A2)')") i_len
  write(pltfile(1)%VarName(i),form) ' Grp',j,'Jx'
  
  i=i+1
  call intlen(j, i_len)
  write(form,"('(A4,I',I0,',A2)')") i_len
  write(pltfile(1)%VarName(i),form) ' Grp',j,'Jy'
  
  i=i+1
  call intlen(j, i_len)
  write(form,"('(A4,I',I0,',A2)')") i_len
  write(pltfile(1)%VarName(i),form) ' Grp',j,'Jz'
 enddo

endif
 

 allocate(pltfile(1)%zne_datafmt(pltfile(1)%NumVar))
!1=float, 
 pltfile(1)%zne_datafmt=1
 
! tecflag=myTecini(1)
! tecflag=myTeczne(1,'xyz',1,2,1,1)
! point(1)=1
! point(2)=2
! point(3)=3
! do i=1,2
!   tecflag=myTecDat(1,3,point)
! enddo
! tecflag=myTecEnd(1)
! deallocate(pltfile)

IsTec=myTecIni(1)
allocate(point(pltfile(1)%NumVar))
point=0


if(IsProcessFlux .eq. 1) then
    
  if(IsRefFlux .eq. 1) then
   write(*,*)  'flux diff. will be written into the tecplot data file' 
   write(READLOG,*)  'flux diff.  written into the tecplot data file'
  else
   write(*,*)  'flux will be written into the tecplot data file'
   write(READLOG,*)  'flux written into the tecplot data file'
  endif

endif

!for curr

if(IsProcessCurr .eq. 1) then
   write(*,*)  'flux&current will be written into the tecplot data file'
   write(READLOG,*)  'flux&current written into the tecplot data file'
endif

do k=1, num_zlev
 do j=1,zlevel(k)%ncy
  do i=1, zlevel(k)%ncx
      ncm=(k-1)*num_cmesh(1)*num_cmesh(2)+(j-1)*num_cmesh(1)+i
   call intlen(ncm,ncm_len)
     write(form,"('(A3,I',I0,')' )")  ncm_len 
     write(zone_name,form) 'CM_', ncm
!... Write the zone header information.
      IsTec=cmTecZne(1,zone_name,1,zlevel(k)%cm_zlev(i,j)%fine(1),&
     zlevel(k)%cm_zlev(i,j)%fine(2), zlevel(k)%cm_zlev(i,j)%fine(3))
!      IsTec=cmTecZne(1,zone_name,1,zlevel(k)%cm_zlev(i,j)%fine(1)*&
!     zlevel(k)%cm_zlev(i,j)%fine(2), zlevel(k)%cm_zlev(i,j)%fine(3),1)

    enddo
  enddo
enddo 
IsTec=cmHeaderEnd(1)

do k=1, num_zlev 
 z0=z_lev_pos(k)
 do j=1,zlevel(k)%ncy
   y0=zlevel(k)%y_cm_pos(j)  
   do i=1, zlevel(k)%ncx
      x0=zlevel(k)%x_cm_pos(i)
      IsTec=cmZoneHeader(1)
      do ink=1, zlevel(k)%cm_zlev(i,j)%fine(3)
        point(3)=z0+(ink-1)*zlevel(k)%cm_zlev(i,j)%finsize(3)+0.5*zlevel(k)%cm_zlev(i,j)%finsize(3)
        do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
        point(2)=y0+(inj-1)*zlevel(k)%cm_zlev(i,j)%finsize(2)+0.5*zlevel(k)%cm_zlev(i,j)%finsize(2)
        do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1)
          point(1)=x0+(ini-1)*zlevel(k)%cm_zlev(i,j)%finsize(1)+0.5*zlevel(k)%cm_zlev(i,j)%finsize(1)     
             
          point(4)=zlevel(k)%cm_zlev(i,j)%mat_matrix(ini,inj,ink)
          if(zlevel(k)%cm_zlev(i,j)%s_flag .ge. 1 ) then
            point(5)=zlevel(k)%cm_zlev(i,j)%src_matrix(ini,inj,ink)
          else
            point(5)=0
          endif
!processing flux
     if(IsProcessFlux .eq. 1) then
       m=5
       do grp=1,num_group
       if(grp_flag(grp,1) .eq. 0) cycle
         m=m+1
         if(IsRefFlux .eq. 1) then
           point(m)=zlevel(k)%cm_zlev(i,j)%fluxset(2)%locg_set(grp)%locg_flux(ini,inj,ink)
         else
           point(m)=zlevel(k)%cm_zlev(i,j)%fluxset(1)%locg_set(grp)%locg_flux(ini,inj,ink)
         endif
       enddo
     endif

!processing curr
  if(IsProcessCurr .eq. 1) then
  m=5
  do grp=1, num_group
    if(grp_flag(grp,1) .eq. 0) cycle
    m=m+1
    point(m)=zlevel(k)%cm_zlev(i,j)%fluxset(1)%locg_set(grp)%locg_flux(ini,inj,ink)
  enddo
 m=6+num_group_eff
 do grp=1,num_group
     if(grp_flag(grp,1) .eq. 0) cycle
        point(m:m+3)=zlevel(k)%cm_zlev(i,j)%fluxset(1)%locg_set(grp)%locg_curr(ini,inj,ink)%fcurr(:)
 m=m+4
    enddo
   endif

!... Write out the field data.
!processing flux
     
     IsTec=mxTecDat(1,pltfile(1)%NumVar,point)   
 
     if(IsHrtDat .eq. 1) then
      if(point(4) .gt. 1) then
        if(IsProcessFlux .eq. 0) then
          write(form,"('(A3,I',I0,')' )")  ncm_len 
          write(11,'(5(ES13.6,2X))') (point(m),m=1,5)
        elseif(IsProcessFlux .eq. 1) then
          write(form,"( '(',I0,'(E12.6,2X)',')' )")  5+num_group 
          write(11,form) (point(m),m=1,5+num_group)
        endif
      endif
     endif

   enddo !ini
  enddo  !inj
enddo !ink
      
   enddo
  enddo
enddo

IsTec=cmTecEnd(1)

if(IsProcessFlux .eq. 1 .or. IsProcessCurr .eq. 1) then
 
 write(*,"('  Total ', I0, ' Grp flux info has been written')")  num_group_eff 
 write(READLOG,"(' Total ', I0, ' Grp flux info has been written')")  num_group_eff
 if(num_group_eff .ne. 0) then
  write(*,"('  Including Grp(s): ', 10(I0,1x))")  effgrp_num 
  write(READLOG,"('  Including Grp(s): ', 10(I0,1x))")  effgrp_num
 endif
  
endif


call OutputMcr
deallocate(pltfile)
return 
end subroutine

!output macro file for tecplot
!grouping CMs in one z level
subroutine OutputMcr
use paraset4
use paraset2, only: num_cmesh
use files
use funs
use ErrControl

curfile%fullname=''
curfile%name=prbname
curfile%ext='mcr'
cur_filename=Getfullname(curfile,0)

open(unit=421,file=cur_filename)
write(421,"('#!MC 900')")
write(421,"('$!VarSet |num_z|=',i3)" ) num_zlev
write(421,"('$!VarSet |cm_per_zlev|=',i4)" ) num_cmesh(1)*num_cmesh(2)
write(421,"('$!FIELDLAYERS SHOWMESH = NO')" )
write(421,"('$!FIELDLAYERS SHOWBOUNDARY = NO')" )
write(421,"('$!FRAMEMODE = THREED')")
write(421,"('$!GLOBALCONTOUR VAR = 4')")
write(421,*) '$!VarSet |cm_tot|=(|cm_per_zlev|*|num_z|)'
write(421,*) '$!FIELD [1-|cm_tot|] SCATTER{COLOR = MULTI}'
write(421,*)  '$!FIELD [1-|cm_tot|] SCATTER{ISFILLED = YES}'
write(421,*)  '$!FIELD [1-|cm_tot|] SCATTER{FILLCOLOR = MULTI}'
if (tot_num_fm .gt. 500000) then
 write(421,*)  '$!FIELD [1-|cm_tot|] SCATTER{FRAMESIZE = 0.5}'
else
 write(421,*)  '$!FIELD [1-|cm_tot|] SCATTER{FRAMESIZE = 1}'
endif
write(421,*)  '$!LOOP |num_z|'
write(421,*)  '$!VARSET |a| = ((|loop|-1)*|cm_per_zlev|+1)'
write(421,*)  '$!VARSET |b| = (|loop|*|cm_per_zlev|)'
write(421,*)  '$!FIELD [|a|-|b|]  GROUP = |loop|'
write(421,*)  '$!ENDLOOP'
write(421,"('$!FIELDLAYERS SHOWSCATTER = YES')" )
close(421)

write(warning_message,"('Macro file for Tecplot, file: ' ,A)") trim(cur_filename)
call DisplayMsg
end subroutine

