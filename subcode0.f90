!parameters set1 - General problem defination paras 
!Reserved for future use
!file subcode0.inp: modules and functions
module constants
real*8:: pi=3.14159265358979
real :: tiny=1e-15
end module

Module paraset1

!timing
real ::  stime
integer :: IsProcessFlux,IsSkipGeometry,IsRefFlux
integer :: IsProcessCurr=0 , IsPlotZlev=1, MidPlot=1, MatSrcFlx=1
integer :: Is3DPlot=0, IsDrawFm=1
real    :: PngSizeFactor=1.0
integer :: ColorMap=3
integer :: IsOffPentran=0, IsOffTitan=1
integer :: IsFidoSrc=1
integer :: IsFluxout=0
integer :: IsFov=0  !detector fov for array adjoint model

character(LEN=30) :: CaseName='' 

integer :: IsZPlot=0,  num_zloc=0
integer :: IsYPlot=0,  num_yloc=0
integer :: IsXPlot=0,  num_xloc=0
integer :: IsLogPlot=0
!3d plot view angle and distance 
real :: view_ang=0, view_pos(3)=4 
integer :: num_vpos=0   
integer, parameter :: maxplot=50
real :: xyzloc(maxplot,3)=0

!If fluxfile name =prb1.flx for each fluxset 
integer :: IsPrbDotFlx(3)=0
integer :: IsNomalizeFlux=0
!output penmsh .out files 
integer :: IsDotOut=0
!output material number hdf5 file
integer :: IsDotH5=0
!ooutput avg. flux for each material
integer :: IsFluxGM=0, IsChiPath=0

!inputfile directory
character*80 input_dir,flux_dir, ref_flux_dir, curr_dir, ChiPath_dir
character*80,dimension(:),allocatable :: fluxset_dir
integer, dimension(:,:), allocatable :: grp_flag
integer :: num_group_eff=0
real :: n_factor=0.0  !global normalization factor
real, allocatable :: n_factor_grp(:)  !group normalization factor
!flux plot range
real :: flx_min=0, flx_max=0
integer :: flx_min_flag=0, flx_max_flag=0

! for heart model
integer :: IsHrtDat=0
integer :: IsProcessHrt=0

character*80  :: hrt_filename=''
character*50 ::  ver_num=""
character(len=256) :: cmd_line=''

!for mox problem
 !IsPosition=0 : upper right
 !IsPosition=1 : lower right
 integer :: IsPosition=1
 !IsDiffPower=1 : calculate power diff. 
 integer :: IsDiffPower=0
 integer :: IsCalPinPower=0
 integer :: IsRefpenmsh=0
!rodded a,b and unrodded
 integer :: IsCase

! for compare flux with MCNP
integer :: IsCompareMCNP=0

!for venus problem
integer :: IsVenusDose=0

type tGlobalfluxset

 real :: max_glb_flux=0, min_glb_flux=0
 real :: avg_glb_flux=0
 real,dimension(:),allocatable :: max_glb_flux_grp, min_glb_flux_grp, avg_glb_flux_grp
 
 character(len=80) :: fluxfile_dir=''
 character(len=80) :: fluxfile_name=''

end type tGlobalfluxset

type(tGlobalfluxset),dimension(:),allocatable :: glb_fluxset

real :: max_src=0
!control spacef value format
!DotBit=3, i.e 45.123, 333.123
!DotBit=-3 i.e  1.123E-01
integer :: DotBit=-5

!pentran varible ndmeth
integer :: ndmeth_global=2
integer :: sec_flag(0:15)=0

!for individual z-level plot
integer :: z_start=0, z_end=0, IsZSingle=0

!for splitting the flux file
integer :: num_flx_out=1

!set nmesh to give a certain mesh size
real :: max_fm_size=-1.0

end module

Module paraset2  ! Gemetry and Source

type cube
character*1 var1 ! x, y,z
character*6 var2 !xmesh=, ymesh=, zmesh=
character*7 var3 !ixfine=
integer     int3 !=6
character*6 var4 !ixmed=,jymed=,kzmed=
integer     int4 !=5
!@cgdim(1,2,3)%boundary: cm position along x, y,z
real, dimension(:), allocatable:: boundary
!pixb: # of pixels along x,y,and z, for heart model 
real, dimension(:), allocatable:: pixb
end type cube

type(cube) :: cgdim(3)

!@num_cm: number of coarse mesh
!@num_cm_zlev: number of coarse mesh every z level
integer :: num_cm,num_cm_zlev
!@num_cmesh(1,2,3): number of x,y,z coarse meshes
integer :: num_cmesh(3)

! number of sources
integer :: num_source
real,dimension(:),allocatable :: x_pos_src, y_pos_src,z_pos_src
real,dimension(:,:,:), allocatable :: src_mesh, src_mesh_vol

!mat_vol: vol for each mat
!mat_fm : num of fm for each mat
!mat_list: num of mat in one cm
real*8, dimension(:), allocatable:: mat_vol
integer, dimension(:), allocatable::mat_fm
integer,dimension(:),allocatable :: mat_list
real*8 :: total_vol=0, total_mat_vol=0,tgt_tot_vol=0,tgt_tot_mass=0,mod_tot_mass=0

!for matbal.dis
type matp

character*24 :: nam='noname'
integer :: num=0
real*8 ::  vol_mod=0
real*8 ::  vol_tgt=0
real ::   den=0

end type matp


type(matp),dimension(:), allocatable :: mat_tgt

!center of material
real, allocatable :: cen_mat(:,:)
end module

Module paraset3 !for pentran 

! PENTRAN VARIABLE NAMES BY BLOCK:
!   BLOCK 1 - 12 Possible:
! ngeom ngroup isn nmatl ixcrs jycrs kzcrs decmpv lodbal timcut tolmgd modadj

!   BLOCK 2 - 12 Possible:
! xmesh ixmed ixfine ymesh jymed jyfine zmesh kzmed kzfine nmattp flxini mathmg

!   BLOCK 3 - 10 Possible:
! lib legord nxtyp ihm iht ihs ihng chig nxcmnt legoxs

!   BLOCK 4 - 12 Possible:
! nprtyp nrdblk tolin tolout maxitr methit
!                                     methac ncoupl ndmeth nzonrb dtwmxw nquit

!   BLOCK 5 - 11 Possible:
!nsdef nscmsh ssnrm sref serg smag rkdef spacpf omegap scalsf kextrp

!   BLOCK 6 -  6 Possible:
!ibback ibfrnt jbeast jbwest kbsout kbnort

!   BLOCK 7 -  9 Possible:
!nxspr ngeopr nsumpr meshpr nfdump nsrcpr nsdump nmatpr nadump


integer :: sn,group
!@@num_mat: number of  materials
integer :: num_mat

!header
!@maxcrs: Max Coarse Meshes/1 axis 
!@maxmed: Max Medium Meshes/1 axis
!@maxfin: Max Fine   Meshes/1 axis  
!@maxlin: maxinum lines in Input dest
!@maxarr Maximum number inputs in fido vector
!@nctlim: Maximum Fido chars per varible
integer :: maxmem
integer :: maxcrs,maxmed,maxfin
integer*8 :: maxlin, maxarr,nctlim

!@maxmmc: Max Medium Meshes/Coarse 
!@maxfmc: Max Fine   Meshes/Coarse  
integer :: maxmmc, maxfmc
!@: Max Leg Scat dimensioned  
integer :: maxleg

!block1
integer :: ngroup(3)
!block2
integer, dimension(:), allocatable :: mathmg
real, dimension(:), allocatable :: flxini

!block 3
!@legord: Legendre scattering order
!@legoxs: legendre scattering order
!@nxtyp: cross section type
!@ihm  number of rows of xsec data
!@iht  total xsec position
!@ihs  whith group scatter xsec position
!@ihng positon of last nertron scatter xsec prior to coupled gamma xsecs
integer :: legord,legoxs
integer :: nxtyp,ihm,iht,ihs,ihng
!@chig  group fission probavilities for each mat and group
!chig(num_mat,group)
real, dimension(:), allocatable :: chig
!@ number of crsss section comment cards
integer :: nxcmnt

!block4
integer :: nprtyp, nrdblk
!real, dimension(:), allocatable :: tolin
!integer :: tolin_len
! tolout tolout(1)=tolerance, tolout(2)=med_grid_multiplier
real :: tolout(2), tolin(2)

! maxitr(1): max_no_of_iter, maxitr(2):criticality_inner_limit
integer :: maxitr(2)
integer :: methit, methac,ncoupl
!nzonrb
integer :: number_of_zones, skip_iter
real  :: damping_fact
real :: dtwmxw
integer :: nquit

! block 5
integer, dimension(:), allocatable :: nscmsh,nsdef

real, dimension(:),allocatable :: smag,serg,sref
real :: rkdef=1.0

! block 6 BCs

character*7 bc_char(6)
! bc_char= ibback, ibfrnt,jbeast,jbwest,kbsout,kbnort
! assignment in subroutine pentran_var

type albedo_typ
 real, dimension(:),allocatable :: albedo_grp 
end type albedo_typ

type(albedo_typ) albedo(6)

!block 7
integer, dimension(:), allocatable:: nadump
!for dump angular flux
integer, dimension(:), allocatable:: front, back, left, right


end module

!for TITAN 
module paraset5  
!TITAN
integer :: num_quad=1
integer :: IsSpect=0
integer :: cmdiff_global=0

!spect
! starting angle , ending angle, number of projection
real :: ang_start=0, ang_end=360
integer :: num_ang=4

!rotation
real :: vec_axis(3)=(/0,0,1/)
real :: pos_axis(3)
real :: r_axis

!collimator effects with circular splitting
integer :: spt_order=0 , spt_circle=1
real :: spt_angle

!detector
real size_det(2)
integer :: pix_det(2)=1


type tSplit

integer :: spid, order,spt_id,topnum
real :: alpha=0.02   !circular splitting derailed angle in radian
end type

type tQuad
integer :: tquad, oquad, num_split
type(tSplit), dimension(:), allocatable :: direc_split

end type

type(tQuad), dimension(:), allocatable :: quad
integer :: num_split

end module

!for penmsh 
Module paraset4
!penmsh.inp 
character*80 prbname,prbname_ref
integer :: prbname_len=1, prbname_len_ref=1

integer :: num_zlev,num_material,flag_mathematica
integer :: num_zlev_ref,num_material_ref,flag_mathematica_ref

real, dimension(:), allocatable::  z_lev_pos
real, dimension(:), allocatable::  z_lev_pos_ref

integer,dimension(:), allocatable:: max_zfine
integer,dimension(:), allocatable:: max_zfine_ref

integer,dimension(:), allocatable:: ratio_x, ratio_y,ratio_z
integer :: s_format=2, xnum_s_mesh=1, ynum_s_mesh=1, znum_s_mesh=1
integer :: num_group=1, sn_order=8,pn_order=0
integer :: xs_format,num_comment,leg_order,len_table
integer ::  bcs(6)
! num_src_mat: number of source mat
! id_s_mat: mat id of the source material id_s_mat(num_src_mat)
! s_intensity(num_src_mat), s_mat_src(num_material)
integer,dimension(:), allocatable:: id_s_mat
real,dimension(:),allocatable::  s_intensity,s_mag_src
integer :: num_src_mat=0

type d3array
integer,dimension(:,:,:), allocatable :: d3_int
end type

!overlay block
type tBlock
 integer overlay_num
 integer,dimension(:),allocatable:: overlay_typ,overlay_mat
 real,dimension(:,:), allocatable:: overlay_bon 
 integer,dimension(:,:), allocatable :: num_rep
 real,dimension(:,:), allocatable :: dist_rep

!repeated block with diff. mat
 type(d3array), dimension(:), allocatable :: mat_rep
end type tBlock

!for cm type 3
!type tCMtype
! real,dimension(:),allocatable :: r_lft, r_rgt, theta_bot, theta_top
! integer, dimension(:),allocatable :: cm_mat
!end type tCMtype

! current vector: magnitude, Jx, Jy, Jz
type tCurr
 real :: fcurr(4) 

end type tCurr

type tFsetg
  real, dimension(:,:,:), allocatable :: locg_flux
! for curr
  type(tCurr), dimension(:,:,:), allocatable :: locg_curr
end type

! flux set 
type tFluxset
  real ::	max_loc_flux=0, min_loc_flux
  real ::  avg_loc_flux
  real, dimension(:), allocatable :: max_loc_flux_grp, min_loc_flux_grp, avg_loc_flux_grp
  type(tFsetg), dimension(:), allocatable :: locg_set
  
end type tFluxset


! prbname1.inp
type zlevcm
 integer :: cm_type=0, num_subregion=0
 integer :: fine(3)=0
 integer cm_mat_num
 integer :: tot_fm=0,uni_zfine=0

 type(tBlock), dimension(:),allocatable :: overlay_block
! type(tCMtype) dimension(:),allocatable :: cm_type_data
 real*8 finsize(3), cm_vol, fine_vol 
!for cm type 3
 real,dimension(:),allocatable :: r_lft, r_rgt, theta_bot, theta_top
 integer, dimension(:),allocatable :: cm_mat

 integer med(3)
 !integer index_cm(3) ! coarse mesh 3D index.
 integer :: ndmeth=1, solver_id=0, qudra_id=1
 real    sum_act
 integer s_flag
 integer cm_num_mat
 integer :: num_block=1

 real,dimension(:), allocatable :: cm_x,cm_y,cm_z
!assigned number for each fine mesh in one coarse mesh
!data_matrix(x,y index, z level index, data source idenifier 
 integer, dimension(:,:,:), allocatable :: mat_matrix
 real, dimension(:,:,:), allocatable :: src_matrix

!num of fm for each mat
 integer, dimension(:), allocatable :: fm_num_mat
! for pentran flux output 
 type(tFluxset), dimension(:), allocatable :: fluxset

end type zlevcm

type zlev
 
 integer ncx,ncy,num_zfine
 real, dimension(:), allocatable:: x_cm_pos,y_cm_pos 
 type(zlevcm),dimension(:,:),allocatable::  cm_zlev
  
end type zlev

type(zlev),dimension(:), allocatable::  zlevel

integer :: max_fm
!maxium num of overlay specifications
integer:: maxos=10 
!overlay shape parameters
!shape id 1

!shape id 2
!put other shape id here
real center_x, center_y, r
!shape id 3
real Ax,Ay, Bx,By, Cx,Cy, Zx, Zy
!Barycentric Coordinates
real b0,b1,b2,b3
!shape id 4 sphere
real center_z
!shape id 5 cone
real theta,h
!shape id 6 hexagon
real vertex_x, vertex_y
real gamma, beta, r_in, r_out, cos_theta, sin_theta, dk
real, parameter :: cos30=0.866025

!reference varibles
integer*8 :: tot_num_fm, maxmem_byte,num_dir

end module

Module files

type fileinfo
 character(len=80) :: fullname=''
 character(Len=50) :: name=''
 character(len=50) :: info=''
 character(len=80) :: dir=''
 character(len=10) :: ext=''
end type fileinfo

!file unit number
integer :: READLOG=499
integer :: TEMPF90=599
integer :: OUTF90=600
integer :: FLUXUNIT=213

type(fileinfo), dimension(:), allocatable :: inputfile, fluxfile,matfile

character(len=100) form
!outputfile(1): log file; outputfile(2): tecplot material and source dist.
!outputfile(3): tecplot binary data file
!outputfile(4): prbname_out.f90 
!outputfile(5): temp.f90 (deleted after done)
!outputfile(6): prbname_titan.out 
type(fileinfo),save ::  outputfile(6)
type(fileinfo),save ::  penmsh_inp, penmsh_inp_ref,prb_xs, prb_inp
type(fileinfo),save ::  curfile
type(fileinfo), save :: nfile
end module


! for generating fido array
Module fido
!@oneline1: output array
!e.g. chig=2R0.3 ...     prevar='chig=', oneline='2R0.3...'
 character, dimension(:), allocatable :: oneline1
 real,dimension(:), allocatable :: fido_input
 integer :: output_len,fido_input_len
 character*40 prevar, appvar
  character*40 fidoform

 integer :: max_input_len, max_fido_char
 integer :: num_fido_char
 integer :: buf_size=10000, buf_lim=50
 
!for spacpf
 integer :: input_len_sp=0, fido_char_sp=0

end module

! error and warning traps
Module ErrControl

character*120 :: cur_filename=''
character*120 :: cur_filedir=''
character*120 err_message,warning_message
character*120 :: err_message_2="" , warning_message_2=""
character*80  cur_var


integer:: cur_fileunit=0
integer:: warning_times(10)=0
integer:: max_warning=5
end module


Module lineread

!maxium length for one line
parameter(MAXLEN=120)
parameter(NUMFLAG=4)
integer :: flag(NUMFLAG)=0
contains
! preread a line, if it's comment line, continue to read next line 
! until a effective line is reached
! return num of comment lines read, and the file position is on the begining
! of the first effective line
!flag: flag(1):
!       on output: number of values in a line
!       on input: =-1 turn on output  
!flag(2)=0
integer function NumCmtLine(unitnum,flag)
use ErrControl
integer unitnum, flag(NUMFLAG)

character(len=1) firstchar
character(len=3) markchar
character(len=MAXLEN+1) line

integer cmt_err,i,j
integer NumVal, IsAnotherNum
integer debug

markchar='/#!'
firstchar='/'
NumCmtLine=-1

! do while(firstchar .eq. markchar .and. i .le. MAXLEN)
!do while(firstchar .eq. markchar)
do while(index(markchar, firstchar) .ne. 0 )

 read(unitnum,'(A)',iostat=cmt_err) line

 if(cmt_err .ne. 0) then
  if(cmt_err .gt. 0) then
   write(err_message, "('read error when reading ', A )") trim(cur_var)
   call TrapReadError(1)
  else
! read to the end of the file
   write(err_message, "('EOF error when reading ', A )") trim(cur_var) 
! flag(2)=1: tolerate EOF   
   if(flag(2) .eq. 0)  call TrapReadError(1)
  endif
  
  NumCmtLine=-1
  return
 else
  write(err_message, "('read error when reading ', A )") trim(cur_var)
 endif

 !skip the space
 i=1
 do while(line(i:i) .eq. ' ' .and. i .le. MAXLEN)
  i=i+1
 enddo
 ! skip a blank line
 if(i .eq. (MAXLEN+1)) then 
   firstchar='/'
   i=1
 else
   firstchar=line(i:i)
 endif
 NumCmtLine=NumCmtLine+1
enddo

!Add your expression recognizing code here
!count number of values in one line
NumVal=flag(1)
if (NumVal .eq. -1) then
   IsAnotherNum=0
   NumVal=0
   do j=1, MAXLEN
 !  if(line(j:j) .eq. '/') exit
   if( index(markchar, line(j:j)) .ne. 0) exit
   if (IsAnotherNum .eq. 0) then
     !check if line(j:j) = 0-9 or .
       debug=index('-0123456789.',line(j:j))
	   if (debug .ne. 0) then ! a number encountered
	    IsAnotherNum=1
		NumVal=NumVal+1 
	   endif
	 elseif (IsAnotherNum .eq. 1) then
     !check if line(j:j)=' ' or ','
       if (line(j:j) .eq. ' ' .or. line(j:j) .eq. ',') then
         IsAnotherNum=0
	   endif
	 endif
   
   enddo

endif
flag(1)=NumVal
backspace(unitnum,err=2001)

goto 2000

2001 write(err_message,*) 'ERR3000: backspace statement error in file: '
call TrapReadError(1)

2000 return
end function NumCmtLine

end module


Module funs

!compare two int array with length len
!if same, return 1, otherwise return 0
CONTAINS
integer function int_compare(array1, array2, len)
integer len
integer array1(len), array2(len)

integer i,value

value=1
i=1
do while(value .ne. 0 .and. i .le. len)
  if(array1(i) .ne. array2(i))  value=0 
  i=i+1
enddo
int_compare=value
return
end function

!same as int_compare, except for real arrays
integer function real_compare(array1, array2, len)
integer len
real array1(len), array2(len)

integer i,value

value=1
i=1
do while(value .ne. 0 .and. i .le. len)
  if(array1(i) .ne. array2(i))  value=0 
  i=i+1
enddo
real_compare=value
return
end function

!used for mox problem, zlev map to lay
integer function map_lev2lay(k,num_of_zlays,num_zlev_lay)
integer k,num_of_zlays
integer num_zlev_lay(num_of_zlays)

integer i, flag,sumzlev(num_of_zlays)

flag=0
sumzlev(1)=num_zlev_lay(1)

do i=2,num_of_zlays
  
 sumzlev(i)=sumzlev(i-1)+num_zlev_lay(i)
 
enddo
i=1
do while(flag .eq. 0 .and. i .le. num_of_zlays)
  
  if(k .le. sumzlev(i) )  flag=1
  
  i=i+1
enddo
map_lev2lay=i-1
return
end function

!map overlay type to the number of boundary specification needed
integer function map_overlay_typ2bon(overlay_typ,cm_num_1z)
use ErrControl
integer overlay_typ,cm_num_1z

select case (overlay_typ)
  !sector
  case(51)
   !center_x, center_y, point1_x,point1_y, point2_x,point2_y, (bottom, top)
   map_overlay_typ2bon=6
  !quarter circle
  case(21,22,23,24,25,26)
   !xleft, xright, ybottom, ytop (zbottom, ztop)
   map_overlay_typ2bon=4
   !right angle triangle
  case(31,32,33,34)
   !xleft, xright, ybottom, ytop (zbottom, ztop)
   map_overlay_typ2bon=4
  !retagular
  case(1)
   !xleft, xright, ybottom, ytop (zbottom, ztop)
   map_overlay_typ2bon=4
  !cylinder
  case(2) 
    map_overlay_typ2bon=4   
  !triangle
  case(3)
    !x1,y1; x2,y2; x3,y3 (bottom, top)
	map_overlay_typ2bon=6
  !sphere
  case(4)
    !center_x, center_y, center_z, ra, (bottom, top)
	map_overlay_typ2bon=4
  !cone
 ! case(5)
    !x0 y0,z0,h,ra (h<0 upside down)
!	map_overlay_typ2bon=5

  !add your overlay type here
  case (6)
   map_overlay_typ2bon=4
  case default
    map_overlay_tpy2bon=4
	write(err_message,"('overlay type ',I0, ' in CM ', I0,' not supported ' )")&
        overlay_typ,cm_num_1z
    call TrapInputError(1)		  
end select

return
end function
 
! get file name
character*120 function Getfullname(efile, pathflag)
use files
type(fileinfo) :: efile
integer pathflag

if(efile%fullname .eq. '') then
 if(pathflag .eq. 0) then
  write(form,"('(A',I0,',A1,A',I0,')' )") len_trim(efile%name),len_trim(efile%ext)
  write(Getfullname,form  ) efile%name, '.', efile%ext
 else
   write(form, "('(A',I0,',A1,A',I0,',A1,A',I0,')' )") len_trim(efile%dir),& 
        len_trim(efile%name), len_trim(efile%ext)
   write(Getfullname,form ) efile%dir,'/',efile%name,'.', efile%ext
 endif

else

 if(pathflag .eq. 0) then
  write(form,"('(A',I0,')')") len_trim(efile%fullname)
  write(Getfullname,form ) efile%fullname 
 else
  write(form,"('(A',I0,',A1,A',I0,')')")len_trim(efile%dir), len_trim(efile%fullname)
  write(Getfullname,form ) efile%dir,'/',efile%fullname
 endif
endif
return 
end function

integer function CharNum(var,num)
character(*) var
integer num

character*20 :: var_in='12'
integer :: i=1
integer :: len_var=0
integer :: sign=1

CharNum=1
num=0

var_in=adjustl(var)
len_var=len_trim(var_in)
do while(i .le. len_var .and. CharNum .ne. 0)
  if (i .eq. 1) then
    CharNum=scan('-0123456789',var_in(i:i))
    if(CharNum .ne. 0) Then
	 if (CharNum .eq. 1) then
	  sign=-1
	 else 
      num=num+(CharNum-2)*10**(len_var-i)
     endif
    endif
  else
    CharNum=scan('0123456789',var_in(i:i))
    if(CharNum .ne. 0) num=num+(CharNum-1)*10**(len_var-i)
  endif
  i=i+1
enddo
num=sign*num

end function 

integer function CharFloat(var,num)
character(*) var
real num

integer ierr
CharFloat=1
num=0

read(var, *, iostat=ierr) num
if (ierr .ne. 0) CharFloat=0

end function 

end module

!to process input deck and decode FIDO
module mIn4deck

integer, parameter :: MAXLENA=5000
integer, parameter :: MAXUN=200
character (LEN=MAXLENA) :: in4line=''
character (LEN=MAXLENA) :: lineseg=''

character (LEN=3) ::  cmtchar='/#!'
character (len=3) ::  cmmarker='cm='
!character (LEN=3) ::  secmark="#"
character (LEN=3) ::  omitchar=",="
character (LEN=54) :: letterchar='/!abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
character (LEN=16)   :: digchar='0123456789.eE+-'
character (LEN=3)  ::  fidochar='RQZ' 

real, allocatable :: fidobank(:)
integer :: num_need=0, num_read=0

end module
