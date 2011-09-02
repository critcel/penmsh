!subroutines to write tecplot binary files


!BINARY FILE FORMAT:
! -----------------------------------------------------------------------
!The binary Datafile format (as produced by the preplot) is
! described below.


!The binary datafile has two main sections.  A header section and a data
! section.

!         +----------------+
!         | HEADER SECTION |
!         +----------------+
!         +---------+
!         |FLOAT32  |               value=357.0 (See below)
!        +---------+
!         +----------------+
!         | DATA SECTION   |
!         +----------------+


! I.  The header section.

!    The header section contains The version number of the file, a title
!    of the file, the names of the varialles to be plotted, the
!    descriptions of all zones to be read in and all text and geometry
!    definitions.
!
!     i.  Magic number, Version number
!         +-----------+
!         | "#!TDV75 "|       8 Bytes, exact characters "#!TDV75 ".
!         +-----------+       Version number is listed after the V.
!
!     ii. Integer value of 1.
!         +-----------+
!         | INT32     |       This is used to determine the byte order
!         +-----------+       of the reader relative to the writer.

!     iii. Title and variable names.
!         +-----------+
!         | INT32*N   |       The TITLE. (See note 1.)
!         +-----------+

!         +-----------+
!         | INT32     |       Number of variables (NumVar) in the datafile.
!         +-----------+
!
!         +-----------+
!         | INT32*N   |       Variable names.  N =  L[1] + L[2] + .... L[NumVar]
!         +-----------+       where:
!                                    L[i] = length of the ith variable name + 1
!                                           (for the terminating 0 value).
!                             (See note 1.)
!     iv.  Zones
!         +-----------+
!         | FLOAT32   |       Zone marker. Value = 299.0
!         +-----------+
!         +-----------+
!         | INT32*N   |       Zone name. (See note 1.) N = length of zone name + 1.
!         +-----------+
!         +-----------+
!         | INT32     |       Format, 0=BLOCK 1=POINT 2=FEBLOCK 3=FEPOINT
!         +-----------+
!         +-----------+
!         | INT32     |       Zone Color (set to -1 if you want tecplot to
!         +-----------+       determine).
!         +-----------+
!         | INT32*3   |       IMax,JMax,KMax  (See Note 2)
!         +-----------+
!      v.  Geometries
!         +-----------+
!         | FLOAT32   |       Geometry marker.  Value = 399.0
!         +-----------+
!         +-----------+
!         | INT32     |       CoordSystem 0=Grid, 1=Frame, 2=FrameOffset(not used)
!         +-----------+
!         +-----------+
!        | INT32     |       Scope 0=Global 1=Local
!         +-----------+
!         +-----------+
!         | FLOAT64*3 |       X,Y,Z Starting Location
!         +-----------+
!         +-----------+
!         | INT32     |       Zone (0=all)
!         +-----------+
!         +-----------+
!         | INT32     |       Color
!         +-----------+
!         +-----------+
!         | INT32     |       FillColor
!         +-----------+
!         +-----------+
!         | INT32     |       IsFilled (0=no 1=yes)
!         +-----------+
!         +-----------+
!         | INT32     |       GeomType  0=Line, 1=Rectangle 2=Square,
!         +-----------+                 3=Circle, 4=ellipse, 5=3DLine,
!         +-----------+
!         | INT32     |       LinePattern  0=Solid 1=Dashed 2=DashDot 3=Dotted...
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Pattern Length
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Line Thickness
!         +-----------+
!         +-----------+
!         | INT32     |       NumEllipsePts
!         +-----------+
!         +-----------+
!         | INT32     |       Arrowhead Style 0=Plain, 1=Filled, 2=Hollow
!         +-----------+
!         +-----------+
!         | INT32     |       Arrowhead Attachment 0=None, 1=Beg, 2=End, 3=Both
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Arrowhead Size
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Arrowhead Angle
!         +-----------+
!         +-----------+
!         | IN32*N    |       Macro Function Command (string: N = Length+1)
!         +-----------+
!         +-----------+
!         | INT32     |       Polyline Field Data Type 1=Float, 2=Double  (GTYPE)
!         +-----------+


! If the geometry type is line or 3dline then:
!         +-----------+
!         | INT32     |       Number of polylines
!         +-----------+
!         +-----------+
!         | INT32     |       Number of points, line 1.
!         +-----------+
!         +-----------+
!         | GTYPE*N   |       X-block geometry points N=NumPts
!         +-----------+
!         +-----------+
!         | GTYPE*N   |       Y-block geometry points N=NumPts
!         +-----------+
!         +-----------+
!         | GTYPE*N   |       Z-block geometry points N=NumPts (3dLine Only)
!         +-----------+
!             .
!             .
!             .

!If the geometry type is Rectangle then
!         +-----------+
!         | GTYPE*2   |       X and Y offset for far corner of rectangle
!         +-----------+

!If the geometry type is Circle then
!         +-----------+
!         | GTYPE     |       Radius
!         +-----------+

!If the geometry type is Square then
!         +-----------+
!         | GTYPE     |       Width
!         +-----------+

! If the geometry type is Ellipse then
!         +-----------+
!         | GTYPE*2   |       X and Y Radii
!         +-----------+

!    vi.   Text
!         +-----------+
!         | FLOAT32   |       Text marker.  Value=499.0
!         +-----------+
!         +-----------+
!         | INT32     |       Position CoordSys 0=Grid, 1=Frame, 2=FrameOffset(not used)
!         +-----------+
!         +-----------+
!         | INT32     |       Scope 0=Global 1=Local
!         +-----------+
!         +-----------+
!         | FLOAT64*2 |       X,Y Starting Location
!         +-----------+
!         +-----------+
!         | INT32     |       FontType
!         +-----------+
!         +-----------+
!         | INT32     |       Character Height Units 0=Grid, 1=Frame, 2=Point
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Height of characters
!         +-----------+
!         +-----------+
!         | INT32     |       Text Box type 0=NoBox 1=Hollow 2=Filled
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Text Box Margin
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Text Box Margin Linewidth
!         +-----------+
!         +-----------+
!         | INT32     |       Text Box Outline Color
!         +-----------+
!         +-----------+
!         | INT32     |       Text Box Fill Color
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Angle
!         +-----------+
!         +-----------+
!         | FLOAT64   |       Line Spacing
!         +-----------+
!         +-----------+
!         | INT32     |       Text Anchor. 0=left,      1=center,     2=right,
!         +-----------+                    3=midleft    4=midcenter   5=midright,
!                                          6=headleft   7=headcenter  8=headright
!         +-----------+
!         | INT32     |       Zone (0=all)
!         +-----------+
!         +-----------+
!         | INT32     |       Color
!         +-----------+
!         +-----------+
!         | INT32*N   |       MacroFunctionCommand (string: N = Length + 1)
!         +-----------+
!         +-----------+
!         | INT32*N   |       Text.  N=Text Length+1
!         +-----------+

!     vii.CustomLabel
!         +-----------+
!         | FLOAT32   |       CustomLabel Marker;  F=599
!         +-----------+
!         +-----------+
!         | INT32     |       Number of labels
!         +-----------+
!         +-----------+
!         | INT32*N   |       Text for label 1.  (N=length of label + 1) See note 1.
!         +-----------+
!         +-----------+
!         | INT32*N   |       Text for label 2.  (N=length of label + 1) See note 1.
!         +-----------+
!            .
!            .
!            .
!         +-----------+
!         | INT32*N   |       Text for label NumLabels.  (N=length of label + 1) See note 1.
!         +-----------+

!    viii.UserRec
!         +-----------+
!         | FLOAT32   |       UserRec Marker;  F=699
!         +-----------+
!         +-----------+
!         | INT32*N   |       Text for UserRec.  See note 1.
!         +-----------+

! II.  Data section
!     The data section contains all of the data associated with the
!     zone definitions in the header.

!     i.  zones
!         +-----------+
!         | FLOAT32   |       Zone marker  Value = 299.0
!         +-----------+
!         +-----------+
!         | INT32     |       Number of repeat vars.
!         +-----------+
!         +-----------+
!         | INT32*N   |       variable numbers to repeat.  N=Number of
!         +-----------+       variables to repeat.
!         +-----------+
!         | INT32*N   |       variable data format, N=Total number of vars
!         +-----------+       1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
!         +-----------+
!         | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
!         +-----------+       specified above.

!     ii. fem zones
!         +-----------+
!         | FLOAT32   |       Zone marker Value = 299.0
!         +-----------+
!         +-----------+
!         | INT32     |       Number of repeat vars.
!         +-----------+
!         +-----------+
!         | INT32*N   |       variable numbers to repeat.  N=Number of
!         +-----------+       variables to repeat.
!         +-----------+
!         | INT32*N   |       variable data format, N=Total number of vars
!         +-----------+       1=Float, 2=Double  DTYPE()
!         +-----------+
!         | xxxxxxxxxx|       Zone Data.  Each variable is in data format as
!         +-----------+       specified above.
!         +-----------+
!         | INT32     |       Repeat adjacency list from previous zone.
!         +-----------+       1=yes, 0=no
!         +-----------+
!         | INT32*N   |       Zone Data N=L*JMax.  This represents
!         +-----------+       JMax sets of adjacency indices where each
!                             set contains L values where L is
!                             3 for TRIANGLES
!                             4 for QUADRILATERALS
!                             4 for TETRAHEDRONS
!                             8 for BRICKS

! NOTES:

! 1.  All character data is represented by INT32 values.

!     Example:  The letter "A" has an ASCII value of 65.  The WORD
!               written to the data file for the letter "A" is then
!               65.
!               In fortran this could be done by doing the following:

!               Integer*32 I
!               .
!               .
!               I = ICHAR('A');
!               WRITE(10) I

!    All character strings are null terminated (i.e. terminated by a zero value)
! 2.  In FE Data I = Number of points, J = Number of elements, and
!    K = Element type where:

!    0 = Triangles;
!    1 = Quadrilaterals;
!    2 = Tetrahedrons.
!    3 = Bricks.

module mMyTecplt
 type tPlt
  integer unit
  integer tempunit
  character*80 dir
  character*80 filename
  integer byte_order
  character*100 title
  integer NumVar
  character*20, dimension(:), allocatable :: VarName
!zone info
  integer zne_format
  integer :: zne_color=-1
  integer zne_imax,zne_jmax,zne_kmax
  integer, dimension(:), allocatable :: zne_datafmt
 end type
 
 type(tPlt), dimension(:), allocatable :: pltfile 
 

character*8 :: MagicNum='#!TDV75'
real :: ZoneMarker=299.0
real :: GeoMarker=399.0
real :: EohMarker=357.0
integer :: ZeroMarker=0
integer :: NumRepVars=0

integer :: IsTec=0

CONTAINS

!initialize the binary file  1.i-1.iii
integer function myTecIni(filenum)
integer filenum

! integer errnum
integer i ! ,j,k

!open the binary file
!defaultfile=trim(dir),status='NEW'
open(UNIT=pltfile(filenum)%unit, FILE=trim(pltfile(filenum)%filename), FORM='BINARY', &
     ERR=1001)

write(pltfile(filenum)%unit) MagicNum,pltfile(filenum)%byte_order
call plt_write(pltfile(filenum)%unit, pltfile(filenum)%title)
write(pltfile(filenum)%unit) ZeroMarker,pltfile(filenum)%NumVar
do i=1, pltfile(filenum)%NumVar
 call plt_write(pltfile(filenum)%unit, pltfile(filenum)%VarName(i))
 write(pltfile(filenum)%unit) ZeroMarker
enddo

myTecini=0
1000 return
! trap open errors
1001 myTecini=-1
     return

end function

!write zone header  1.iv
!pltunit: file unit
! zne_format:
! 0=block, 1=point,2=FEBLOCK 3=FEPOINT
integer function  myTecZne(filenum,zne_name,zne_format,imax,jmax,kmax)
integer filenum
character(*) zne_name
integer zne_format
integer imax,jmax,kmax

integer :: zne_color=-1

integer i 

write(pltfile(filenum)%unit) ZoneMarker

call plt_write(pltfile(filenum)%unit, zne_name)
write(pltfile(filenum)%unit) ZeroMarker

write(pltfile(filenum)%unit) zne_format
pltfile(filenum)%zne_format=zne_format

write(pltfile(filenum)%unit) zne_color
pltfile(filenum)%zne_color=zne_color

write(pltfile(filenum)%unit) imax,jmax,kmax
pltfile(filenum)%zne_imax=imax
pltfile(filenum)%zne_jmax=jmax
pltfile(filenum)%zne_kmax=kmax

!write zone data section header to a temp file

open(UNIT=pltfile(filenum)%tempunit, FILE='plt_temp', FORM='BINARY', &
     ERR=1001)

write(pltfile(filenum)%tempunit) ZoneMarker
write(pltfile(filenum)%tempunit) NumRepVars
write(pltfile(filenum)%tempunit) (pltfile(1)%zne_datafmt(i),i=1,pltfile(1)%NumVar)


1000 myTecZne=0
     return

1001 myTecZne=-1
     return

end function


!write zone data section, to temp unit
integer function myTecDat(filenum,num_val,zne_data)
integer filenum
integer num_val
real  zne_data(num_val)

write(pltfile(filenum)%tempunit) zne_data

myTecDat=0

return
end function

!write zone data section, to unit
integer function mxTecDat(filenum,num_val,zne_data)
integer filenum
integer num_val
real  zne_data(num_val)

write(pltfile(filenum)%unit) zne_data

mxTecDat=0

return
end function

!get the file
integer function myTecEnd(filenum)
integer filenum

! integer i,j
integer :: ivar=0
integer a(256)

write(pltfile(filenum)%unit) EohMarker
rewind(pltfile(filenum)%tempunit) 

do while(ivar .eq. 0)
 read(pltfile(filenum)%tempunit,IOSTAT=ivar) a
 write(pltfile(filenum)%unit) a
end do

if (ivar .gt. 0) then
 myTecEnd=0
else
 myTecEnd=-1
endif
close(pltfile(filenum)%tempunit,status='DELETE')
close(pltfile(filenum)%unit)

return

end function

!cm lib
!used to write cm data
integer function  cmTecZne(filenum,zne_name,zne_format,imax,jmax,kmax)

integer filenum
character(*) zne_name
integer zne_format
integer imax,jmax,kmax

integer :: zne_color=-1

! integer i,j

write(pltfile(filenum)%unit) ZoneMarker

call plt_write(pltfile(filenum)%unit, zne_name)
write(pltfile(filenum)%unit) ZeroMarker

write(pltfile(filenum)%unit) zne_format
pltfile(filenum)%zne_format=zne_format

write(pltfile(filenum)%unit) zne_color
pltfile(filenum)%zne_color=zne_color

write(pltfile(filenum)%unit) imax,jmax,kmax
pltfile(filenum)%zne_imax=imax
pltfile(filenum)%zne_jmax=jmax
pltfile(filenum)%zne_kmax=kmax

1000 cmTecZne=0
     return

1001 cmTecZne=-1
     return

end function

!write EohMarker
integer function cmHeaderEnd(filenum)
integer filenum

write(pltfile(filenum)%unit) EohMarker

cmHeaderEnd=0
end function

!write zone data header
integer function cmZoneHeader(filenum)
integer filenum

write(pltfile(filenum)%unit) ZoneMarker
write(pltfile(filenum)%unit) NumRepVars
write(pltfile(filenum)%unit) (pltfile(1)%zne_datafmt(i),i=1,pltfile(1)%NumVar)


cmZoneHeader=1
end function

!write cm zone data 
!for test only , not used
integer function cmZoneData(filenum,i,j,k)
use paraset1
use paraset4
integer filenum
integer i,j,k

integer ini,inj,ink

if(j .eq. 2) then 
 ini=1
endif
do ink=1,zlevel(k)%cm_zlev(i,j)%fine(3)
  do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
   write(pltfile(filenum)%unit)  zlevel(k)%cm_zlev(i,j)%cm_x
  enddo
enddo

do ink=1,zlevel(k)%cm_zlev(i,j)%fine(3)
  do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
   do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1) 
    write(pltfile(filenum)%unit)  zlevel(k)%cm_zlev(i,j)%cm_y(inj)
   enddo
  enddo
enddo


do ink=1,zlevel(k)%cm_zlev(i,j)%fine(3)
  do inj=1,zlevel(k)%cm_zlev(i,j)%fine(2)
   do ini=1, zlevel(k)%cm_zlev(i,j)%fine(1) 
    write(pltfile(filenum)%unit)  zlevel(k)%cm_zlev(i,j)%cm_z(ink)
   enddo
  enddo
enddo

write(pltfile(filenum)%unit)  zlevel(k)%cm_zlev(i,j)%mat_matrix
if(zlevel(k)%cm_zlev(i,j)%s_flag .ge. 1 ) then
  write(pltfile(filenum)%unit)  zlevel(k)%cm_zlev(i,j)%src_matrix
else
  write(pltfile(filenum)%unit) (0.0, ini=1,zlevel(k)%cm_zlev(i,j)%fine(3)* &
              zlevel(k)%cm_zlev(i,j)%fine(2)*zlevel(k)%cm_zlev(i,j)%fine(1) )

endif

cmZoneData=0
end function

!close file
integer function  cmTecEnd(filenum)
integer filenum


close(pltfile(filenum)%unit)
cmTecEnd=0

end function

end module


!sample plt routine
subroutine testwrite
use mMyTecplt

integer tecflag,i
real point(3)

allocate(pltfile(1))
 pltfile(1)%unit=1
 pltfile(1)%tempunit=pltfile(1)%unit+111
 pltfile(1)%dir=''
 pltfile(1)%filename='testec3.plt'
 pltfile(1)%byte_order=1
 pltfile(1)%title='test tecplot'
 pltfile(1)%NumVar=3
 
 allocate(pltfile(1)%VarName(pltfile(1)%NumVar))
 pltfile(1)%VarName(1)='X'
 pltfile(1)%VarName(2)='Y'
 pltfile(1)%VarName(3)='Z'
 
 allocate(pltfile(1)%zne_datafmt(pltfile(1)%NumVar))
!1=float, 
 pltfile(1)%zne_datafmt=1
 
 tecflag=myTecIni(1)
 tecflag=myTecZne(1,'xyz',1,2,1,1)
 point(1)=1
 point(2)=2
 point(3)=3
 do i=1,2
   tecflag=myTecDat(1,3,point)
 enddo
 tecflag=myTecEnd(1)
 deallocate(pltfile)
end subroutine


!write a char as an integer
subroutine plt_write(unitnum, plt_string)
integer unitnum
character(*) plt_string

integer i
integer  plt_len, plt_char

plt_len=len_trim(plt_string)

do i=1, plt_len
  
  plt_char=ichar(plt_string(i:i))
  write(unitnum) plt_char

enddo

end subroutine