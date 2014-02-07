!read penmsh input files
!Generate inputfile for pentran
! main program
! complied by CVF 6.6 or ifc, pgf90 in windows 
! complied by ifc 8.0 pgf90 v6 or higher on linux

!********************************************

!Main Program starting here
!******************************************************************
program inpred
use paraset1
use paraset2
use paraset4
use files
use ErrControl

ver_num="version 2.73b (Feb 2014)"
!*****************************************
!control varibles in this section

IsProcessFlux=0
IsProcessCurr=0
IsDotOut=0 
IsOffTitan=0
!IsOffPentran=0

if( IsProcessFlux .eq. 1) then
 IsRefFlux=0
 IsCalPinPower=0
 IsPrbDotFlx(1)=0
 IsNomalizeFlux=0
endif

if(IsRefFlux .eq. 1) then
 IsPrbDotFlx(2)=0
endif


!IsSkipGeometry: 0, output geometry; 1, skip geometry output
IsSkipGeometry=0  
!IsProcessFlux: 0, no; 1, read heart atn and act input files instead of penmsh input

! ****************************************
! define input files directory

input_dir='.'
flux_dir='.'
curr_dir='.'
ref_flux_dir='.'

!unit READLOG, logfile
outputfile(1)%fullname='read.log'

!*****************************************
!read input files 

call GetArgument
call DisplayShortHeader

stime=SECNDS(0.0)

open(unit=READLOG, file=outputfile(1)%fullname)
write (READLOG,"('penmshxp ', A)" ) ver_num
write(form,"('read.log generated at : ')") 
call TimeStamp(READLOG)

write(*,*)  'processing inputfiles.....'
if(IsProcessHrt .eq. 1) then
!in supplement01 
 call ReadHrtInput
 
! stop 'give me a break, dude'
endif

if(IsProcessHrt .eq. 0) then
call ReadInput
 write(*,*)  'inputfiles done, logfile: ',trim(outputfile(1)%fullname)
 write(READLOG,*) 'PENMSH input files(.inp) completed'
 write(READLOG,*) '*************************************************'
 
 write(READLOG,*) 'Begin processing log'
 write(READLOG,*)

!***************************************
!processing input files in subcode2.f90
call Meshing
 write(*,*) 'Assigning Mat Num to fine meshes.....'
 call AssignMatNum
 write(*,*) 'Mat Num Assignment Done'
 write(READLOG,*) 'Assigning Mat Num to fine meshes.....done'
 if(s_format .ne. 2) then
   write(*,*) 'Assigning Src Mag. to fine meshes.....'
   call AssignSrcMag
   write(*,*) 'Src Mag Assignment Done'
   write(READLOG,*) 'Assigning Src Mag to fine meshes.....done'
 endif

endif
!***************************************
!processing flux if IsProcessFlux =1
if(IsProcessFlux .eq. 1) then
 write(*,*) 'allocating memory for flux...'
 call AllocateFluxArray
 write(*,*)  'allocation done'
 write(READLOG,*) 'allocating memory for flux...done'

 write(*,*) 'reading flux info...'
 call ReadInFluxset(1)
! call ReadInFlux
 write(*,*) 'flux done'
 write(READLOG,*) 'reading flux info...done'

 !avg. flux for each material
 if(IsFluxGM .ne. 0) call AvgMatFlux

 !flux output file
 if(IsFluxout .ne. 0 .and. num_group .eq. num_group_eff) &
  call MxOutputFlx
 !fov info for array project model
 if(IsFov .eq. 1) call WriteFov
 if(IsRefFlux .eq. 1) then
   write(*,*) 'reading Ref flux...'
   call ReadInFluxset(2)
   write(*,*) 'Ref flux done'
   write(READLOG,*) 'reading Ref. flux info...done'
 endif

endif

!***************************************
!processing current if IsProcessCurr =1
if(IsProcessCurr .eq. 1) then
 write(*,*) 'allocating memory for flux and current...'
 call AllocateFluxArray
 write(*,*)  'allocation done'
 write(READLOG,*) 'allocating memory for flux and current...done'

 write(*,*) 'reading flux&current info...'
 call ReadInCurrset(1)
! call ReadInFlux
 write(*,*) 'flux&current done'
 write(READLOG,*) 'reading flux&current info...done'
endif

!***************************************

!output geometry
if( IsSkipGeometry .eq. 0) then
  write(*,*) 'Dumping binary tecplot file...(-offplt option to turn off .plt file)'
  call MxOutputBinGemetry
  write(*,*) 'Geometry Done  ', trim(outputfile(3)%fullname)
  write(READLOG,*) 'Dumping binary tecplot file...done ', trim(outputfile(3)%fullname)
endif

!**************************************
! Write Dot f90 pentran input file
if (IsZSingle .eq. 0) then
if (IsProcessFlux .eq. 0) then

  if(IsOffPentran*IsOffTitan .eq. 0) then
    !assign transport parameters in subcode3.f90
     write(*,*) 'Assigning Transport Calculation Varibles...'
    call AssignPentranVaribles
    write(*,*) 'Assignment Done'
    write(READLOG,*) 'Assigning Pentran Varibles...done'
  
  endif
  
  if(IsOffPentran .eq. 0) then
    !Generate the pentran input file in subcode5.f90
    write(*,*) 'Generating pentran input deck...(-offf90 option to turn it off )'
    write(*,*) '(-nofido option to turn off fido in source section if it takes too long)'
    call WriteDotF90
    write(*,*) 'PENTRAN input file done ', trim(outputfile(4)%fullname)
    write(READLOG,*) 'Generating .pen pentran inputfile...done ',trim(outputfile(4)%fullname)
  else
    write(*,*) 'IsOffPentran=1: skipped generating PENTRAN input file' 
    write(READLOG,*) 'skipped generating PENTRAN input file'
  endif
  
  if(IsOffTitan .eq. 0) then
    !Generate the TITAN input file in subcode5.f90
    write(*,*) 'Generating TITAN input deck...(-offtan option to turn off TITAN input )'
    call WriteTitanInp
    write(*,*) 'TITAN input file done: ', trim(outputfile(6)%fullname)
    write(READLOG,*) 'Generating TITAN input file...done ',trim(outputfile(6)%fullname)
  else
    write(*,*) 'IsOffTitan=1: skipped generating TITAN input file' 
    write(READLOG,*) 'skipped generating TITAN input file'
  endif

elseif(IsProcessFlux .eq. 1) then
  write(*,*) 'skipped generating PENTRAN/TITAN input files' 
  write(READLOG,*) 'skipped generating PENTRAN/TITAN input files'
endif
endif !IsZSingle plot

!output material .out file 
if (IsDotOut .eq. 1) then
  write(*,*) 'generating .out file (penmsh outputfile) ...'
  call WriteDotOut
  write(*,*) '.out files done'
  write(READLOG,*) 'generating .out file (penmsh outputfile) ...done'
endif

if(IsPlotZlev .eq. 1) then 
  write(*,*) 'Plotting z levels...(use -offpng option to turn off plotting)'
  call OutputPlot
  write(*,*) 'z lev graphic file(s)  done '
  write(READLOG,*) 'Plotting z level...  z lev graphic file(s) done '
  
endif

write(*,*)
call SubWarning
write(*,*)
write(*,"('All done, logfile : ',A )") trim(outputfile(1)%fullname)
write(READLOG,"('All done, logfile : ',A )") trim(outputfile(1)%fullname)
write(READLOG,*) '***************************************************'
stime=SECNDS(stime)
write(*,"('Total Runtime : ', f7.1, 'sec' )") stime
write(READLOG,"('Total Runtime : ', f7.1, 'sec' )") stime

write(READLOG,*)
write(form,"('run finished at : ')")
call TimeStamp(READLOG)
close(READLOG)

end

!***********************************************************************
!end of main flow

! command -i input_dir -f flux_dir_rel
! flux_dir=./input_dir/flux_dir_rel
subroutine GetArgument
use paraset1
use files
use funs
use ErrControl

parameter(ndash=2)
character*80,dimension(:), allocatable :: buf
integer count
integer i
integer :: IsSuccess=0
!dash_index :  flag     flag_argument
!      1        -i      input dir
!      2        -f      [rel flux dir]
!      3        -l      [logfile name]
integer,dimension(:),allocatable :: id_arg  

real rarg
!count = iargc( )
count = IARGC( )
!buf(0): command itself
allocate(id_arg(0:count), buf(0:count))
id_arg=0
buf=''
call  GETARG(0,buf(0))

do i = 1, count
  
  call GETARG(i, buf(i))
  
  Flags: select case (trim(buf(i)))
  
    case ('-i','i','--input')
      id_arg(i)=1
    case ('-f', '-f1','-fl','--flux') 
      id_arg(i)=2
      IsProcessFlux=1
    case ('-l','l','--log')
      id_arg(i)=3
    case ('-dot')
      id_arg(i)=4 
    case('-w')
      id_arg(i)=5
    case('-nd','--ndmesh')
      id_arg(i)=6	
	case ('-f2')
      id_arg(i)=7
      IsRefFlux=1
    case ('-f3', '-fjn')
     IsProcessCurr=1
     id_arg(i)=8
	case ('-n')
     IsNomalizeFlux=1
     id_arg(i)=9

	 case ('-ff')   !generate prbname.flx.out : flux summary
      IsFluxout=1
      IsProcessFlux=1
      id_arg(i)=11
    case ('-fa')   !generate prbname.adj.out : adjoint flux summary
      IsFluxout=2
      IsProcessFlux=1
      id_arg(i)=11
	  
	case ('-chipath')   !generate prbname.adj.out : adjoint flux summary
      IsChiPath=1
      id_arg(i)=12
	  
    case ('-offf90','-offpen')
      IsOffPentran=1
    case ('-offtan', '-offtitan')
      IsOffTitan=1
    case ('-prbflx','-prb')
      IsPrbDotFlx=1
    case('-nofido')
      IsFidoSrc=0
    
	case ('-offplt')
      IsSkipGeometry=1
    case ('-offpng')
      IsPlotZlev=0
    
	case ('-ming')
      flx_min_flag=2
    case ('-maxg')
      flx_max_flag=2
    
	case ('-min')
      id_arg(i)=28
      flx_min_flag=1
    case ('-max')
      id_arg(i)=29
      flx_max_flag=1

    case ('-plot', '-plotzyx')
      id_arg(i)=30
    case ('-xplot', '-plotx', '-xp')
      id_arg(i)=31
      IsXPlot=1
    case ('-yplot','-ploty', '-yp')
      id_arg(i)=32
      IsYPlot=1
    case ('-zplot','-plotz', '-zp')
      id_arg(i)=33
      IsZPlot=1
    case ('-pngsize','-png', '-size')
      id_arg(i)=34
    case('-plotmsf','-msf')
      id_arg(i)=35
    case ('-z')
      IsZSingle=1
      id_arg(i)=36
    case ('-color')
      id_arg(i)=37
    case ('-plot3d', '-3d','-3D')
      id_arg(i)=38
      Is3DPlot=1
    case ('-3dpos')
      id_arg(i)=39
      Is3DPlot=1
    case ('-nofm')
      IsDrawFm=0
    case('-log')
      IsLogPlot=1
    case ('-titan','-t')
      IsOffTitan=0
    case ('-gm', '-fgm')
      IsFluxGM=1
      IsProcessFlux=1
    case ('-agm')
      IsFluxGM=2
      IsProcessFlux=1
    case ('-fov')
      IsFov=1
 !for binary heart-like input    
    case ('-hrt')
      IsProcessHrt=1
      id_arg(i)=40
 
    case ('-V', '-v')
      call ShowVersion

    case ('help','-h','-H','h','H','--help')
      call DisplayHelp
    case ('colormap','-helpcm','helpcm')
      call DisplayColorMap
 
!support for fm size-based nmesh allocation
    case ('-maxmesh')
	  id_arg=41

    case default
        
      Flag_arg: select case (id_arg(i-1))
       
      case(1)  !-i
!         input_dir='./'//trim(buf(i))
         input_dir=trim(buf(i))
      case(2)  !-f1
         flux_dir=trim(buf(i))
      case(3)  !logfile name
         outputfile(1)%fullname=trim(buf(i))
      case(4)
         IsSuccess=CharNum(buf(i), DotBit)
         if(IsSuccess .eq. 0 ) then
            DotBit=-5
            write(*,"('GetArg: -dot flag argument fail','ignored -dot ', A,' using default Dotbit=-5 ')") &
              trim(buf(i))
         endif
      case(5)
         IsSuccess=CharNum(buf(i),max_warning)
         if(IsSuccess .eq. 0 ) then
            max_warning=5
            write(*,"('GetArg: -w flag argument fail, ignored -w ', A, ' using default max_warning=5')")& 
              trim(buf(i))
         endif
      case(6)
         IsSuccess=CharNum(buf(i),ndmeth_global)
         if(IsSuccess .eq. 0 ) then
            ndmeth_global= 2
           write(*,"('GetArg: -nd flag argument fail, ignored -nd ',A, ' using default ndmeth_global=2 ')")& 
              trim(buf(i))
        endif 
      case(7) !-f2
         ref_flux_dir=trim(buf(i))
      case(8) !-f3
         curr_dir=trim(buf(i))
      case(9)
	    IsSuccess=CharFloat(buf(i),n_factor)
        if(IsSuccess .eq. 0 ) then
		  nfile%fullname=trim(buf(i))
		  n_factor=-1.0
		elseif (n_factor .lt. 0) then
          write(*,"('GetArg: -n flag argument fail, normalization factor has to be positive ',A, ' using default n_factor=max global flux ')")& 
            trim(buf(i))
		  n_factor=0.0
		endif
	  
	  case(11)
        IsSuccess=CharNum(buf(i),num_flx_out)
        if(IsSuccess .eq. 0 ) then
          num_flx_out=1
          write(*,"('GetArg: multi volumn flux file off ')")
        endif
       
	  case(12)
	    ChiPath_dir=trim(buf(i))
	  case(28)
	    IsSuccess=CharFloat(buf(i),flx_min)
		flx_min_flag=3
        if(IsSuccess .eq. 0 ) then
         write(*,"('GetArg: -min flag argument not present ',A, 'using default = min global flux')") & 
           trim(buf(i))  
		 flx_min_flag=1
!		elseif (flx_min .lt. 0) then
!          write(*,"('GetArg: -min flag argument fail, mininum flux has to be positive, ignored ',A)")& 
!            trim(buf(i))
!		  flx_min=0.0
!		  flx_min_flag=1
		endif
      case(29)
	    IsSuccess=CharFloat(buf(i),flx_max)
		flx_max_flag=3
        if(IsSuccess .eq. 0 ) then
         write(*,"('GetArg: -max flag argument not present ',A,' using default = max global flux ')") & 
           trim(buf(i))  
		   flx_max_flag=1
!		elseif (flx_max .le. 0) then
!          write(*,"('GetArg: -max flag argument fail, maxium flux has to be positive, ignored ',A, ' using default n_factor=max global flux ')")& 
!            trim(buf(i))
!		  flx_max=0.0
!		  flx_max_flag=1
		endif
	  case(30)
        IsSuccess=CharNum(buf(i),MidPlot)
        if(IsSuccess .eq. 0 ) then
          MidPlot=1
          write(*,"('GetArg: -plot flag argument fail, ignored ',A, ' using default plotting z levels ')")& 
          trim(buf(i))
        endif
      case (31)
        call Str2Val(buf(i), rarg, IsSuccess)
        if(IsSuccess .eq. 0) then 
          num_xloc=num_xloc+1
          id_arg(i)=31
          if(num_xloc .gt. maxplot) then
            write(*,"('GetArg: -xplot flag argument, too many positions, after ', I0, ' ignored')")   maxplot
          else
            xyzloc(num_xloc,1)=rarg
          endif 
       else
          write(*,"('GetArg: -xplot flag argument, unregnized number: ',A, ' ignored')")   trim(buf(i))
       endif
     case (32)
       call Str2Val(buf(i), rarg, IsSuccess)
       if(IsSuccess .eq. 0) then 
         num_yloc=num_yloc+1
         id_arg(i)=32
         if(num_yloc .gt. maxplot) then
           write(*,"('GetArg: -yplot flag argument, too many positions, after ', I0, ' ignored')")   maxplot
         else
           xyzloc(num_yloc,2)=rarg
         endif
       else
         write(*,"('GetArg: -yplot flag argument, unregnized number: ',A, ' ignored')")   trim(buf(i))
       endif
     case (33)
       call Str2Val(buf(i), rarg, IsSuccess)
       if(IsSuccess .eq. 0) then 
         num_zloc=num_zloc+1
         id_arg(i)=33
         if(num_zloc .gt. maxplot) then
           write(*,"('GetArg: -zplot flag argument, too many positions, after ', I0, ' ignored')")   maxplot
         else
           xyzloc(num_zloc,3)=rarg
         endif
       else
          write(*,"('GetArg: -zplot flag argument, unregnized number: ',A, ' ignored')")   trim(buf(i))
       endif
    case(34)
      call Str2Val(buf(i), rarg, IsSuccess)
      if(IsSuccess .eq. 0) then 
         if(rarg .gt. 9) then
            write(*,"('GetArg: -plotsize factor too big, using 9 ')")   
     
         elseif(rarg .lt. 0.1) then
           write(*,"('GetArg: -plotsize factor too small , using 0.1 ')")
         else
           PngSizeFactor=rarg
         endif
      else
        write(*,"('GetArg: -pngsize flag argument, unregnized number: ',A, ' ignored')")   trim(buf(i))
      endif
    case (35)
     IsSuccess=CharNum(buf(i),MatSrcFlx)
     if(IsSuccess .eq. 0 ) then
       MatSrcFlx=1
       write(*,"('GetArg: -plot flag argument fail, ignored ',A, ' using default plotting z levels ')")& 
       trim(buf(i))
     endif
     if (abs(MatSrcFlx) .ge. 4) IsProcessFlux=1
   case (36)
       IsSuccess=CharNum(buf(i),z_start)
       z_end=z_start
       if(IsSuccess .eq. 0 ) then
          z_start=1
          z_end=1
          write(*,"('GetArg: -z option without an argument, plotting the frist z level only ')")
       endif
   case(37)
     IsSuccess=CharNum(buf(i),ColorMap)
     if(IsSuccess .eq. 0 ) then
        ColorMap=3
        write(*,"('GetArg: invalid colormap id, use <penmsh colormap> for more help, set to default')")
     endif
   case (38) 
      IsSuccess=CharFloat(buf(i),view_ang)
	  if(Issuccess .eq. 0) then
	  view_ang=0
	  write(*,"('GetArg: -3d option without an argument, 3d plotting from 0 degree angle (default) ')")
      endif
   case (39)
    call Str2Val(buf(i), rarg, IsSuccess)
	if(IsSuccess .eq. 0) then 
      num_vpos=num_vpos+1
      id_arg(i)=39
      if(num_vpos .gt. 3) then
        write(*,"('GetArg: -3dpos flag argument, too many values, first three entries are used')")   
      else
       view_pos(num_vpos)=rarg
      endif
    else
      write(*,"('GetArg: -3dpos without argument, using view point (3,3,3) default')")   
    endif
 
   case (40)
       hrt_filename=trim(buf(i))  
   case (41)
      call Str2Val(buf(i), rarg, IsSuccess)
      if(IsSuccess .eq. 0) then 
         if(rarg .gt. 0) then
          max_fm_size=rarg
         else
           max_fm_size=1.0
         endif
      else
        write(*,"('GetArg: -maxfm is not valid: ',A, ' ignored')")   trim(buf(i))
      endif

   case default
         write(*,*) "Warning: ", trim(buf(i)), ": unknown argument ignored"
         write(*,*) "Use inpred -h to display help"
         stop
   end select Flag_arg

  end select  Flags  
  
enddo

if (IsProcessFlux .eq. 0 .and. IsRefFlux .eq. 1) then
  write(*,*) "Commandline option conflict: ",  "-f2  option used only if -f1 presents"
  write(*,*) "Use inpred -h to display help"
  stop
endif

if (IsProcessFlux .eq. 0 .and. IsRefFlux .eq. 1) then
  write(*,*) "Commandline option conflict: ",  "-f2  option used only if -f1 presents"
  write(*,*) "Use inpred -h to display help"
  stop
endif

if (IsZSingle .eq. 1) then

if(IsProcessHrt .eq. 1) then
  write(*,*) "Commandline option conflict: ",  "-z  and -hrt'"
  write(*,*) "Use inpred -h to display help"
  stop
endif

endif

If(IsProcessCurr .eq. 1) then
 if(IsProcessFlux .ne. 0 .or. IsRefFlux .ne. 0) then
  write(*,*) "Commandline option conflict: ",  "-f3 (fjn) option is used, -f1 or -f2 will be ignored "
  write(*,*) "Use inpred -h to display help"
  IsProcessFlux=0
  IsRefFlux=0
 endif
endif


end subroutine

!Display the help screen and terminate the program
subroutine DisplayHelp
use paraset1, only: ver_num

write(*,"(A)") "NAME"
write(*,"(A)") 
write(*,"(A)") "PENMSHXP - generate PENTRAN inputfile"
write(*,'("           (inpred)  " ,A)')  trim(ver_num)
write(*,"(A)") "   input files : penmsh.inp, prb_name[z_level#].inp"
write(*,"(A)") "      optional : prb_name[grp#].flx or prb[grp#].flx : "
write(*,"(A)") "                    pentran serial/parallel  flux files "
write(*,"(A)") "                 prb_name.spc, prb_name.chi   : "
write(*,"(A)") "                    src spectrum and fission spectrum file "
write(*,"(A)") "                 prb_name.mba (mat balance file) "
write(*,"(A)") "                 prb_name.src (src density grid file) "
write(*,"(A)")
write(*,"(A)") "  output files : prb_name_out.pen (PENTRAN input); "    
write(*,"(A)") "                 prb_name_mix.plt (tecplot data file) "
write(*,"(A)") "                 prb_name.mcr (tecplot macro file) " 
write(*,"(A)") "                 prb_name_out.mba (Material vol balance ) "    
write(*,"(A)") "                 read.log         (read log file        ) "    
write(*,"(A)")
write(*,"(A)") "SYNOPSIS"
write(*,"(A)")
write(*,"(A)") "    penmshxp [flag1] [flag1 argument] [flag2] [flag2 argument] ..." 
write(*,"(A)")
write(*,"(A)") "DESCRIPTION"
write(*,"(A)")
write(*,"(A)")  'Flag          Flag_argument'
write(*,"(A)")
write(*,"(A)")  '-i,--input    input_dir'
write(*,"(A)")  '       input file directory, default current dir '
write(*,"(A)")
write(*,"(A)")  '-f1,-f,--flux     flux_dir'
write(*,"(A)")  '       flux file dir, default: current dir '
write(*,"(A)")  '       ./[flux_dir]'
write(*,"(A)")  '-f2           ref_flux_dir'
write(*,"(A)")  '       reference flux file dir, default: current dir '
write(*,"(A)")  '       ./[flux_dir]'
write(*,"(A)")  '-f3, -fjn         curr_dir'
write(*,"(A)")  '       flux current file dir, default: current dir '
write(*,"(A)")  '       ./[flux_dir]'
write(*,"(A)")  ' Note: .fjn files are geneated by PENDATA (Option 9) with 5 rows only:  '
write(*,"(A)")  '       Phi   J    Jx    Jy     Jz (Select Varibles No.10-14 in PENDATA) '
write(*,"(A)")  '-n     <factor or a filename  > '
write(*,"(A)")  '   To nomalized fluxes flux=flux/factor, if no argument given, normalized to global max flux' 
write(*,"(A)")  '   filename: a file containing a factor array for each group (fido supported) ' 
write(*,"(A)")  '             size=num_group, in increasing order (fwd) or decreasing order (adj) '
write(*,"(A)")  '   e.g. -n 100.0  all fluxes are divided by 100,'
write(*,"(A)")      
write(*,"(A)")  '-l,--log      logfile name'
write(*,"(A)")  '       log file name' 
write(*,"(A)")  '-dot   DotBit   '
write(*,"(A)")  '       FIDO control on spacef value format'
write(*,"(A)")  '       DotBit=3, i.e 45.123, 333.123   '
write(*,"(A)")  '       DotBit=-5 i.e  1.12345E-01 ,default        ' 
write(*,"(A)")  '-nd --ndmeth  ndmeth   '
write(*,"(A)")  '       GLOBLE ndmeth varible: differencing scheme '
write(*,"(A)")  '       ndmeth=2 : default, adaptive DTW, upgradable to EDW   '
write(*,"(A)")  '       ndmeth=-2  fixed DTW       ' 
write(*,"(A)")  '-offplt                         '
write(*,"(A)")  '     Skip tecplot Geometry file  '
write(*,"(A)")  '-offpng                         '
write(*,"(A)")  '     Skip z level plots   '
write(*,"(A)")  '-plotmsf -msf MatSrcFlx                        '
write(*,"(A)")  '    MatSrcFlx=1, 2 or 4: plotting mat, src, or flx  '
write(*,"(A)")  '    default=1 , e.g. -plotmsf 5 will plot mat and flx  '
write(*,"(A)")  '    e.g -plotmsf -4 will plot flux, and write flux to a .cat file  '
write(*,"(A)")  '-plotzyx MidPlot                        '
write(*,"(A)")  '    MidPlot=1,2 or 4: plotting z, y, or x mid-levels  '
write(*,"(A)")  '    default=1 , e.g. -plotzyx 3 will plot z and y  '
write(*,"(A)")  '-plotx  pos1 pos2...                         '
write(*,"(A)")  '    x-y plot only at z= <pos1 pos2 ...>     '
write(*,"(A)")  '-ploty  pos1 pos2 ...                         '
write(*,"(A)")  '    x-z plot only at y= <pos1 pos2 ...>     '
write(*,"(A)")  '-plotz  pos1 pos2...                         '
write(*,"(A)")  '    y-z plot only at x= <pos1 pos2 ...>     '
write(*,"(A)")  '-size  factor                         '
write(*,"(A)")  '    png file resolution multiplication factor     '
write(*,"(A)")  '-3d  view_angle                         '
write(*,"(A)")  '    generate 3d plots, and set view angle in degree     '
write(*,"(A)")  '    default angle is 0 degree and top view, negative entry for bottom view  '
write(*,"(A)")  '-3dpos  x  y  z                         '
write(*,"(A)")  '    generate 3d plots, and set viewpoint position at (x, y, z) '
write(*,"(A)")  '    default position is (4 4 4),  plotting area coordinates  '

write(*,"(A)")  '-max  <flux max>                         '
write(*,"(A)")  '    flux plot scale : maxium, default max among all groups   '

write(*,"(A)")  '-min  <flux min>                         '
write(*,"(A)")  '    flux plot scale : minmum, warning lower than the minmum points will black out '

write(*,"(A)")  '-maxg                          '
write(*,"(A)")  '    flux plot scale : max for one group   '

write(*,"(A)")  '-ming                           '
write(*,"(A)")  '    flux plot scale : minmum for one group '

write(*,"(A)")  '-nofm                         '
write(*,"(A)")  '    turn off Fine Mesh lines in all coarse meshes   '
write(*,"(A)")  '    Note: FM lines automatically off if FM size < 8 pixels in a CM '
write(*,"(A)")  '-offf90                         '
write(*,"(A)")  '    Skip generate pentran input file  '
write(*,"(A)")  '-z  z_level_num                       '
write(*,"(A)")  '    plot one z-level only  '
write(*,"(A)")  '-color ColorMap number'
write(*,"(A)")  '   type penmshxp -helpcm for more help'


write(*,"(A)")  '-fgm                         '
write(*,"(A)")  '    generate file flux.fgm : avg. flux per material zone for each group'
write(*,"(A)")  '    Note: require load the flux using -f option '
write(*,"(A)")  '-agm                         '
write(*,"(A)")  '    generate file flux.agm: avg. adjoint flux per mat zone for each group'
write(*,"(A)")  '    Note: groups are flipped  '
write(*,"(A)")  '-ff                         '
write(*,"(A)")  '    generate file prbname.flx.out : flux distribution file'
write(*,"(A)")  '    Note: require load the all group flux using -f option'
write(*,"(A)")  '-fa                         '
write(*,"(A)")  '    generate file prbname.adj.out : adjoint flux distribution file '
write(*,"(A)")  '    Note: groups are flipped  '

write(*,"(A)")  '-titan                         '
write(*,"(A)")  '     generate titan input file  '
!write(*,"(A)")  '-offtitan                         '
!write(*,"(A)")  '     Skip generate titan input file  '
write(*,"(A)")  '-nofido                         '
write(*,"(A)")  '     Turn off FIDO generation in source section of PENTRAN input deck '
write(*,"(A)")  '-w   number                     '
write(*,"(A)")  '     Max warnings                '
write(*,"(A)")  '-hrt   phantom input file name    '
write(*,"(A)")  '    handle binary phantom data files                 '
write(*,"(A)")  '-maxmesh number '
write(*,"(A)")  '     automatic number of meshes per coarse mesh, based on a maximum '
write(*,"(A)")  '     allowable mesh size of (number) this number is then multiplied '
write(*,"(A)")  '     by the number in the prbZ.inp file to allow for variable mesh size '
write(*,"(A)")  '-V       '
write(*,"(A)")  '    show version number                 '
write(*,"(A)")  

stop  'Report bugs to CE YI<yice@gatech.edu>'

end subroutine

subroutine DisplayColorMap

write(*,"(A)") 'Use -color <colormap id> option to set a color map'
write(*,"(A)") 'colormap id: '
write(*,"(A)") ' 0 : dislin id=SMALL, defines a small colour table with the 8 colours'
write(*,"(A)") '     1 = BLACK, 2 = RED, 3 = GREEN, 4 = BLUE, 5 = YELLOW, 6 = ORANGE, 7 = CYAN and 8 = MAGENTA'
write(*,"(A)") ' 1 : dislin id= VGA, defines the 16 standard colours of a VGA graphics card'
write(*,"(A)") ' 2 : dislin id=RAIN, defines 256 colours arranged in a rainbow where 0 means black and 255 means white'
write(*,"(A)") ' 3 : dislin id=SPEC, defines 256 colours arranged in a rainbow where 0 means black and 255 means white' 
write(*,"(A)") '     This colour table uses more violet colours than RAIN'
write(*,"(A)") ' 4 : dislin id=GREY, defines 256 grey scale colours where 0 means black and 255 is white'
write(*,"(A)") ' 5 : dislin id=RRAIN, is the reverse colour table of RAIN '
write(*,"(A)") ' 6 : dislin id=RSPEC, is the reverse colour table of SPEC'
write(*,"(A)") ' 7 : dislin id=RGREY, is the reverse colour table of GREY'
write(*,"(A)") ' 8 : dislin id=TEMP, defines a temperature colour table'
write(*,"(A)") ' The default colour table is 3 (SPEC) '
write(*,"(A)")

stop  'Report bugs to CE YI<yice@gatech.edu>'

return
end subroutine
!Display the help screen and terminate the program
subroutine DisplayHeader
use paraset1, only: ver_num

write(*,"(A)") 
write(*,"(A)") 
write(*,"(A)") "                              PENMSHXP                             " 
write(*,"(A,A)") "                PENMSH Express ", trim(ver_num)     
write(*,"(A)") "           A Mesh Generator to Build PENTRAN Input Deck            "
write(*,"(A)") "     with compatibility to PENMSH developed by Dr. A. Haghighat    " 
write(*,"(A)") "            (DISLIN graphic library by Helmut Michels)             " 
write(*,"(A)") 
write(*,"(A)")
write(*,"(A)") "                               CE YI                              " 
write(*,"(A)") "                        email: yice@gatech.edu                    "
write(*,"(A)")
write(*,"(A)")
write(*,"(A)") "                    Use 'penmshxp -h' for more info               "
write(*,"(A)")
write(*,"(A)")
write(*,"(A)")
write(*,"(A)")
write(*,"(A)")  'Begin processing...'


end subroutine

!Display the help screen and terminate the program
subroutine DisplayShortHeader
use paraset1, only: ver_num


write(*,"('PENMSH ', A)")  trim(ver_num)
write(*,"('use -h option for help, and check the read.log file')")     
write(*,"(A)") 
write(*,"(A)")  'Begin processing...'


end subroutine

!Display Version and copyright info
Subroutine ShowVersion
use paraset1, only : ver_num

write(*, "(' PENMSHXP ', A)")  ver_num
write(*, "(' Copyright 2011 Ce Yi at GaTech')" )
if( index(ver_num, 'b ') .ne. 0)   &
write(*, "(' DISLIN Graphic Library by Helmut Michels at Max Planck Institute' )") 
stop '      '
end subroutine

!put a tiem stamp on unit_num
!if unit_num=0, stamp on the screen
!a message delivered too at the same time
!the message is stored in 'form' @ module 'files'
subroutine TimeStamp(unit_num)
use files

integer unit_num



character(LEN =8) str_date
character(LEN =12) str_time

call date_and_time (str_date, str_time)

if(unit_num .eq. 0) then

write(*, "(A, 4x,A2,':',A2,':',A2,2x,'on ',A2,'/',A2,'/', A4)" )  trim(form), &
   str_time(1:2),str_time(3:4),str_time(5:6),str_date(5:6),str_date(7:8), str_date(1:4)
write(*,*)

else
write(unit_num, "(A, 4x,A2,':',A2,':',A2,2x,'on ',A2,'/',A2,'/', A4)" )  trim(form), &
   str_time(1:2),str_time(3:4),str_time(5:6),str_date(5:6),str_date(7:8), str_date(1:4)
write(unit_num,*)
endif

end subroutine

!convert string to value
subroutine Str2Val(str, val, ierr)

use ErrControl

character(*) str
real val
integer ierr

read(str,*, iostat=ierr) val
!if(ierr .ne. 0) then
!write(err_message,"('error in converting string to val ')") 
 
!call TrapInputError(1)
!endif

return
end subroutine
