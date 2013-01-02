! subcode6.f90   Generating FIDO array
! subroutines in this file
!          
!          FidoIntArray:  Fido out an integer array
!          FidoRealArray: Fidoing out a real array
!          
!

! write online1(fido output array) to file unit=filehandler
! num_space: number of space written before the array
! with prevar
! be careful  oneline1 is deallocated after the sub
subroutine fidowriting(filehandler,num_space)
use fido
integer filehandler, num_space

character, dimension(:), allocatable :: oneline2  !, inline
! character (len=20) :: line_last='',line_second='' 

!character :: bufline(78)=''
!integer :: num_line
!integer :: bufleng

integer last,lastspace,line_pointer, i,j
integer len1,lenadd, lenadd_pre, lenadd_app


lenadd_pre=len_trim(prevar)
lenadd_app=len_trim(appvar)
lenadd=lenadd_pre+lenadd_app+num_space 
len1=output_len+lenadd

!num_line=0
!do while(len1 .ne. 0)
! num_line=num_line+1
! if(num_line .eq. 1) then
!   bufline(1:lenadd_pre)=prevar(1:lenadd_pre)
!   do j=1, num_space
!    bufline(lenadd_pre+j)=' '  
!   enddo
! endif
!end do

if( lenadd .ne. 0) then
 allocate(oneline2(output_len))
 oneline2(1:output_len)=oneline1(1:output_len)
 deallocate(oneline1)
 allocate(oneline1(len1))
 do i=1,len1
   
   if(i .ge. 0 .and. i .le. lenadd_pre) then 
      oneline1(i)=prevar(i:i)
   else if(i .ge. (lenadd_pre+1) .and. i .le. (lenadd_pre+num_space) ) then
      oneline1(i)=' '
   else if(i .ge. (lenadd_pre+num_space+1) .and. i .le. (lenadd_pre+num_space+output_len) ) then
     oneline1(i)=oneline2(i-lenadd_pre-num_space)
   else if(i .ge. (lenadd_pre+num_space+output_len+1) .and. i .le. len1) then
    oneline1(i)=appvar(i-(len1-lenadd_app):i-(len1-lenadd_app))
   endif
 enddo
 deallocate(oneline2)
endif

line_pointer=0
do while(len1 .ne. 0)                         
   if(len1 .le. 78) then
     write(filehandler,'(78A)') (oneline1(line_pointer+j),j=1,len1)
     len1=0
   else
     last=79   !if oneline1(81)=' ', it will be omitted. (next line will not begin with a space)
     do while(oneline1(line_pointer+last) .ne. ' ')
      last=last-1
     enddo

!to avoid 'T'-alone line
! len1 : length of the remainder of oneline1
! last : last space location for current line
	 if(len1-last .gt. 0) then
     j=line_pointer+last+1 
	 !'T' -alone line
     if(oneline1(j) .eq. 'T' ) then
        last=last-1
	   !skip space
	   do while(oneline1(line_pointer+last) .eq. ' ')
        last=last-1
       enddo
	   !skip the last fido entry, which will be moved to the next line, which otherwise would be a T-alone line
       do while(oneline1(line_pointer+last) .ne. ' ')
        last=last-1
       enddo

	 endif
	 endif


	 lastspace=last
     write(filehandler,'(79A)') (oneline1(line_pointer+j),j=1,min(79,lastspace))
  
  if(len1-lastspace .gt. 0) then
   line_pointer=line_pointer+lastspace
!   allocate(inline(len1-lastspace))
!   do j=1,len1-lastspace
!        inline(j)=oneline1(lastspace+j)
!   enddo
!   deallocate(oneline1)
!   allocate(oneline1(len1-lastspace))
!      do j=1,len1-lastspace
!        oneline1(j)=inline(j)
!   enddo
!   deallocate(inline)
   len1=len1-lastspace
  else
   len1=0
  endif   
   
   endif
enddo

!to avoid 'T'-alone line
!if(line_pointer .ne. 0 .and. index(appvar,'T') .ne. 0 ) then  !avoid

!backspace(filehandler)
!read(filehandler,"(A)" ) line_last
!line_last=adjustl(line_last)
!if(line_last(1:1) .eq. 'T') then
!   backspace(filehandler)
!   backspace(filehandler)
!   read(filehandler,"(A)" ) line_second
!   last=len_trim(line_second)  
!   len1=last
!   do while(line_second(last:last) .ne. ' ')
!     last=last-1
!    enddo
!   if (last .ge. 1 .and. last .lt. len1 ) then
!     backspace(filehandler) 
!     write(filehandler,"(A)" ) line_second(1:last)
!     write(filehandler,"(A,2x, A)" ) line_second(last+1:len1), line_last(1:len_trim(line_last))
!!forward one lines
!   else
!     read(filehandler,*) 
!   endif
!endif
!endif   !avoid

deallocate(oneline1) 

end subroutine


!generating one FIDO input
! input: fidonum1, fidochar2, fidomask, output: fidoset, lenset
!fidnum1: integer; fidochar2: character; fidochar: character
!e.g.  if fidonum1=10, fidochar2='-123.45    ', fidomask='R' then
!      fidoset='10R-123.45            ', lenset=10, and
! fidoset(1:lenset)='10R-123.45'  meaning repeating -123.45 10 times.

subroutine fido_generating(fidonum1, fidochar2,fidomask, fidoset,lenset)
use fido
integer fidonum1,lenset
character*(*) fidochar2,  fidoset 
character*1 fidomask

integer len1, len2

call intlen(fidonum1,len1)
fidochar2=adjustl(fidochar2)
len2=len_trim(fidochar2)
lenset=len1+len2+1
if (len(fidoset) .ge. lenset) then
  write(fidoform,"('(I',I0,',A,A',I0,')' )") len1,len2 
  write(fidoset, fidoform) fidonum1,fidomask, fidochar2
else
  stop 'trouble in subroutine fido_generating, length for fidoset too small' 
endif

num_fido_char=num_fido_char+1

return
end subroutine  

!calculating the length of an integer
!input: i, output: length
!e.g.  if i=-1234 then length=5
subroutine intlen(i, length)
integer i, length

character*20  tempchar

write(tempchar, '(I20)') i
length=len_trim(adjustl(tempchar))

return
end subroutine

!generating fido array from an real input array
! be careful: before calling this sub, the input array should be loaded into 
! input array(input_array), after this call, data in fido output array(oneline1,global) 
! should be extracted, after that, deallocating oneline1 is recommended(not necessary).
! or during the next call, the old data in oneline1 will be lost.
! typeflag: id of output format
! typeflag=positive or 0: f<total>.<typeflag>   
! typeflag=negtive      : ES<total>.<typeflag> 
! e.g.  typeflag=3  a fido output would be: 12R123.123
!       typeflag=0  a fido output would be: 12R123.
!       typeflag=-3 a fido output would be: 12R1.231E+02
! qfidoflag: same as the subroutine fidoarray_int
subroutine FidoRealArray(input_array, input_len,typeflag,qfidoflag, rfidoflag)
use fido
use funs
use ErrControl

!optional :: rfidoflag
integer input_len, typeflag, qfidoflag, rfidoflag
real input_array(input_len)

integer i,j,m
integer sum, len1

integer status
real point_pre, point_cur
integer i_value
integer point_bit, lenadd, lenset
character*20 :: fidochar2=''
character*60 :: fidoset=''
character, dimension(:), allocatable :: oneline2

integer len_before, len_after
integer scan, scan_scope, scan_begin
integer qbingo,pos,flag, q_repeat,qlength,afterq
integer typeflag_in, ES_bit
integer skip, scan_skip, qfidoflag_in
integer, dimension(:), allocatable :: qposition, qposition_temp

if(input_len .gt. max_input_len) max_input_len=input_len
num_fido_char=0

i=1
sum=0
len1=0
skip=0
qfidoflag_in=abs(qfidoflag)
typeflag_in=abs(typeflag)

allocate(oneline1(len1), stat=status)
if( status .ne. 0) then
 deallocate(oneline1)
 allocate(oneline1(buf_size))
endif

!if(present(rfidoflag) .and. rfidoflag .eq. 0) then
if(rfidoflag .eq. 0) then
do m=1, input_len

   point_cur=input_array(m)
   if(point_cur .lt. 0) then
     ES_bit=7
   else
     ES_bit=7
   endif
   if (typeflag .gt. 0) then
    i_value=int(point_cur)
    call intlen(i_value, point_bit)
    point_bit=point_bit+typeflag+1
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform,"('(f',I0,'.',I0,')' )") point_bit, typeflag
    write(fidochar2,fidoform) point_cur
   elseif(typeflag .eq. 0) then
    i_value=nint(point_cur)
    call intlen(i_value, point_bit)
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidochar2,'(I0)') i_value
   else
    point_bit=abs(typeflag)+Es_bit
 
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform, "('(ES',I0,'.',I0,')' )") point_bit,typeflag_in 
    write(fidochar2,fidoform) point_cur 
   endif
   fidochar2=adjustl(fidochar2)

    lenadd=1+point_bit            !plus a space  
!   allocate(oneline1(len1+lenadd))
!   do j=1,len1
!   oneline1(j)=oneline2(j)
! enddo
    if( (size(oneline1) - len1) .lt. lenadd) then 
     allocate(oneline2(len1))
     oneline2=oneline1
     deallocate(oneline1)
     allocate(oneline1(len1+buf_size))
     oneline1(1:len1)=oneline2(1:len1)
  deallocate(oneline2)
    endif  
    do j=1, lenadd
  oneline1(len1+j)=fidochar2(j:j)
 enddo
 
   len1=lenadd+len1
!  deallocate(oneline2)
enddo
output_len=len1
return
endif

if(input_len .eq. 1) then
   sum=sum+1
   point_cur=input_array(1)
   if(point_cur .lt. 0) then
     ES_bit=7
   else
     ES_bit=7
   endif
   if (typeflag .gt. 0) then
    i_value=int(point_cur)
    call intlen(i_value, point_bit)
    point_bit=point_bit+typeflag+1
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform,"('(f',I0,'.',I0,')' )") point_bit, typeflag
    write(fidochar2,fidoform) point_cur
   elseif(typeflag .eq. 0) then
    i_value=nint(point_cur)
    call intlen(i_value, point_bit)
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidochar2,'(I0)') i_value
   else
    point_bit=abs(typeflag)+Es_bit
 
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform, "('(ES',I0,'.',I0,')' )") point_bit,typeflag_in 
    write(fidochar2,fidoform) point_cur 
   endif
   fidochar2=adjustl(fidochar2)
   lenadd=point_bit+1
   deallocate(oneline1)
   allocate(oneline1(lenadd))
   do i=1,lenadd
     oneline1(i)=fidochar2(i:i)
   enddo
   len1=lenadd+len1
endif


do m=2, input_len          !do_label_2:   point to point sweep 
   if(skip .ne. 0) then    !if_label skip
     skip=skip-1
   else  
      
   point_pre=input_array(m-1)  
   point_cur=input_array(m)

   if(point_pre .lt. 0) then
     ES_bit=7
   else
     ES_bit=6
   endif

! check whether current  num is the same as the previous one
   if(point_cur .eq. point_pre  )   i=i+1     !if same, add 1
   if(point_cur .ne. point_pre .or. m .eq. input_len) then !if_label_1  

!   allocate(oneline2(len1))
!   oneline2=oneline1
!   deallocate(oneline1)
      
   if (typeflag .gt. 0) then
    i_value=int(point_pre)
    call intlen(i_value, point_bit)
    point_bit=point_bit+typeflag+1
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform,"('(f',I0,'.',I0,')' )") point_bit, typeflag
       write(fidochar2,fidoform) point_pre
   elseif(typeflag .eq. 0) then
       i_value=nint(point_pre)
       call intlen(i_value, point_bit)
       if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
       write(fidochar2,'(I0)') i_value
   
   else
    point_bit=abs(typeflag)+ES_bit
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform, "('(ES',I0,'.',I0,')' )") point_bit,typeflag_in 
           write(fidochar2,fidoform) point_pre       
   
   endif
   sum=sum+i
   
   
   if(i .eq. 1) then               !if_label_3: see whether it only repeats 1 time 
     fidochar2=adjustl(fidochar2)
     lenadd=1+point_bit            !plus a space  

!       allocate(oneline1(len1+lenadd))
!     do j=1,len1
!    oneline1(j)=oneline2(j)
!     enddo
    if( (size(oneline1) - len1) .lt. lenadd) then 
     allocate(oneline2(len1))
     oneline2=oneline1
     deallocate(oneline1)
     allocate(oneline1(len1+buf_size))
     oneline1(1:len1)=oneline2(1:len1)
  deallocate(oneline2)
    endif
        do j=1, lenadd
   oneline1(len1+j)=fidochar2(j:j)
  enddo
   else 
     call fido_generating(i, fidochar2,'R', fidoset,lenset)  
        if(lenset .gt. 59) stop 'FIDO(len>59)'
  lenadd=lenset+1
!  allocate(oneline1(len1+lenadd))
!       do j=1,len1
!    oneline1(j)=oneline2(j)
!     enddo
   if( (size(oneline1) - len1) .lt. lenadd) then 
     allocate(oneline2(len1))
     oneline2=oneline1
     deallocate(oneline1)
     allocate(oneline1(len1+buf_size))
     oneline1(1:len1)=oneline2(1:len1)
     deallocate(oneline2)
   endif
   do j=1,lenadd
    oneline1(len1+j)=fidoset(j:j)
   enddo

   endif !if_label_3 end
                           
      len1=len1+lenadd
!   deallocate(oneline2)
      i=1     
         
   endif   !if_label_1: end
 
!haddle with the last pixel in one z layer if it is different with the previous one
   if(point_cur .ne. point_pre .and. m .eq. input_len) then !if_label 4
      
   sum=sum+1
   if (typeflag .gt. 0) then
    i_value=int(point_cur)
    call intlen(i_value, point_bit)
    point_bit=point_bit+typeflag+1
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform,"('(f',I0,'.',I0,')' )") point_bit, typeflag
    write(fidochar2,fidoform) point_cur
   elseif(typeflag .eq. 0) then
       i_value=nint(point_cur)
       call intlen(i_value, point_bit)
       if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
       write(fidochar2,'(I0)') i_value
   else
    point_bit=abs(typeflag)+ES_bit
    if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
    write(fidoform, "('(ES',I0,'.',I0,')' )") point_bit,typeflag_in 
    write(fidochar2,fidoform ) point_cur  
   
   endif
   fidochar2=adjustl(fidochar2)  
   lenadd=point_bit+1

    if( (size(oneline1) - len1) .lt. lenadd) then 
     allocate(oneline2(len1))
     oneline2=oneline1
     deallocate(oneline1)
     allocate(oneline1(len1+buf_size))
     oneline1(1:len1)=oneline2(1:len1)
  deallocate(oneline2)
    endif
     ! allocate(oneline2(len1))
 ! do j=1, len1 
 !   oneline2(j)=oneline1(j)
 ! enddo
 ! deallocate(oneline1)
 !  allocate(oneline1(lenadd+len1))
    !  do j=1, len1 
 !    oneline1(j)=oneline2(j)
 !  enddo
      do j=1,lenadd
     oneline1(len1+j)=fidochar2(j:j)
      enddo
   
   len1=len1+lenadd
!   deallocate(oneline2)
   
 endif         !if_label_4 end
    
 ! handeling Q    
 
 if(abs(qfidoflag) .gt. 0 .and. point_cur .ne. point_pre) then !if_label_5     
      
   len_before=m-1
   len_after=input_len-m
    qbingo=0
      if (qfidoflag_in .eq. 1) qfidoflag_in=input_len     
   scan_scope=min(qfidoflag_in, min(len_before,len_after+1))
   if(scan_scope .ge. 2) then     !at least go back 2
    allocate(qposition(qbingo),stat=status)
    if( status .ne. 0) then
          deallocate(qposition)
          allocate(qposition(qbingo))
       endif
    
    scan_begin=m-scan_scope
    if(qfidoflag .lt. 0) then
      scan_skip=scan_scope
    else
      scan_skip=1 
       endif
    
    do scan=scan_begin,m-2, scan_skip     !label_do: back scan
    
   flag=real_compare(input_array(scan:m-1),input_array(m:2*m-1-scan), m-scan)
              
      if(flag .eq. 1) then     
     allocate(qposition_temp(qbingo))
           qposition_temp=qposition
     deallocate(qposition)
     qbingo=qbingo+1
     allocate(qposition(qbingo))
     do pos=1, qbingo-1
       qposition(pos)=qposition_temp(pos)
     enddo  
     qposition(qbingo)=scan         
     deallocate(qposition_temp)
   endif
   
       enddo                            !lable_do: back scan
    if( qbingo .ne. 0) then          !if_label:Q
      q_repeat=1
   qlength=m-qposition(1)
         flag=1
   afterq=m+qlength
   do while (afterq+qlength-1 .le. input_len .and. flag .ne. 0)
    flag=real_compare(input_array(qposition(1):m-1),input_array(&
         afterq:afterq+qlength-1),qlength)
    if(flag .eq. 1) then 
      q_repeat=q_repeat+1  
      afterq=afterq+qlength
    endif
   enddo

!        allocate(oneline2(len1))
!      oneline2=oneline1
!      deallocate(oneline1)
   call intlen(qlength, point_bit)
   if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
   write(fidochar2,'(I0)') qlength
   call fido_generating(q_repeat, fidochar2,'Q', fidoset,lenset)  
   if(lenset .gt. 59) stop 'FIDO(len>59)'
   lenadd=lenset+1
   if( (size(oneline1) - len1) .lt. lenadd) then 
     allocate(oneline2(len1))
     oneline2=oneline1
     deallocate(oneline1)
     allocate(oneline1(len1+buf_size))
     oneline1(1:len1)=oneline2(1:len1)
     deallocate(oneline2)
   endif
   do j=1,lenadd
     oneline1(len1+j)=fidoset(j:j)
   enddo 
   sum=sum+qlength*q_repeat
   len1=len1+lenadd
!  deallocate(oneline2)
   deallocate(qposition)
   if(m+qlength*q_repeat-1 .lt. input_len-1) then   !if_label: qend
     skip=qlength*q_repeat
   elseif(m+qlength*q_repeat-1 .eq. input_len-1) then
          
   if (typeflag .gt. 0) then
      i_value=int(input_array(input_len))
      call intlen(i_value, point_bit)
      point_bit=point_bit+typeflag+1
      if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
      write(fidoform,"('(f',I0,'.',I0,')' )") point_bit, typeflag
      write(fidochar2,fidoform ) input_array(input_len)
   elseif(typeflag .eq. 0) then
      i_value=nint(input_array(input_len))
      call intlen(i_value, point_bit)
      if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
        write(fidochar2,'(I0)') i_value
      else
        if(input_array(input_len) .lt. 0) then
          point_bit=abs(typeflag)+7
        else
          point_bit=abs(typeflag)+6
        endif
   if (point_bit .gt. 19) stop 'too many digits(>19) for the repeating number in fido'
      write(fidoform, "('(ES',I0,'.',I0,')' )") point_bit,typeflag_in 
      write(fidochar2,fidoform ) input_array(input_len)     
   endif
   fidochar2=adjustl(fidochar2)                
   lenadd=point_bit+1
   if( (size(oneline1) - len1) .lt. lenadd) then 
     allocate(oneline2(len1))
     oneline2=oneline1
     deallocate(oneline1)
     allocate(oneline1(len1+buf_size))
     oneline1(1:len1)=oneline2(1:len1)
     deallocate(oneline2)
   endif

!    allocate(oneline2(len1))
!    do j=1, len1 
!      oneline2(j)=oneline1(j)
!    enddo
!       deallocate(oneline1)
!       allocate(oneline1(lenadd+len1))
!          do j=1, len1 
!         oneline1(j)=oneline2(j)
!       enddo
    do j=1,lenadd
      oneline1(len1+j)=fidochar2(j:j)
    enddo
   
    len1=len1+lenadd
!       deallocate(oneline2) 
    sum=sum+1
    goto 100
   else
    goto 100
  
   endif    !if_label: qend
      
    endif  !if_lable: Q    
      
   endif                         !at lease go back 2
   
 endif                                                !if_label_5 end 

endif  !if_label skip
120   enddo          ! end m               

if(input_len .eq. 0) then
  write(warning_message,*) 'fido routine called on a null input array'
  call TrapError(-5)
endif

100 output_len=len1

if(num_fido_char .gt. max_fido_char) max_fido_char=num_fido_char

if(sum .ne. input_len) then
  write(*,*) ' error in converting fido'
endif

return
end subroutine


!Subroutine FidoIntArray
!generating fido array from a integer input array
! be careful: before calling this sub, the input array should be loaded into 
! input array(input_array), after this call, data in fido output array(oneline1,global) 
! should be extracted, after that, deallocating oneline1 is recommended(not necessary).
! or during  next call, the old data in oneline1 will be lost. 
! qfidoflag: handling FIDO character Q
! qfidoflag=0, only generate R
! qfidoflag>=2, Q may be generated in the output array. the number is the scan scope
! e.g. if qfidoflag=100, that will force the sub only scans at most the last 100 numbers 
! be carefull, if input_array and qfidoflag are both big, it may take some time to
! finish scan. but if qfidoflag is too small, it may be not enough to generate Q 
! qfidoflag <= -2 (negative) that will force this sub only to scan for abs(qfidoflag) numbers
! e.g. if qfidoflag=100 , the input array has 600 numbers, when the sub processes No.102,it 
! may generate FIDO like 2Q100, 4Q99 ... (no fido like 4Q101)
! if qfidoflag=-100, it only could generate fido like 3Q100...(no fido like 2Q99 or 2Q101) 
! qfidoflag =1  : this will force the sub to scan from the first possible leftmost number in
! the input array 
! qfidoflag=-1 : like qfidoflag negative, only scan once on the first possible leftmost number
! the ability to generate Q is limited. R has higher priority than Q
! only R and Q in FIDO character set can be generated now, and Q generating is limited 
 
subroutine FidoIntArray(input_array,input_len,qfidoflag,rfidoflag)
use fido
!optional :: rfidoflag
integer input_len
integer input_array(input_len)
integer qfidoflag, rfidoflag

allocate(fido_input(input_len))
call ConvertIntToReal(input_array,fido_input,input_len)
!if(present(rfidoflag) ) then
call FidoRealArray(fido_input,input_len,0,qfidoflag, rfidoflag)
!else
!call FidoRealArray(fido_input,input_len,0,qfidoflag)
!endif
deallocate(fido_input)

return
end subroutine

!Subroutine ConvertIntToReal
!convert an integer array to a real(4) array
subroutine ConvertIntToReal(int_array, real_array,input_len)
integer input_len
integer int_array(input_len)
real    real_array(input_len)

integer i

do i=1,input_len
 real_array(i)=real(int_array(i))
end do

return
end subroutine
