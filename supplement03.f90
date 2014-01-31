!write CM mat number to file in HDF5 format

module MyHDF5
  use HDF5
    
  integer ::err_h5  ! Error flag
  integer(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id       ! Dataset identifier
  integer(HID_T) :: dspace_id     ! Dataspace identifier

!  integer(HID_T), allocatable :: dset_id(:,:,:)       ! Dataset identifier
!  integer(HID_T), allocatable :: dspace_id(:,:,:)     ! Dataspace identifier
  integer(HSIZE_T), dimension(3) :: ddims
end module 

! write mat number to a h5 file
subroutine WriteDotH5
  use MyHDF5
  use files
  use ErrControl
  use paraset1
  use paraset2
  use paraset3
  use paraset4
! Initialize h5 FORTRAN interface.
  
  implicit none
  integer cmi, cmj, cmk
  integer i, j, k
  integer rank, cm_num
  
  character (LEN=20) :: dsetname=''
  
  cm_num=0  
  rank=3
  
  call h5open_f(err_h5)
  cur_filename=trim(prbname)//".h5"

! Create a new file using default properties.  
  call h5fcreate_f(cur_filename, H5F_ACC_TRUNC_F, file_id, err_h5)
  
!  allocate (dspace_id(num_cmesh(1), num_cmesh(2), num_cmesh(3))
!  allocate (dset_id((num_cmesh(1), num_cmesh(2), num_cmesh(3))
  
  do cmk=1, num_cmesh(3) 
    do cmj=1,num_cmesh(2)
      do cmi=1, num_cmesh(1)
         cm_num=cm_num+1
         ! Create the dataspace.
         ddims=zlevel(cmk)%cm_zlev(cmi,cmj)%fine
         call h5screate_simple_f(rank, ddims, dspace_id, err_h5)
         ! Create the dataset with default properties.
         write (dsetname,"('cm', I3.3)") cm_num
         call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, dset_id, err_h5)
         call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, zlevel(cmk)%cm_zlev(cmi,cmj)%mat_matrix, ddims, err_h5)
         ! End access to the dataset and release resources used by it.
         call h5dclose_f(dset_id, err_h5)
         ! Terminate access to the data space.
         call h5sclose_f(dspace_id, err_h5)

      enddo
    enddo
  enddo
  
  
  ! Close the file.
  call h5fclose_f(file_id, err_h5)

  ! Close FORTRAN interface.
  call h5close_f(err_h5)

end subroutine

    
    
