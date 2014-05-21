# - Try to find dislin
# Once done this will define
#  DISLIN_FOUND - System has dislin
#  DISLIN_INCLUDE_DIRS - The dislin include directories
#  DISLIN_LIBRARIES - The libraries needed to use dislin
#  DISLIN_DEFINITIONS - Compiler switches required for using dislin
 
find_package(PkgConfig)
pkg_check_modules(PC_LIBDISLIN QUIET libdislin)
set(DISLIN_DEFINITIONS ${PC_LIBDISLIN_CFLAGS_OTHER})
 
find_path(DISLIN_INCLUDE_DIR dislin/dislin.h
          HINTS ${PC_DISLIN_INCLUDEDIR} ${PC_DISLIN_INCLUDE_DIRS} $ENV{DISLIN_INCDIR}
          PATH_SUFFIXES dislin )
 
find_library(DISLIN_LIBRARY NAMES dislin
             HINTS ${PC_DISLIN_LIBDIR} ${PC_DISLIN_LIBRARY_DIRS} $ENV{DISLIN_LIBDIR} 
             PATH_SUFFIXES dislin lib )
 
set(DISLIN_LIBRARIES ${DISLIN_LIBRARY} )
set(DISLIN_INCLUDE_DIRS ${DISLIN_INCLUDE_DIR} )
 
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set DISLIN_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(dislin  DEFAULT_MSG
                                  DISLIN_LIBRARY DISLIN_INCLUDE_DIR)
 
mark_as_advanced(DISLIN_INCLUDE_DIR DISLIN_LIBRARY )
