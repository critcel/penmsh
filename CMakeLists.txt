cmake_minimum_required (VERSION 2.6)
project (penmsh)
enable_language (Fortran)

# FFLAGS depend on the compiler
# message ("FC = " $ENV{FC})
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    set (CMAKE_Fortran_FLAGS "-O3 -ffree-form -ffree-line-length-none")
    message ("Fortran compiler options: " ${CMAKE_Fortran_FLAGS})

elseif (Fortran_COMPILER_NAME MATCHES "nagfor.*")
    # nagfor
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    set (CMAKE_Fortran_FLAGS "-O3")
    message ("Fortran compiler options: " ${CMAKE_Fortran_FLAGS})

endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# set(extern_LIB /usr/local/dislin)
# set(extern_LIB ${CMAKE_CURRENT_SOURCE_DIR})
# set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
message ("CMAKE_MODULE_PATH = " ${CMAKE_MODULE_PATH})
find_package (DISLIN REQUIRED)
include_directories(${DISLIN_INCLUDE_DIRS})
set(LIBS ${LIBS} ${DISLIN_LIBRARIES})
set(LIBS ${LIBS} user32.lib gdi32.lib)

add_executable (penmsh dislin.f90
                inpred.f90 subcode0.f90 subcode1.f90 subcode2.f90
                subcode3.f90 subcode4.f90 subcode5.f90 subcode6.f90 
                subcode7.f90 subcode8.f90 mytecio.f90  
                supplement01.f90 supplement02-plot.f90)

target_link_libraries (penmsh ${LIBS})
