# FindFFTW.cmake -- locate FFTW3 (double precision) headers and library.
# Sets FFTW_FOUND, FFTW_INCLUDE_DIRS, FFTW_LIBRARIES. Honours FFTW_ROOT.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_FFTW QUIET fftw3)
endif()

find_path(FFTW_INCLUDE_DIR fftw3.f03
  HINTS ${FFTW_ROOT} $ENV{FFTW_ROOT} ${PC_FFTW_INCLUDEDIR} ${PC_FFTW_INCLUDE_DIRS}
        /opt/homebrew/opt/fftw /usr/local/opt/fftw
  PATH_SUFFIXES include)
find_library(FFTW_LIBRARY NAMES fftw3
  HINTS ${FFTW_ROOT} $ENV{FFTW_ROOT} ${PC_FFTW_LIBDIR} ${PC_FFTW_LIBRARY_DIRS}
        /opt/homebrew/opt/fftw /usr/local/opt/fftw
  PATH_SUFFIXES lib lib64)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDE_DIR)
if(FFTW_FOUND)
  set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
  set(FFTW_LIBRARIES ${FFTW_LIBRARY})
endif()
mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)
