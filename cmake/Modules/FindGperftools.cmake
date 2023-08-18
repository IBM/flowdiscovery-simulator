# FindGperftools.cmake exports:
#
# LD_LIBRARY: library used to perform profiling
#

include(FindPackageHandleStandardArgs)

find_library(LD_LIBRARY
  NAMES profiler
  PATHS
      /usr/local/lib
      /usr/lib
)

find_package_handle_standard_args(Gperftools DEFAULT_MSG LD_LIBRARY)
