cmake_minimum_required(VERSION 3.8)
project(flow-simulator)

# auxiliary cmake files should be placed in /cmake/Modules/ .
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default C++ flags
set(COMMON "-fdiagnostics-color")
set(WARNINGS "-Wall -Wextra -Wshadow -Wsign-conversion -Wsign-promo -Wpedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON} ${WARNINGS}")

# Task-specific C++ flags
set(PROFILE "-lprofiler")
set(DEBUG "-ggdb3")
set(COVERAGE "--coverage")
set(MACHINE "-march=native")
set(OPENMP "-fopenmp")

# Set CMake build types
# Under development, one may add flag -D_VALIDATION
set(CMAKE_CXX_FLAGS_DEBUG   "-Og ${DEBUG} -Wno-unknown-pragmas")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${MACHINE} ${OPENMP} -DNDEBUG -Werror")
set(CMAKE_CXX_FLAGS_PROFILE "-Og ${MACHINE} ${DEBUG} ${OPENMP} ${PROFILE}")

# output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# additional cmake options
include(CMakeToolsHelpers OPTIONAL)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

# additional include directories (-I in g++).
include_directories(./)

##### Dependencies
# required packages
find_package(Armadillo REQUIRED)
find_package(Glog REQUIRED)
find_package(HDF5 REQUIRED)
find_package(RapidJSON REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(TCLAP REQUIRED)
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
include_directories(SYSTEM ${TCLAP_INCLUDE_DIRS})
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
include_directories(SYSTEM ${RAPIDJSON_INCLUDE_DIRS})
list(APPEND LIBS ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES} ${SUNDIALS_LIBRARIES})


##### sources
# production
file(GLOB         MAIN    "src/main.*")
file(GLOB_RECURSE HEADERS "src/*/*.h")
file(GLOB_RECURSE SOURCES "src/*/*.cc")

# test
set(MAIN_TEST                   "test/main_test.cc")
set(MAIN_STRESS_TEST            "test/main_stress_test.cc")
file(GLOB_RECURSE HEADERS_TEST  "test/src/*.h")

##### executables
# production
add_executable(flow-simulator.x ${MAIN} ${HEADERS} ${SOURCES})
target_link_libraries(flow-simulator.x ${LIBS} ${GLOG_LIBRARIES})

# test
find_package(GTest)
find_package(Threads)
IF (Threads_FOUND AND GTest_FOUND)
    include_directories(SYSTEM ${GTest_INCLUDE_DIR})
    include_directories(SYSTEM ${Threads_INCLUDE_DIR})
    add_executable(flow-simulator-test.x ${MAIN_TEST} ${SOURCES} ${HEADERS} ${HEADERS_TEST})
    target_link_libraries(flow-simulator-test.x ${LIBS} ${GLOG_LIBRARIES} ${GTEST_LIBRARIES} ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    set_target_properties(flow-simulator-test.x PROPERTIES  COMPILE_FLAGS "${COVERAGE}"
                                                       LINK_FLAGS    "${COVERAGE}")
    add_executable(flow-simulator-stress-test.x ${MAIN_STRESS_TEST} ${SOURCES} ${HEADERS} ${HEADERS_TEST})
    target_link_libraries(flow-simulator-stress-test.x ${LIBS} ${GLOG_LIBRARIES} ${GTEST_LIBRARIES} ${THREADS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF ()

# profile
IF (CMAKE_BUILD_TYPE MATCHES Profile)
    find_package(Gperftools REQUIRED)
    add_executable(flow-simulator-profile.x ${MAIN} ${HEADERS} ${SOURCES})
    target_link_libraries(flow-simulator-profile.x ${LIBS} ${GLOG_LIBRARIES} ${LD_LIBRARY})
ENDIF()

##### auxiliary make directives
# cpplint
add_custom_target(cpplint
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND cpplint --recursive --quiet --linelength=100 --extensions=h,cc --filter=-runtime/references,-whitespace/empty_loop_body src/ test/
  USES_TERMINAL
)

# docs
add_custom_target(docs
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND doxygen
  COMMAND gitstats . docs/stats
  USES_TERMINAL
)

# coverage
set(COVERAGE_INFO_FILE "docs/coverage/coverage.info")
add_custom_target(coverage
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  DEPENDS flow-simulator-test.x
  COMMAND flow-simulator-test.x &>/dev/null
  COMMAND lcov --quiet --no-external --capture --directory . --output-file ${COVERAGE_INFO_FILE}
  COMMAND lcov --quiet --remove ${COVERAGE_INFO_FILE} 'test/*' -o ${COVERAGE_INFO_FILE}
  COMMAND genhtml ${COVERAGE_INFO_FILE} --quiet --show-details --legend --output-directory docs/coverage/
  COMMAND lcov --list ${COVERAGE_INFO_FILE}
  USES_TERMINAL
)

# clean-docs
add_custom_target(clean-docs
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/clean-docs.cmake
)

# clean-results
add_custom_target(clean-results
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/clean-results.cmake
)

# clean-debug
add_custom_target(clean-debug
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/clean-debug.cmake
)

# clean-cmake
add_custom_target(clean-cmake
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/clean-cmake.cmake
)

# clean-all
add_custom_target(clean-all
  COMMAND ${CMAKE_BUILD_TOOL} clean
  COMMAND ${CMAKE_BUILD_TOOL} clean-docs
  COMMAND ${CMAKE_BUILD_TOOL} clean-results
  COMMAND ${CMAKE_BUILD_TOOL} clean-debug
  COMMAND ${CMAKE_BUILD_TOOL} clean-cmake
)
