#!For finding Boost specify only root directory-DBOOST_ROOT=...Software/boost_1_67_0/
#!For finding CERN Root libraries execution of 'source .../root_v6.28.06/build/bin/thisroot.sh'
#is required either at the start of each bash session or at system startup

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16)
project(Data_processing)
cmake_policy(SET CMP0144 NEW)

find_package(Boost 1.33 REQUIRED)
find_package(ROOT 6.04 REQUIRED)
# Add pthread library:
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

#----------------------------------------------------------------------------
# Setup include directories for this project
#
#include(${Boost_USE_FILE}) does not work
#include(${Geant4_USE_FILE})
include_directories(${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file (GLOB_RECURSE sources1 ${PROJECT_SOURCE_DIR}/source/*.cpp)
list (APPEND sources ${sources1})
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB_RECURSE headers1 ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE headers2 ${PROJECT_SOURCE_DIR}/include/*.hpp)
list (APPEND headers ${headers1})
list (APPEND headers ${headers2})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

#----------------------------------------------------------------------------
# Set -rpath for binary. Allows it to find libraries in non-stardard locations
# (e.g. not in /usr/) without using LD_LIBRARY_PATH
#
#set(CMAKE_SKIP_BUILD_RPATH FALSE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
#list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#if("${isSystemDir}" STREQUAL "-1")
#    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#endif("${isSystemDir}" STREQUAL "-1")

#----------------------------------------------------------------------------
#use the flags used in ROOT's compilation
set(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS})
#use the flags used in ROOT's compilation
set(CMAKE_SHARED_LINKER_FLAGS ${ROOT_EXE_LINKER_FLAGS})
message(STATUS "CXX compilation flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "CXX linking flags: ${CMAKE_SHARED_LINKER_FLAGS}")

#----------------------------------------------------------------------------
# Remove some of the warnings for gcc as
# boost has quite a lot of these particular ones
#
if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # additional warnings
    add_compile_options(-Wno-pedantic -Wno-shadow)
endif()

#----------------------------------------------------------------------------
# Add the executable, and link it to the libraries
#
add_executable(${PROJECT_NAME} ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} PUBLIC Threads::Threads ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
message(STATUS "ROOT libraries: ${ROOT_LIBS}")

#----------------------------------------------------------------------------
# No installaion is required
#
#install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
