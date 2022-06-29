cmake_minimum_required(VERSION 3.14)
set(CMAKE_TOOLCHAIN_FILE "C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake")
include( "C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake")
if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET} AND NOT DEFINED VCPKG_TARGET_TRIPLET)
    set(VCPKG_TARGET_TRIPLET "$ENV{VCPKG_DEFAULT_TRIPLET}" CACHE STRING "" FORCE)
endif()

project(GBTStraightPipeLinear)

set(CMAKE_CXX_STANDARD 20)


#=================== INCLUSION OF Project Files ====================#
set(FORMS_DIR "${CMAKE_SOURCE_DIR}/forms")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")

include_directories(${FORMS_DIR})
include_directories(${INCLUDE_DIR})
include_directories(${SOURCE_DIR})

file(GLOB_RECURSE SOURCES
    "${FORMS_DIR}/*.qrc"
    "${INCLUDE_DIR}/*.h"
    "${SOURCE_DIR}/*.cpp"
)
#=================== SETUP EXECTUABLE ====================#

# Add the executable
add_executable(GBTStraightPipeLinear EXCLUDE_FROM_ALL ${SOURCES})

# Add the target includes for MY_PROJECT 
target_include_directories(GBTStraightPipeLinear PRIVATE ${FORMS_DIR})
target_include_directories(GBTStraightPipeLinear PRIVATE ${INCLUDE_DIR})
target_include_directories(GBTStraightPipeLinear PRIVATE ${SOURCE_DIR})


#===================== LINKING LIBRARIES =======================#

find_package(matplot++ REQUIRED)
target_link_libraries(GBTStraightPipeLinear  PRIVATE Matplot++::cimg std::filesystem Matplot++::matplot Matplot++::nodesoup)

