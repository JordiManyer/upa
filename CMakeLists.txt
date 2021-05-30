cmake_minimum_required(VERSION 3.11)

# Set project name
project(UPA)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set installation directories
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR})

set(UPA_LIBS_INSTALL ${CMAKE_INSTALL_PREFIX}/libs)
set(UPA_INCLUDE_INSTALL ${CMAKE_INSTALL_PREFIX}/include)
set(UPA_BIN_INSTALL ${CMAKE_INSTALL_PREFIX}/bin)

# Setup placeholders for project libraries
set(UPA_LIBRARIES CACHE INTERNAL "UPA_LIBRARIES")
set(EXTERNAL_LIBRARIES CACHE INTERNAL "EXTERNAL_LIBRARIES")

# External libraries
set(CMAKE_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/build/myCmakeModules)

message(STATUS ">> Finding MPI library: " )
find_package(MPI REQUIRED)
include_directories(SYSTEM "${MPI_CXX_INCLUDE_DIRS}")
set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${MPI_CXX_LIBRARIES} CACHE INTERNAL "EXTERNAL_LIBRARIES")
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "-std=c++11")

message(STATUS ">> Finding METIS library: " )
find_package(METIS REQUIRED)
include_directories(SYSTEM "${METIS_INCLUDE_DIRS}")
set(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${METIS_LIBRARIES} CACHE INTERNAL "EXTERNAL_LIBRARIES")

# Includes
file(GLOB UPA_INC_DIRS ${CMAKE_CURRENT_LIST_DIR}/src/*)
include_directories(${UPA_INC_DIRS})

# Build UPA
add_subdirectory(src/Integration)     # Integration
add_subdirectory(src/IO)              # Input/Output routines
add_subdirectory(src/Math)            # Basic Math
add_subdirectory(src/Preconditioners) # Preconditioners and Iterators
add_subdirectory(src/Solvers)         # Sparse Solvers
add_subdirectory(src/DataStructures)  # Data Structures (Sparse matrices, graphs, maps, ...)
add_subdirectory(src/Mesh)  # Meshes

add_subdirectory(src/Tests)           # Executables




