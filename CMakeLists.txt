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

# Libraries
set(UPA_LIBRARIES CACHE INTERNAL "UPA_LIBRARIES")

# Includes
file(GLOB UPA_INC_DIRS ${CMAKE_CURRENT_LIST_DIR}/src/*)
include_directories(${UPA_INC_DIRS})

# Build UPA
add_subdirectory(src/Integration)     # Integration
add_subdirectory(src/IO)              # Input/Output routines
add_subdirectory(src/Math)            # Basic Math
add_subdirectory(src/Preconditioners) # Preconditioners and Iterators
add_subdirectory(src/Solvers)         # Sparse Solvers
add_subdirectory(src/Sparse)          # Sparse Matrices
add_subdirectory(src/StructuredMesh)  # Meshes

add_subdirectory(src/Tests)           # Executables




