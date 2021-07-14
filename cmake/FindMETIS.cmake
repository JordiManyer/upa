
# Defines:
#     METIS_FOUND        - system has METIS
#     METIS_INCLUDE_DIRS - include directories for METIS
#     METIS_LIBRARIES    - libraries for METIS
#=============================================================================

find_path(METIS_INCLUDE_DIR
        NAMES metis.h
        PATHS ${CMAKE_INCLUDE_PATH} ${METIS_INCLUDE_DIR} ENV METIS_INCLUDE_DIR ${METIS_DIR} ENV METIS_DIR
        PATH_SUFFIXES include
        DOC "Directory where the METIS header files are located"
        )

find_library(METIS_LIBRARY
        NAMES metis metis${METIS_LIB_SUFFIX} libmetis
        HINTS ${CMAKE_LIBRARY_PATH} ${METIS_LIB_DIR} ENV METIS_LIB_DIR ${METIS_DIR} ENV METIS_DIR
        PATH_SUFFIXES lib lib64
        DOC "Directory where the METIS library is located"
        )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS DEFAULT_MSG METIS_LIBRARY METIS_INCLUDE_DIR)

if(METIS_FOUND)
    set(METIS_LIBRARIES ${METIS_LIBRARY})
    set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
endif()

mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY)