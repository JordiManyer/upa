# CMAKE for the Solvers library

# Get source files
file(GLOB SOURCES_H   *.h)
file(GLOB SOURCES_CPP *.cpp)
set(SOURCE_FILES ${SOURCES_H} ${SOURCES_CPP})
set(LIBNAME upa_solvers)

# Build library
add_library(${LIBNAME} ${SOURCE_FILES})
target_link_libraries(${LIBNAME} ${UPA_LIBRARIES})

install(TARGETS ${LIBNAME} DESTINATION ${UPA_LIBS_INSTALL})
install(FILES ${SOURCES_H} DESTINATION ${UPA_INCLUDE_INSTALL})

# Store library
set(UPA_LIBRARIES ${UPA_LIBRARIES} ${LIBNAME} CACHE INTERNAL "UPA_LIBRARIES")



