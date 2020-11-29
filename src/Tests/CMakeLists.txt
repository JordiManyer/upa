# CMake file for executables

#Get sources for tests
file(GLOB PROGRAMS_CPP *.cpp)


#Configure tests
message(STATUS "Configuring Tests")
foreach(program ${PROGRAMS_CPP})

    get_filename_component(binary ${program} NAME)
    string(REPLACE ".cpp" "" binary ${binary})
    message(STATUS "  ${binary}")

    # Build
    add_executable(${binary} ${program})
    target_link_libraries(${binary} ${UPA_LIBRARIES})
    set_target_properties(${binary} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${UPA_BIN_INSTALL})

endforeach(program ${PROGRAMS_CPP})