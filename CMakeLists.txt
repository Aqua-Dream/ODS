cmake_minimum_required(VERSION 3.6)
project(ODS LANGUAGES CXX)
set(ODS_SOURCE_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCE_FILES main.cpp Feistel.cpp memory.cpp sampling.cpp global.cpp ODS.cpp argsolver.cpp)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
	message("OPENMP FOUND")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    if(NUM_THREADS)
        add_definitions(-DNUM_THREADS=${NUM_THREADS})
    else()
        add_definitions(-DNUM_THREADS=4)
    endif()
endif()

# Your-external "mylib", add GLOBAL if the imported library is located in directories above the current.
add_library( blake3 SHARED IMPORTED)
# You can define two import-locations: one for debug and one for release.
set_target_properties(blake3 PROPERTIES IMPORTED_LOCATION ${ODS_SOURCE_ROOT}/libblake3.so)

add_executable(main ${SOURCE_FILES})

target_link_libraries(main blake3)
