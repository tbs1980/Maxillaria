cmake_minimum_required(VERSION 2.6)

#set compiler flags
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -pedantic -Wfatal-errors -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -pedantic -Wfatal-errors -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -pedantic -Wall -Wextra -Wfatal-errors -O3")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

find_library(M_LIB m)


add_executable(unit_gauss.exe unit_gauss.c priors.c)
target_link_libraries(unit_gauss.exe maxillaria ${M_LIB})



