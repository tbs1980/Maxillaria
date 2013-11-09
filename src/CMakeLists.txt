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


add_library(maxillaria nested_sampler.c live_point.c mt19937.c nested_diagnostic.c utils.c post_equal_weights.c post_stats.c explore.c)
