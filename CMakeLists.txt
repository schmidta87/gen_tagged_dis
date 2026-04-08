# CMakeLists for the e2a analysis repository

# Set up
cmake_minimum_required (VERSION 3.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (tagged_dis_eic)

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

#Add a subdirectory for cross section library
add_subdirectory(Cross_Sections)


add_executable(gen_tagged_dis gen_tagged_dis.cpp)
target_link_libraries(gen_tagged_dis ${ROOT_LIBRARIES} )
