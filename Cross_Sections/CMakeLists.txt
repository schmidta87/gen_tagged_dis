# CMakeLists for the cross section libary

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add Cross_Section library
add_library(Cross_Sections Cross_Sections.cpp)
target_link_libraries(Cross_Sections ${ROOT_LIBRARIES})
