cmake_minimum_required(VERSION 3.6)

PROJECT(proj3)
SET(VTK_DIR /home/himani/Desktop/SciVIs/VTK7_Build)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(proj3 proj3)

#SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(proj3 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj3 vtkHybrid)
endif()


