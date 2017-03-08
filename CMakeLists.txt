cmake_minimum_required(VERSION 2.6)
project(qphot C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(phypp)

include_directories(${PHYPP_INCLUDE_DIRS})

add_executable(qphot qphot.cpp)
target_link_libraries(qphot ${PHYPP_LIBRARIES})
install(TARGETS qphot DESTINATION bin)

add_executable(qphot-psfsub qphot-psfsub.cpp)
target_link_libraries(qphot-psfsub ${PHYPP_LIBRARIES})
install(TARGETS qphot-psfsub DESTINATION bin)
