cmake_minimum_required(VERSION 3.10)

project(forward_kinematics_task)

set(BUILD_MODE RELEASE)
set(PLATFORM PLATFORM_DESKTOP)

add_executable(forward_kinematics src/controller.cpp)

add_custom_target(
    build_raylib ALL
    COMMAND mkdir raylib && cd raylib && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=RELEASE ${CMAKE_SOURCE_DIR}/3rdParty/raylib && make || echo "raylib is already built"
)
add_dependencies(forward_kinematics build_raylib)

add_definitions(-D_DEFAULT_SOURCE -DRAYLIB_BUILD_MODE=${BUILD_MODE} -D${PLATFORM} -DNDEBUG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -Wall -Wextra -Wpedantic -O3 -march=native -ffast-math")

target_link_directories(forward_kinematics PUBLIC ${CMAKE_BINARY_DIR}/raylib/raylib)

target_link_libraries(forward_kinematics PUBLIC m pthread rt dl raylib)

target_include_directories(forward_kinematics PUBLIC 3rdParty/raylib/src 3rdParty/raygui/src)

