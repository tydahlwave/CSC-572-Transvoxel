cmake_minimum_required(VERSION 3.0)
project(LOGLEngine)

option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
add_subdirectory(LOGLEngine/Vendor/glfw)

option(ASSIMP_BUILD_ASSIMP_TOOLS OFF)
option(ASSIMP_BUILD_SAMPLES OFF)
option(ASSIMP_BUILD_TESTS OFF)
add_subdirectory(LOGLEngine/Vendor/assimp)

option(BUILD_BULLET2_DEMOS OFF)
option(BUILD_CPU_DEMOS OFF)
option(BUILD_EXTRAS OFF)
option(BUILD_OPENGL3_DEMOS OFF)
option(BUILD_UNIT_TESTS OFF)
add_subdirectory(LOGLEngine/Vendor/bullet)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -std=c++11")
    if(NOT WIN32)
        set(GLAD_LIBRARIES dl)
    endif()
endif()

include_directories(LOGLEngine/Headers/
                    LOGLEngine/Vendor/assimp/include/
                    LOGLEngine/Vendor/bullet/src/
                    LOGLEngine/Vendor/glad/include/
                    LOGLEngine/Vendor/glfw/include/
                    LOGLEngine/Vendor/glm/
                    LOGLEngine/Vendor/stb/)

file(GLOB VENDORS_SOURCES LOGLEngine/Vendor/glad/src/glad.c)
file(GLOB PROJECT_HEADERS LOGLEngine/Headers/*.hpp
                          LOGLEngine/Headers/*.h)
file(GLOB PROJECT_SOURCES LOGLEngine/Sources/*.cpp
                          LOGLEngine/Sources/*.c)
file(GLOB PROJECT_SHADERS LOGLEngine/Shaders/*.comp
                          LOGLEngine/Shaders/*.frag
                          LOGLEngine/Shaders/*.geom
                          LOGLEngine/Shaders/*.vert)
file(GLOB PROJECT_CONFIGS CMakeLists.txt
                          Readme.md
                         .gitattributes
                         .gitignore
                         .gitmodules)

source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Shaders" FILES ${PROJECT_SHADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})
source_group("Vendors" FILES ${VENDORS_SOURCES})

add_definitions(-DGLFW_INCLUDE_NONE
                -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS}
                               ${PROJECT_SHADERS} ${PROJECT_CONFIGS}
                               ${VENDORS_SOURCES})
target_link_libraries(${PROJECT_NAME} assimp glfw
                      ${GLFW_LIBRARIES} ${GLAD_LIBRARIES}
                      BulletDynamics BulletCollision LinearMath)
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
