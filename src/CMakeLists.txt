cmake_minimum_required( VERSION 3.5.0 )
project( tutorials C CXX )

set( CMAKE_BUILD_TYPE Debug )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_FLAGS_DEBUG "-g" )
set( CMAKE_CXX_FLAGS_RELEASE "-O3" )

add_definitions( -DRESOURCES_PATH="${CMAKE_SOURCE_DIR}/res/" 
                 -DRESULTS_PATH="${CMAKE_SOURCE_DIR}/results/" )

#message( "resources path: ${CMAKE_SOURCE_DIR}/res/" )
#message( "results path: ${CMAKE_SOURCE_DIR}/results/" )

#set(OpenCV_DIR "/data/leticia/.ucsp/lib/python3.6/site-packages/cv2")
#set(OpenCV_DIR) #/data/leticia/.ucsp/lib/python3.6/site-packages/cv2/

find_package(OpenCV REQUIRED)

message(${OpenCV_DIR})
message(OPENCV_INCLUDE_DIRS)
message(OpenCV_FOUND)
include_directories( ${OpenCV_INCLUDE_DIRS} )

###############################################################

## set( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/libs" )

# generate helpers library
add_library( leti_funcs SHARED
             functions.cpp )
target_link_libraries( leti_funcs
                       ${OpenCV_LIBS} )

###############################################################

add_executable( main
                main.cpp )
target_link_libraries( main
                       ${OpenCV_LIBS}
                       leti_funcs )