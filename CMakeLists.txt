cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(vtkdemo)
add_executable(cubic2d cubic2d.cpp)
add_executable(cubic3d cubic3d.cpp)
#add_executable(tetr3d tetr3d.cpp)
add_executable(pikachu2 pikachu2.cpp)
#add_executable(porshe porshe.cpp)

# Найдём VTK
find_package(VTK REQUIRED)
# ... и подключим её со всеми её компонентами её нативными средствами
include(${VTK_USE_FILE})

# gmsh укажем явно по кусочкам (потому что он был просто распакован руками)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include"
)
set(GMSH_LIBRARY gmsh)

# если бы gmsh был поставлен прилично, то можно было бы примерно вот так
# find_package(gmsh REQUIRED)
# include_directories(${GMSH_INCLUDE_DIR})
# target_link_libraries(tetr3d ${GMSH_LIBRARY})

# Укажем библиотеки для линковки
# ... тут только VTK
target_link_libraries(cubic2d ${VTK_LIBRARIES})
target_link_libraries(cubic3d ${VTK_LIBRARIES})

# ... а тут ещё и gmsh
target_link_libraries(pikachu2 ${VTK_LIBRARIES})
target_link_libraries(pikachu2 ${GMSH_LIBRARY})


#target_link_libraries(porshe ${VTK_LIBRARIES})
#target_link_libraries(porshe ${GMSH_LIBRARY})
