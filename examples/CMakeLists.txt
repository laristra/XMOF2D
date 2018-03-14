#[[
This file is part of the Ristra XMOF2D project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/XMOF2D/blob/master/LICENSE

Created by Evgeny Kikinzon.
Copyright © 2018, Los Alamos National Security, LLC. 
All rights reserved.
]]

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(test_line_interface src/test_line_interface.cpp)
add_executable(single_cell src/single_cell.cpp)
if (${ARCHOS} STREQUAL "x86_64_Darwin")
  target_link_libraries(test_line_interface ${XMOF2DLIB} ${EXTRA_LIBS} c++)
  target_link_libraries(single_cell ${XMOF2DLIB} ${EXTRA_LIBS} c++)
else ()
  target_link_libraries(test_line_interface ${XMOF2DLIB} ${EXTRA_LIBS})
  target_link_libraries(single_cell ${XMOF2DLIB} ${EXTRA_LIBS})
endif()
install(TARGETS test_line_interface DESTINATION ${PROJECT_SOURCE_DIR}/examples)
install(TARGETS single_cell DESTINATION ${PROJECT_SOURCE_DIR}/examples)
