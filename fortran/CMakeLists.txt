#[[
This file is part of the Ristra XMOF2D project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/XMOF2D/blob/master/LICENSE

Created by Evgeny Kikinzon.
Copyright © 2018, Triad National Security, LLC. 
All rights reserved.
]]

add_executable(single_cell_fortran $<TARGET_OBJECTS:fxmof2D_single_mmc> src/single_cell_fortran.f)
set_property(TARGET single_cell_fortran PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(single_cell_fortran ${XMOF2DLIB})
install(TARGETS single_cell_fortran DESTINATION ${PROJECT_SOURCE_DIR}/fortran)
