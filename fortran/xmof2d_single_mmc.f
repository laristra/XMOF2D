      ! This file is part of the Ristra XMOF2D project.
      ! Please see the license file at the root of this repository, or at:
      !     https://github.com/laristra/XMOF2D/blob/master/LICENSE

      ! Created by Evgeny Kikinzon.
      ! Copyright © 2018, Triad National Security, LLC.
      ! All rights reserved.

      module xmof2d_single_mmc

      use iso_c_binding
      INTERFACE

        SUBROUTINE xmof2d_set_distance_tolerance(dist_eps) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), VALUE, INTENT(IN) :: dist_eps
        END SUBROUTINE xmof2d_set_distance_tolerance

        SUBROUTINE xmof2d_get_distance_tolerance(dist_eps) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: dist_eps
        END SUBROUTINE xmof2d_get_distance_tolerance

        SUBROUTINE xmof2d_set_area_tolerance(area_eps) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), VALUE, INTENT(IN) :: area_eps
        END SUBROUTINE xmof2d_set_area_tolerance

        SUBROUTINE xmof2d_get_area_tolerance(area_eps) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: area_eps
        END SUBROUTINE xmof2d_get_area_tolerance

        SUBROUTINE xmof2d_set_angle_tolerance(ang_eps) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), VALUE, INTENT(IN) :: ang_eps
        END SUBROUTINE xmof2d_set_angle_tolerance

        SUBROUTINE xmof2d_get_angle_tolerance(ang_eps) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          REAL(KIND=C_DOUBLE), INTENT(OUT) :: ang_eps
        END SUBROUTINE xmof2d_get_angle_tolerance

        SUBROUTINE xmof2d_set_max_iter_num(mof_max_iter) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: mof_max_iter
        END SUBROUTINE xmof2d_set_max_iter_num

        SUBROUTINE xmof2d_get_max_iter_num(mof_max_iter) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), INTENT(OUT) :: mof_max_iter
        END SUBROUTINE xmof2d_get_max_iter_num

        SUBROUTINE xmof2d_set_mmc_vertices(nvrts, crds) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN) :: nvrts
          REAL(KIND=C_DOUBLE), DIMENSION(*), INTENT(IN) :: crds
        END SUBROUTINE xmof2d_set_mmc_vertices

        SUBROUTINE xmof2d_set_materials_data(nmaterials, mat_ids, 
     &vol_fractions, centroids) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN) :: nmaterials
          INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN) :: mat_ids
          REAL(KIND=C_DOUBLE), DIMENSION(*), INTENT(IN) :: vol_fractions
          REAL(KIND=C_DOUBLE), DIMENSION(*), INTENT(IN) :: centroids
        END SUBROUTINE xmof2d_set_materials_data

        SUBROUTINE xmof2d_initialize_reconstructor() BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
        END SUBROUTINE xmof2d_initialize_reconstructor

        SUBROUTINE xmof2d_free_reconstructor() BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
        END SUBROUTINE xmof2d_free_reconstructor

        SUBROUTINE xmof2d_perform_reconstruction(succeeded) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          LOGICAL(KIND=C_BOOL), INTENT(OUT) :: succeeded          
        END SUBROUTINE xmof2d_perform_reconstruction

        SUBROUTINE xmof2d_mesh_get_ncells(ncells) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), INTENT(OUT) :: ncells
        END SUBROUTINE xmof2d_mesh_get_ncells

        SUBROUTINE xmof2d_mesh_get_nfaces(nfaces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), INTENT(OUT) :: nfaces
        END SUBROUTINE xmof2d_mesh_get_nfaces

        SUBROUTINE xmof2d_mesh_get_nnodes(nnodes) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), INTENT(OUT) :: nnodes
        END SUBROUTINE xmof2d_mesh_get_nnodes

        SUBROUTINE xmof2d_cell_get_mat_id(icell, mat_id) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: icell
          INTEGER(KIND=C_INT),        INTENT(OUT) :: mat_id
        END SUBROUTINE xmof2d_cell_get_mat_id

        SUBROUTINE xmof2d_cell_get_nfaces(icell, nfaces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: icell
          INTEGER(KIND=C_INT),        INTENT(OUT) :: nfaces
        END SUBROUTINE xmof2d_cell_get_nfaces

        SUBROUTINE xmof2d_cell_get_face_ids(icell, ifaces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN)  :: icell
          INTEGER(KIND=C_INT), DIMENSION(*), INTENT(OUT) :: ifaces
        END SUBROUTINE xmof2d_cell_get_face_ids

        SUBROUTINE xmof2d_cell_get_node_ids(icell, inodes) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN)  :: icell
          INTEGER(KIND=C_INT), DIMENSION(*), INTENT(OUT) :: inodes
        END SUBROUTINE xmof2d_cell_get_node_ids

        SUBROUTINE xmof2d_cell_get_size(icell, area) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: icell
          REAL(KIND=C_DOUBLE),        INTENT(OUT) :: area
        END SUBROUTINE xmof2d_cell_get_size

        SUBROUTINE xmof2d_cell_get_center(icell, cen_coords) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,          INTENT(IN)  :: icell
          REAL(KIND=C_DOUBLE), DIMENSION(0:1), INTENT(OUT) :: cen_coords
        END SUBROUTINE xmof2d_cell_get_center

        SUBROUTINE xmof2d_face_get_cell_ids(iface, 
     &left_cell_id, right_cell_id) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: iface
          INTEGER(KIND=C_INT),        INTENT(OUT) :: left_cell_id
          INTEGER(KIND=C_INT),        INTENT(OUT) :: right_cell_id
        END SUBROUTINE xmof2d_face_get_cell_ids

        SUBROUTINE xmof2d_face_get_node_ids(iface,
     &tail_id, head_id) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: iface
          INTEGER(KIND=C_INT),        INTENT(OUT) :: tail_id
          INTEGER(KIND=C_INT),        INTENT(OUT) :: head_id
        END SUBROUTINE xmof2d_face_get_node_ids

        SUBROUTINE xmof2d_face_is_boundary(iface, is_boundary) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: iface
          LOGICAL(KIND=C_BOOL),       INTENT(OUT) :: is_boundary
        END SUBROUTINE xmof2d_face_is_boundary

        SUBROUTINE xmof2d_face_get_parent_face_id(iface, 
     &iparent) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: iface
          INTEGER(KIND=C_INT),        INTENT(OUT) :: iparent
        END SUBROUTINE xmof2d_face_get_parent_face_id

        SUBROUTINE xmof2d_face_get_size(iface, length) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: iface
          REAL(KIND=C_DOUBLE),        INTENT(OUT) :: length
        END SUBROUTINE xmof2d_face_get_size

        SUBROUTINE xmof2d_face_get_center(iface, cen_coords) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,          INTENT(IN)  :: iface
          REAL(KIND=C_DOUBLE), DIMENSION(0:1), INTENT(OUT) :: cen_coords
        END SUBROUTINE xmof2d_face_get_center

        SUBROUTINE xmof2d_node_get_coords(inode, coords) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,          INTENT(IN)  :: inode
          REAL(KIND=C_DOUBLE), DIMENSION(0:1), INTENT(OUT) :: coords
        END SUBROUTINE xmof2d_node_get_coords

        SUBROUTINE xmof2d_node_get_ncells(inode, ncells) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: inode
          INTEGER(KIND=C_INT),        INTENT(OUT) :: ncells
        END SUBROUTINE xmof2d_node_get_ncells

        SUBROUTINE xmof2d_node_get_cell_ids(inode, icells) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN)  :: inode
          INTEGER(KIND=C_INT), DIMENSION(*), INTENT(OUT) :: icells
        END SUBROUTINE xmof2d_node_get_cell_ids

        SUBROUTINE xmof2d_node_get_nfaces(inode, nfaces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: inode
          INTEGER(KIND=C_INT),        INTENT(OUT) :: nfaces
        END SUBROUTINE xmof2d_node_get_nfaces

        SUBROUTINE xmof2d_node_get_face_ids(inode, ifaces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN)  :: inode
          INTEGER(KIND=C_INT), DIMENSION(*), INTENT(OUT) :: ifaces
        END SUBROUTINE xmof2d_node_get_face_ids

        SUBROUTINE xmof2d_node_is_boundary(inode, is_boundary) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: inode
          LOGICAL(KIND=C_BOOL),       INTENT(OUT) :: is_boundary
        END SUBROUTINE xmof2d_node_is_boundary

        SUBROUTINE xmof2d_node_get_parent_node_id(inode, 
     &iparent) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: inode
          INTEGER(KIND=C_INT),        INTENT(OUT) :: iparent
        END SUBROUTINE xmof2d_node_get_parent_node_id

        SUBROUTINE xmof2d_node_get_nparent_faces(inode,
     &nparent_faces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE, INTENT(IN)  :: inode
          INTEGER(KIND=C_INT),        INTENT(OUT) :: nparent_faces
        END SUBROUTINE xmof2d_node_get_nparent_faces

        SUBROUTINE xmof2d_node_get_parent_face_ids(inode, 
     &iparent_faces) BIND(C)
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE,        INTENT(IN)  :: inode
          INTEGER(KIND=C_INT), DIMENSION(*), INTENT(OUT) :: 
     &iparent_faces
        END SUBROUTINE xmof2d_node_get_parent_face_ids

      END INTERFACE

      end module xmof2d_single_mmc
