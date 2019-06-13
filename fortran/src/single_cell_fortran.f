      ! This file is part of the Ristra XMOF2D project.
      ! Please see the license file at the root of this repository, or at:
      !     https://github.com/laristra/XMOF2D/blob/master/LICENSE

      ! Created by Evgeny Kikinzon.
      ! Copyright © 2018, Triad National Security, LLC.
      ! All rights reserved.

      program single_cell_fortran
      use xmof2d_single_mmc

      implicit none

      integer ivrt, imat, id
      integer nbase_vrts, nmats, max_iter_num
      real*8  dist_tol, ddot_tol, area_tol, ang_tol
      real*8, dimension(12)  :: base_cell_vrts
      integer, dimension(4)  :: mat_ids
      real*8, dimension(4)   :: mat_vfracs
      real*8, dimension(8)   :: mat_centroids

      integer ncells, nfaces, nnodes
      integer icell, iface, inode
      integer cell_mat_id, cell_nfaces, iref, iparent
      integer itail, ihead, ileft_cell, iright_cell
      integer node_ncells, node_nfaces 
      integer iparent_node, node_nparent_faces
      real*8  cell_area, cen_diff, face_length
      logical*1 rec_success, is_boundary
      real*8, dimension(1:2) :: cell_center, face_center
      real*8, dimension(1:2) :: node_crd
      integer, dimension (:), allocatable :: cell_ifaces
      integer, dimension (:), allocatable :: cell_inodes
      integer, dimension (:), allocatable :: node_icells
      integer, dimension (:), allocatable :: node_ifaces
      integer, dimension (:), allocatable :: node_iparent_faces

      ! Initialize base cell vertices
      base_cell_vrts = (/
     & 8.517743d-01, 3.781231d-01,
     & 8.506274d-01, 3.789720d-01,
     & 8.486559d-01, 3.770629d-01,
     & 8.490196d-01, 3.748427d-01,
     & 8.505602d-01, 3.748087d-01,
     & 8.519694d-01, 3.762168d-01 /)
      nbase_vrts = 6

      ! Initialize material IDs
      mat_ids = (/ 2, 3, 5, 7 /)
      nmats = 4

      ! Initialize volume fractions
      mat_vfracs = (/ 1.862369d-01, 4.742934d-01, 
     & 8.186527d-02, 2.576044d-01 /)

      ! Initialize reference centroids
      mat_centroids = (/
     & 8.501062d-01, 3.779714d-01,
     & 8.501696d-01, 3.767483d-01,
     & 8.493403d-01, 3.752849d-01,
     & 8.511479d-01, 3.761662d-01 /)

      print*, "****************************************"
      print*, "Input data:"
      print*, "Vertices of the base cell"
      do ivrt = 1, nbase_vrts
        print*, "Vertex #", ivrt, " at (", 
     & base_cell_vrts(2*ivrt - 1),
     & ", ", base_cell_vrts(2*ivrt), ")"
      end do     
      !write(10,*) "" 
      print*, ""
      print*, "IDs of materials inside the base cell: ";
      do imat = 1, nmats
        print*, mat_ids(imat), " "
      end do
      print*, ""
      print*, "Volume fractions of materials inside the base cell"
      do imat = 1, nmats
        print*, "Material #", mat_ids(imat), ": ", mat_vfracs(imat)
      end do
      print*, ""
      print*, "Reference centroids of materials inside the base cell"
      do imat = 1, nmats
        print*, "Material #", mat_ids(imat), ": (",
     & mat_centroids(2*imat - 1), ", ", mat_centroids(2*imat), ")"
      end do
      print*, ""

      open(unit=10, file="mmc_node_crds.txt")
      do ivrt = 1, nbase_vrts + 1
        id = mod(ivrt - 1, nbase_vrts) + 1
        write(10,*) base_cell_vrts(2*id - 1), base_cell_vrts(2*id)
      end do
      close(unit=10)
      open(unit=20, file="gnuplot_script_mmc.gnu")
  210 format ('set label "', I0, '" at ', E21.16, ', ', E21.16
     & ' tc lt 4')
  211 format ('set label "', I0, '" at ', E21.16, ', ', E21.16 
     & ' right tc lt 3')
      do ivrt = 1, nbase_vrts
        write(20,210) ivrt - 1,  
     & base_cell_vrts(2*ivrt - 1), base_cell_vrts(2*ivrt)
        id = mod(ivrt, nbase_vrts) + 1 
        face_center(1) = 0.5*(base_cell_vrts(2*ivrt - 1) + 
     & base_cell_vrts(2*id - 1))
        face_center(2) = 0.5*(base_cell_vrts(2*ivrt) +
     & base_cell_vrts(2*id))
        write(20,211) ivrt - 1, face_center(1), face_center(2)
      end do
      write(20,*) "set size square"
      write(20,*) "plot 'mmc_node_crds.txt' w lp notitle"
      close(unit=20)

      print*, "****************************************"
      print*, "Setting coordinates of multi-material cell's vertices..."
      call xmof2d_set_mmc_vertices(nbase_vrts, base_cell_vrts)
      print*, "Setting IDs, volume fractions, and centroids
     & for materials..."
      call xmof2d_set_materials_data(nmats, mat_ids, mat_vfracs, 
     & mat_centroids)

      call xmof2d_get_distance_tolerance(dist_tol)
      print*, "Default distance tolerance: ", dist_tol
      call xmof2d_get_dot_product_tolerance(ddot_tol)
      print*, "Default dot product tolerance: ", ddot_tol
      call xmof2d_get_area_tolerance(area_tol)
      print*, "Default area tolerance: ", area_tol
      call xmof2d_get_angle_tolerance(ang_tol)
      print*, "Default angle tolerance: ", ang_tol
      call xmof2d_get_max_iter_num(max_iter_num)
      print*, "Default max number of iterations: ", max_iter_num
  
      area_tol = 1.0d-14
      print*, "Relaxing area tolerance to ", area_tol
      call xmof2d_set_area_tolerance(area_tol);

      print*, "Creating a reconstructor instance..."
      call xmof2d_initialize_reconstructor()
      print*, "Constructing a minimesh..."
      call xmof2d_perform_reconstruction(rec_success)
      if (.not.rec_success) then
        print*, "Reconstruction has failed, exiting..."
        error stop
      endif      
      print*, "****************************************"

      call xmof2d_mesh_get_ncells(ncells)
      call xmof2d_mesh_get_nfaces(nfaces)
      call xmof2d_mesh_get_nnodes(nnodes)
      print*, "Resulting minimesh has ", ncells, " cells, ",
     & nfaces, " faces, and ", nnodes, " nodes."
      print*, "****************************************"

      open(unit=30, file="matpolys_node_crds.txt")
      open(unit=40, file="gnuplot_script_minimesh.gnu")
  410 format ('set label "#', I0, ',MID ', I0, '" at ',
     & E21.16, ', ', E21.16, ' center tc lt 1')
      print*, "Cells data:"
      do icell = 1, ncells
        print*, "Cell #", icell
        call xmof2d_cell_get_mat_id(icell - 1, cell_mat_id)
        print*, "Material ID: ", cell_mat_id
        call xmof2d_cell_get_size(icell - 1, cell_area)
        print*, "Area: ", cell_area
        call xmof2d_cell_get_center(icell - 1, cell_center)
        print*, "Actual Centroid: (", cell_center(1), ", ",
     & cell_center(2), ")"
        write(40,410) icell, cell_mat_id, 
     & cell_center(1), cell_center(2)
        iref = 0
        do imat = 1, nmats
          if (mat_ids(imat) == cell_mat_id) then
            iref = imat
            exit
          endif
        end do
        cen_diff = sqrt((cell_center(1) - mat_centroids(2*iref - 1))**2 
     & + (cell_center(2) - mat_centroids(2*iref))**2)
        print*, "Discrepancy in centroids: ", cen_diff
        call xmof2d_cell_get_nfaces(icell - 1, cell_nfaces)
        print*, "Number of faces: ", cell_nfaces
    
        allocate(cell_ifaces(cell_nfaces))
        call xmof2d_cell_get_face_ids(icell - 1, cell_ifaces)
        cell_ifaces = cell_ifaces + 1
        print*, "Indices of faces: ", cell_ifaces(1:cell_nfaces)
        deallocate(cell_ifaces)
      
        allocate(cell_inodes(cell_nfaces))
        call xmof2d_cell_get_node_ids(icell - 1, cell_inodes)
        cell_inodes = cell_inodes + 1
        print*, "Indices of nodes: ", cell_inodes(1:cell_nfaces)
        print*, ""

        do ivrt = 1, cell_nfaces + 1
          id = mod(ivrt - 1, cell_nfaces) + 1
          call xmof2d_node_get_coords(cell_inodes(id) - 1, node_crd)
          write(30,*) node_crd(1), node_crd(2)
        end do
        write(30,*) ""
        deallocate(cell_inodes)      
      end do
      close(unit=30)

      print*, "****************************************"
      print*, "Faces data:"
      do iface = 1, nfaces
        print*, "Face #", iface
        call xmof2d_face_get_node_ids(iface - 1, itail, ihead)
        itail = itail + 1
        ihead = ihead + 1
        print*, "Index of the tail node: ", itail
        print*, "Index of the head node: ", ihead
        call xmof2d_face_get_cell_ids(iface - 1, 
     & ileft_cell, iright_cell)
        ileft_cell = ileft_cell + 1
        iright_cell = iright_cell + 1
        print*, "Index of the cell to the left: ", ileft_cell
        print*, "Index of the cell to the right: "
        if (iright_cell == 0) then
          print*, "NONE"
        else
          print*, iright_cell
        endif
        call xmof2d_face_get_size(iface - 1, face_length)
        print*, "Length: ", face_length
        call xmof2d_face_get_center(iface - 1, face_center)
        print*, "Center: (", face_center(1), ", ",
     & face_center(2), ")"

        write(40,211) iface, face_center(1), face_center(2)

        call xmof2d_face_is_boundary(iface - 1, is_boundary)
        print*, "This face is ";
        if (.not.is_boundary) then
          print*, "NOT "
        endif
        print*, "boundary"
    
        call xmof2d_face_get_parent_face_id(iface - 1, iparent)
        iparent = iparent + 1
        print*, "Parent face ID: "
        if (iparent == 0) then
          print*, "NONE"
        else
          print*, iparent
        endif
        print*, ""
      end do
      print*, "****************************************"
      print*, "Nodes data:"
      do inode = 1, nnodes
        print*, "Node #", inode
        call xmof2d_node_get_coords(inode - 1, node_crd)
        print*, "Coordinates: (", node_crd(1), ", ",
     & node_crd(2), ")"

        write(40,210) inode, node_crd(1), node_crd(2)

        call xmof2d_node_get_ncells(inode - 1, node_ncells)
        print*, "This is a node of ", node_ncells, " cell(s)"
        allocate(node_icells(node_ncells))
        call xmof2d_node_get_cell_ids(inode - 1, node_icells)
        node_icells = node_icells + 1
        print*, "Indices of cells: ", node_icells(1:node_ncells)
        deallocate( node_icells)

        call xmof2d_node_get_nfaces(inode - 1, node_nfaces)
        print*, "This is a node of ", node_nfaces, " face(s)"
        allocate(node_ifaces(node_nfaces))
        call xmof2d_node_get_face_ids(inode - 1, node_ifaces)
        node_ifaces = node_ifaces + 1
        print*, "Indices of faces: ", node_ifaces(1:node_nfaces)
        deallocate(node_ifaces)
    
        call xmof2d_node_is_boundary(inode - 1, is_boundary)
        print*, "This node is "
        if (.not.is_boundary) then
          print*, "NOT "
        endif
        print*, "boundary"
    
        call xmof2d_node_get_parent_node_id(inode - 1, iparent_node)
        iparent_node = iparent_node + 1
        if (iparent_node == 0) then
          print*, "Parent node ID: ", "NONE"
        else
          print*, "Parent node ID: ", iparent_node
        endif
    
        call xmof2d_node_get_nparent_faces(inode - 1, 
     & node_nparent_faces)
        if (node_nparent_faces == 0) then
          print*, "This node belongs to", "NO", " base face(s)"
        else
          print*, "This node belongs to", node_nparent_faces,
     & " base face(s)"
        endif
    
        if (node_nparent_faces /= 0) then
          allocate(node_iparent_faces(node_nparent_faces))
          call xmof2d_node_get_parent_face_ids(inode - 1, 
     & node_iparent_faces)
          node_iparent_faces = node_iparent_faces + 1
          print*, "Indices of those base faces: ", 
     & node_iparent_faces(1:node_nparent_faces)
          deallocate( node_iparent_faces)
        endif
      end do

      write(40,*) "set size square"
      write(40,*) "plot 'matpolys_node_crds.txt' w l notitle"
      close(unit=40)

      print*, "****************************************"
      print*, "Removing the reconstructor instance..."
      call   xmof2d_free_reconstructor()
      print*, "All done!"

      stop
      end program single_cell_fortran
