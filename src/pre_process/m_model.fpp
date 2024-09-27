!>
!! @file   m_model.fpp
!! @author Henry Le Berre <hberre3@gatech.edu>
!! @brief  Contains module m_model

#:include 'macros.fpp'

module m_model

    use m_helper
    use m_mpi_proxy
    use m_derived_types

    use iso_c_binding, only: c_char, c_int32_t, c_int16_t, c_float

    implicit none

    private

    public :: f_model_read, s_model_write, s_model_free, f_model_is_inside, & 
              f_check_boundary, f_distance, register_edge, f_normals, & 
              f_distance_normals_3D, f_interpolate_2D, f_interpolate_3D, &
              f_interpolated_distance, & 
              f_check_interpolation_2D, &
              f_check_interpolation_3D

contains

    !> This procedure reads a binary STL file.
    subroutine s_read_stl_binary(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(OUT) :: model

        integer :: i, j, iunit, iostat

        character(kind=c_char, len=80) :: header
        integer(kind=c_int32_t) :: nTriangles

        real(kind=c_float) :: normal(3), v(3, 3)
        integer(kind=c_int16_t) :: attribute

        open (newunit=iunit, file=filepath, action='READ', &
              form='UNFORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open Binary STL file ", filepath

            call s_mpi_abort()
        end if

        read (iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not read header from Binary STL file ", filepath

            call s_mpi_abort()
        end if

        model%ntrs = nTriangles

        allocate (model%trs(model%ntrs))

        do i = 1, model%ntrs
            read (iunit) normal(:), v(1, :), v(2, :), v(3, :), attribute

            model%trs(i)%v = v
            model%trs(i)%n = normal
        end do

        close (iunit)

    end subroutine s_read_stl_binary

    !> This procedure reads an ASCII STL file.
    subroutine s_read_stl_ascii(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(OUT) :: model

        integer :: i, j, iunit, iostat

        character(80) :: line

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open ASCII STL file ", filepath

            call s_mpi_abort()
        end if

        model%ntrs = 0
        do
            if (.not. f_read_line(iunit, line)) exit

            if (line(1:6) == "facet ") then
                model%ntrs = model%ntrs + 1
            end if
        end do

        allocate (model%trs(model%ntrs))

        rewind (iunit)

        i = 1
        do
            if (.not. f_read_line(iunit, line)) exit

            if (line(1:5) == "solid") cycle
            if (line(1:8) == "endsolid") exit

            if (line(1:12) /= "facet normal") then
                print *, "Error: expected facet normal in STL file ", filepath

                call s_mpi_abort()
            end if

            call s_skip_ignored_lines(iunit)
            read (line(13:), *) model%trs(i)%n

            call s_skip_ignored_lines(iunit)
            read (iunit, '(A)') line

            do j = 1, 3
                if (.not. f_read_line(iunit, line)) exit

                if (line(1:6) /= "vertex") then
                    print *, "Error: expected vertex in STL file ", filepath

                    call s_mpi_abort()
                end if

                call s_skip_ignored_lines(iunit)
                read (line(7:), *) model%trs(i)%v(j, :)
            end do

            if (.not. f_read_line(iunit, line)) exit
            if (.not. f_read_line(iunit, line)) exit

            if (line(1:8) /= "endfacet") then
                print *, "Error: expected endfacet in STL file ", filepath

                call s_mpi_abort()
            end if

            i = i + 1
        end do

    end subroutine s_read_stl_ascii

    !> This procedure reads an STL file.
    subroutine s_read_stl(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(OUT) :: model

        integer :: iunit, iostat

        character(80) :: line

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath

            call s_mpi_abort()
        end if

        read (iunit, '(A)') line

        close (iunit)

        if (line(1:5) == "solid") then
            call s_read_stl_ascii(filepath, model)
        else
            call s_read_stl_binary(filepath, model)
        end if

    end subroutine

    !> This procedure reads an OBJ file.
    subroutine s_read_obj(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(OUT) :: model

        integer :: i, j, k, l, iunit, iostat, nVertices

        t_vec3, allocatable :: vertices(:, :)

        character(80) :: line

        open (newunit=iunit, file=filepath, action='READ', &
              form='FORMATTED', status='OLD', iostat=iostat, &
              access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open model file ", filepath

            call s_mpi_abort()
        end if

        nVertices = 0
        model%ntrs = 0
        do
            if (.not. f_read_line(iunit, line)) exit

            select case (line(1:2))
            case ("v ")
                nVertices = nVertices + 1
            case ("f ")
                model%ntrs = model%ntrs + 1
            end select
        end do

        rewind (iunit)

        allocate (vertices(nVertices, 1:3))
        allocate (model%trs(model%ntrs))

        i = 1
        j = 1

        do
            if (.not. f_read_line(iunit, line)) exit

            select case (line(1:2))
            case ("g ")
            case ("vn")
            case ("vt")
            case ("l ")
            case ("v ")
                read (line(3:), *) vertices(i, :)
                i = i + 1
            case ("f ")
                read (line(3:), *) k, l, j
                model%trs(j)%v(1, :) = vertices(k, :)
                model%trs(j)%v(2, :) = vertices(l, :)
                model%trs(j)%v(3, :) = vertices(j, :)
                j = j + 1
            case default
                print *, "Error: unknown line type in OBJ file ", filepath
                print *, "Line: ", line

                call s_mpi_abort()
            end select
        end do

        deallocate (vertices)

        close (iunit)

    end subroutine

    !> This procedure reads a mesh from a file.
    !! @param filepath Path to the file to read.
    !! @return The model read from the file.
    function f_model_read(filepath) result(model)

        character(LEN=*), intent(IN) :: filepath

        type(t_model) :: model

        select case (filepath(len(trim(filepath)) - 3:len(trim(filepath))))
        case (".stl")
            call s_read_stl(filepath, model)
        case (".obj")
            call s_read_obj(filepath, model)
        case default
            print *, "Error: unknown model file format for file ", filepath

            call s_mpi_abort()
        end select

    end function f_model_read

    !> This procedure writes a binary STL file.
    subroutine s_write_stl(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(IN) :: model

        integer :: i, j, iunit, iostat

        character(kind=c_char, len=80), parameter :: header = "Model file written by MFC."
        integer(kind=c_int32_t) :: nTriangles
        real(kind=c_float) :: normal(3), v(3)
        integer(kind=c_int16_t) :: attribute

        open (newunit=iunit, file=filepath, action='WRITE', &
              form='UNFORMATTED', iostat=iostat, access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open STL file ", filepath

            call s_mpi_abort()
        end if

        nTriangles = model%ntrs
        write (iunit, iostat=iostat) header, nTriangles

        if (iostat /= 0) then
            print *, "Error: could not write header to STL file ", filepath

            call s_mpi_abort()
        end if

        do i = 1, model%ntrs
            normal = model%trs(i)%n
            write (iunit) normal

            do j = 1, 3
                v = model%trs(i)%v(j, :)
                write (iunit) v(:)
            end do

            attribute = 0
            write (iunit) attribute
        end do

        close (iunit)

    end subroutine s_write_stl

    !> This procedure writes an OBJ file.
    subroutine s_write_obj(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(IN) :: model

        integer :: iunit, iostat

        integer :: i, j

        open (newunit=iunit, file=filepath, action='WRITE', &
              form='FORMATTED', iostat=iostat, access='STREAM')

        if (iostat /= 0) then
            print *, "Error: could not open OBJ file ", filepath

            call s_mpi_abort()
        end if

        write (iunit, '(A)') "# Model file written by MFC."

        do i = 1, model%ntrs
            do j = 1, 3
                write (iunit, '(A, " ", (f30.20), " ", (f30.20), " ", (f30.20))') &
                    "v", model%trs(i)%v(j, 1), model%trs(i)%v(j, 2), model%trs(i)%v(j, 3)
            end do

            write (iunit, '(A, " ", I0, " ", I0, " ", I0)') &
                "f", i*3 - 2, i*3 - 1, i*3
        end do

        close (iunit)

    end subroutine s_write_obj

    !> This procedure writes a binary STL file.
    !! @param filepath  Path to the file to write.
    !! @param triangles Triangles to write.
    subroutine s_model_write(filepath, model)

        character(LEN=*), intent(IN) :: filepath
        type(t_model), intent(IN) :: model

        select case (filepath(len(trim(filepath)) - 3:len(trim(filepath))))
        case (".stl")
            call s_write_stl(filepath, model)
        case (".obj")
            call s_write_obj(filepath, model)
        case default
            print *, "Error: unknown model file format for file ", filepath

            call s_mpi_abort()
        end select

    end subroutine s_model_write

    !> This procedure frees the memory allocated for an STL mesh.
    subroutine s_model_free(model)

        type(t_model), intent(INOUT) :: model

        deallocate (model%trs)

    end subroutine s_model_free

    function f_read_line(iunit, line) result(bIsLine)

        integer, intent(IN) :: iunit
        character(80), intent(OUT) :: line
        logical :: bIsLine

        integer :: iostat

        bIsLine = .true.

        do
            read (iunit, '(A)', iostat=iostat) line

            if (iostat < 0) then
                bIsLine = .false.
                exit
            end if

            line = adjustl(trim(line))

            if (len(trim(line)) == 0) cycle
            if (line(1:5) == "solid") cycle
            if (line(1:1) == "#") cycle

            exit
        end do

    end function f_read_line

    subroutine s_skip_ignored_lines(iunit)

        integer, intent(IN) :: iunit

        character(80) :: line

        if (f_read_line(iunit, line)) then
            backspace (iunit)
        end if

    end subroutine s_skip_ignored_lines

    !> This procedure, recursively, finds whether a point is inside an octree.
    !! @param model    Model to search in.
    !! @param point    Point to test.
    !! @param spacing  Space around the point to search in (grid spacing).
    !! @param spc      Number of samples per cell.
    !! @return True if the point is inside the octree, false otherwise.
    function f_model_is_inside(model, point, spacing, spc) result(fraction)

        type(t_model), intent(in) :: model
        t_vec3, intent(in) :: point
        t_vec3, intent(in) :: spacing
        integer, intent(in) :: spc

        real(kind(0d0)) :: fraction

        type(t_ray) :: ray
        integer :: i, j, nInOrOut, nHits

        real(kind(0d0)), dimension(1:spc, 1:3) :: ray_origins, ray_dirs

        do i = 1, spc
            call random_number(ray_origins(i, :))
            ray_origins(i, :) = point + (ray_origins(i, :) - 0.5)*spacing(:)

            call random_number(ray_dirs(i, :))
            ray_dirs(i, :) = ray_dirs(i, :) - 0.5
            ray_dirs(i, :) = ray_dirs(i, :)/sqrt(sum(ray_dirs(i, :)*ray_dirs(i, :)))
        end do

        nInOrOut = 0
        do i = 1, spc
            ray%o = ray_origins(i, :)
            ray%d = ray_dirs(i, :)

            nHits = 0
            do j = 1, model%ntrs
                if (f_intersects_triangle(ray, model%trs(j))) then
                    nHits = nHits + 1
                end if
            end do

            nInOrOut = nInOrOut + mod(nHits, 2)
        end do

        fraction = real(nInOrOut)/real(spc)

    end function f_model_is_inside

    ! From https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
    !> This procedure checks if a ray intersects a triangle.
    !! @param ray      Ray.
    !! @param triangle Triangle.
    !! @return         True if the ray intersects the triangle, false otherwise.
    function f_intersects_triangle(ray, triangle) result(intersects)

        type(t_ray), intent(in) :: ray
        type(t_triangle), intent(in) :: triangle

        logical :: intersects

        real(kind(0d0)) :: v0v1(3), v0v2(3), N(3), P(3), C(3), edge(3), vp(3)
        real(kind(0d0)) :: area2, d, t, NdotRayDirection

        intersects = .false.

        N = triangle%n
        area2 = sqrt(sum(N(:)*N(:)))

        NdotRayDirection = sum(N(:)*ray%d(:))

        if (abs(NdotRayDirection) < 0.0000001) then
            return
        end if

        d = -sum(N(:)*triangle%v(1, :))
        t = -(sum(N(:)*ray%o(:)) + d)/NdotRayDirection

        if (t < 0) then
            return
        end if

        P = ray%o + t*ray%d

        edge = triangle%v(2, :) - triangle%v(1, :)
        vp = P - triangle%v(1, :)
        C = f_cross(edge, vp)
        if (sum(N(:)*C(:)) < 0) then
            return
        end if

        edge = triangle%v(3, :) - triangle%v(2, :)
        vp = P - triangle%v(2, :)
        C = f_cross(edge, vp)
        if (sum(N(:)*C(:)) < 0) then
            return
        end if

        edge = triangle%v(1, :) - triangle%v(3, :)
        vp = P - triangle%v(3, :)
        C = f_cross(edge, vp)
        if (sum(N(:)*C(:)) < 0) then
            return
        end if

        intersects = .true.

    end function f_intersects_triangle

    subroutine f_distance_normals_3D(model, point, normals, distance)
        type(t_model), INTENT(IN) :: model        
        t_vec3, intent(in) :: point
        t_vec3, intent(out) :: normals
        real(kind(0d0)), intent(out) :: distance

        real(kind(0d0)), dimension(1:model%ntrs, 1:3) :: tri_normals
        real(kind(0d0)), dimension(1:3) :: v1, v2, v3, cross
        real(kind(0d0)), dimension(1:3) :: face_v1, face_v2, face_v3
        real(kind(0d0)) :: xcc, ycc, zcc
        real(kind(0d0)) :: v1_x, v1_y, v1_z, &
                            v2_x, v2_y, v2_z, &
                            v3_x, v3_y, v3_z
        real(kind(0d0)) :: dist_min, dist_buffer, dist_buffer1, dist_buffer2, dist_buffer3
        integer :: i, tri_idx

        dist_min = 1d12
        xcc = point(1); ycc = point(2); zcc = point(3)
        distance = 0d0

        do i = 1, model%ntrs
            v1_x = model%trs(i)%v(1,1)
            v1_y = model%trs(i)%v(1,2)
            v1_z = model%trs(i)%v(1,3)
            dist_buffer1 = dsqrt((xcc-v1_x)**2 + &
                                 (ycc-v1_y)**2 + &
                                 (zcc-v1_z)**2)

            v2_x = model%trs(i)%v(2,1)
            v2_y = model%trs(i)%v(2,2)
            v2_z = model%trs(i)%v(2,3)
            dist_buffer2 = dsqrt((xcc-v2_x)**2 + &
                                (ycc-v2_y)**2 + &
                                (zcc-v2_z)**2)

            v3_x = model%trs(i)%v(3,1)
            v3_y = model%trs(i)%v(3,2)
            v3_z = model%trs(i)%v(3,3)
            dist_buffer3 = dsqrt((xcc-v3_x)**2 + &
                                (ycc-v3_y)**2 + &
                                (zcc-v3_z)**2)

            dist_buffer = MINVAL((/dist_buffer1, dist_buffer2, dist_buffer3/))
    
            if (dist_buffer < dist_min) then
                dist_min = dist_buffer
                tri_idx = i
            end if
        end do

        normals(1:3) = model%trs(tri_idx)%n(1:3)
        distance = dist_min
        
    end subroutine f_distance_normals_3D

    function f_distance(boundary_v, boundary_vertex_count, boundary_edge_count, point, spacing) result(distance)
        integer, intent(in) :: boundary_vertex_count, boundary_edge_count
        real(kind(0d0)), intent(in), dimension(1:boundary_edge_count, 1:3, 1:2) :: boundary_v 
        t_vec3, intent(in) :: point
        t_vec3, intent(in) :: spacing

        integer :: i, j, k
        real(kind(0d0)) :: v_x, v_y, v_z
        real(kind(0d0)) :: xcc, ycc, zcc, &
                         & v1_x, v1_y, v1_z, &
                         & v2_x, v2_y, v2_z, &
                         & v3_x, v3_y, v3_z
        
        real(kind(0d0)) :: dist_buffer1, dist_buffer2, dist_buffer3
        real(kind(0d0)), dimension(1:boundary_edge_count) :: dist_buffer
        real(kind(0d0)) :: distance

        xcc = point(1); ycc = point(2); zcc = point(3)
        distance = 0d0

        do i = 1, boundary_edge_count
            v1_x = boundary_v(i, 1, 1)
            v1_y = boundary_v(i, 1, 2)
            v1_z = 0d0

            dist_buffer1 = dsqrt((xcc-v1_x)**2 + &
                                & (ycc-v1_y)**2)

            v2_x = boundary_v(i, 2, 1)
            v2_y = boundary_v(i, 2, 2)
            v2_z = 0d0

            dist_buffer2 = dsqrt((xcc-v2_x)**2 + &
                                & (ycc-v2_y)**2)

            dist_buffer(i) = MINVAL((/dist_buffer1, dist_buffer2/))
        end do

        distance = MINVAL(dist_buffer)

    end function f_distance

    subroutine f_normals(boundary_v, boundary_vertex_count, boundary_edge_count, point, spacing, normals)
        integer, intent(in) :: boundary_vertex_count, boundary_edge_count
        real(kind(0d0)), intent(in), dimension(1:boundary_edge_count, 1:3, 1:2) :: boundary_v
        t_vec3, intent(in) :: point
        t_vec3, intent(in) :: spacing
        t_vec3, intent(out) :: normals
        integer :: i, j, k, idx_2comp, idx_buffer
        real(kind(0d0)) :: dist_min, dist_buffer, dist_buffer1, dist_buffer2, v_norm
        real(kind(0d0)) :: v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, xcc, ycc, zcc

        xcc = point(1); ycc = point(2); zcc = point(3)
        dist_buffer = 0d0
        dist_min = 1d08
        idx_buffer = 0
        idx_2comp = 0

        do i = 1, boundary_edge_count
            v1_x = boundary_v(i, 1, 1)
            v1_y = boundary_v(i, 1, 2)
            v1_z = 0d0

            dist_buffer1 = dsqrt((xcc-v1_x)**2 + &
                                & (ycc-v1_y)**2)

            v2_x = boundary_v(i, 2, 1)
            v2_y = boundary_v(i, 2, 2)
            v2_z = 0d0

            dist_buffer2 = dsqrt((xcc-v2_x)**2 + &
                                & (ycc-v2_y)**2)

            dist_buffer = MINVAL((/dist_buffer1, dist_buffer2/)) 

            if (dist_buffer1 > dist_buffer2) then
                idx_2comp = 2
            else
                idx_2comp = 1
            end if

            if (dist_buffer < dist_min) then
                dist_min = dist_buffer
                idx_buffer = i
            end if
        end do

        normals(1) = boundary_v(idx_buffer, 3, 1)
        normals(2) = boundary_v(idx_buffer, 3, 2)
        normals(3) = 0d0

    end subroutine f_normals

    ! Find all the boundary edges/vertices of the STL file
    subroutine f_check_boundary(model, boundary_v, boundary_vertex_count, boundary_edge_count)
        type(t_model), INTENT(IN) :: model
        real(kind(0d0)), allocatable, INTENT(OUT), dimension(:, :, :) :: boundary_v  ! Output boundary vertices
        integer, INTENT(OUT) :: boundary_vertex_count, boundary_edge_count
    
        integer :: i, j, edge_count, edge_index, store_index
        real(kind(0d0)), dimension(1:2, 1:2) :: edge
        real(kind(0d0)), dimension(1:2) :: boundary_edge
        real(kind(0d0)), dimension(1:(3*model%ntrs), 1:2, 1:2) :: temp_boundary_v
        integer, dimension(1:(3*model%ntrs)) :: edge_occurrence
        real(kind(0d0)) :: edgetan, initial, v_norm, xnormal, ynormal
    
        ! Total number of edges in 2D STL
        edge_count = 3 * model%ntrs
        
        ! Initialize edge_occurrence array to zero
        edge_occurrence = 0
        edge_index = 0
    
        ! Collect all edges of all triangles and store them
        do i = 1, model%ntrs
            ! First edge (v1, v2)
            edge(1, :) = model%trs(i)%v(1, 1:2)
            edge(2, :) = model%trs(i)%v(2, 1:2)
            call register_edge(temp_boundary_v, edge, edge_index, edge_count)
    
            ! Second edge (v2, v3)
            edge(1, :) = model%trs(i)%v(2, 1:2)
            edge(2, :) = model%trs(i)%v(3, 1:2)
            call register_edge(temp_boundary_v, edge, edge_index, edge_count)
    
            ! Third edge (v3, v1)
            edge(1, :) = model%trs(i)%v(3, 1:2)
            edge(2, :) = model%trs(i)%v(1, 1:2)
            call register_edge(temp_boundary_v, edge, edge_index, edge_count)
        end do
    
        ! Check all edges and count repeated edges
        do i = 1, edge_count
            do j = 1, edge_count
                if (i /= j) then
                    if (((abs(temp_boundary_v(i, 1, 1) - temp_boundary_v(j, 1, 1)) < 1.0d-8) .and. &
                        (abs(temp_boundary_v(i, 1, 2) - temp_boundary_v(j, 1, 2)) < 1.0d-8) .and. &
                        (abs(temp_boundary_v(i, 2, 1) - temp_boundary_v(j, 2, 1)) < 1.0d-8) .and. &
                        (abs(temp_boundary_v(i, 2, 2) - temp_boundary_v(j, 2, 2)) < 1.0d-8)) .or. &
                        ((abs(temp_boundary_v(i, 1, 1) - temp_boundary_v(j, 2, 1)) < 1.0d-8) .and. &
                        (abs(temp_boundary_v(i, 1, 2) - temp_boundary_v(j, 2, 2)) < 1.0d-8) .and. &
                        (abs(temp_boundary_v(i, 2, 1) - temp_boundary_v(j, 1, 1)) < 1.0d-8) .and. &
                        (abs(temp_boundary_v(i, 2, 2) - temp_boundary_v(j, 1, 2)) < 1.0d-8))) then
                        
                        edge_occurrence(i) = edge_occurrence(i) + 1

                    end if
                end if
            end do
        end do
    
        ! Count the number of boundary vertices/edges
        boundary_vertex_count = 0
        boundary_edge_count = 0

        do i = 1, edge_count
            if (edge_occurrence(i) == 0) then
                boundary_vertex_count = boundary_vertex_count + 2
                boundary_edge_count = boundary_edge_count + 1
            end if
        end do
    
        ! Allocate the boundary_v array based on the number of boundary edges
        allocate(boundary_v(boundary_edge_count, 1:3, 1:2))
    
        ! Store boundary vertices
        store_index = 0
        do i = 1, edge_count
            if (edge_occurrence(i) == 0) then
                store_index = store_index + 1
                boundary_v(store_index, 1, :) = temp_boundary_v(i, 1, :)
                boundary_v(store_index, 2, :) = temp_boundary_v(i, 2, :)
            end if
        end do

         do i = 1, boundary_edge_count
            print*, 'edge:',i, 'start', boundary_v(i, 1, 1:2)
            print*, 'edge:',i, 'end', boundary_v(i, 2, 1:2)
            print*, '=============================================='
        end do


        ! Find the normal vector of the boundary edges
        do i = 1, boundary_edge_count
            boundary_edge(1) = boundary_v(i, 2, 1) - boundary_v(i, 1, 1)
            boundary_edge(2) = boundary_v(i, 2, 2) - boundary_v(i, 1, 2)

            edgetan =  boundary_edge(1)/boundary_edge(2)

            if (abs(boundary_edge(2)) < 1d-10) then
                if (edgetan > 0d0) then
                    ynormal = -1
                    xnormal = 0d0
                else
                    ynormal = 1
                    xnormal = 0d0
                end if
            else
                initial = boundary_edge(2)
                ynormal = -edgetan*initial
                xnormal = initial
            end if

            v_norm = dsqrt(xnormal**2 + ynormal**2)

            boundary_v(i, 3, 1) = xnormal/v_norm
            boundary_v(i, 3, 2) = ynormal/v_norm

            ! if ((boundary_edge(1) > 0d0) .and. &
            !     (boundary_edge(2) > 0d0)) then

            !     boundary_v(i, 3, 1) = -abs(boundary_v(i, 3, 1))
            !     boundary_v(i, 3, 2) = abs(boundary_v(i, 3, 2))

            ! else if ((boundary_edge(1) < 0d0) .and. &
            !     (boundary_edge(2) < 0d0)) then

            !     boundary_v(i, 3, 1) = abs(boundary_v(i, 3, 1))
            !     boundary_v(i, 3, 2) = -abs(boundary_v(i, 3, 2))

            ! else if ((boundary_edge(1) < 0d0) .and. &
            !     (boundary_edge(2) > 0d0)) then

            !     boundary_v(i, 3, 1) = -abs(boundary_v(i, 3, 1))
            !     boundary_v(i, 3, 2) = -abs(boundary_v(i, 3, 2))     

            ! else if ((boundary_edge(1) > 0d0) .and. &
            !     (boundary_edge(2) < 0d0)) then

            !     boundary_v(i, 3, 1) = abs(boundary_v(i, 3, 1))
            !     boundary_v(i, 3, 2) = abs(boundary_v(i, 3, 2))  

            ! end if

            ! print*, 'boundary edge', boundary_edge
            ! print*, 'edge:', i, 'normals', boundary_v(i, 3, 1:2)
        end do

    end subroutine f_check_boundary
    
    subroutine register_edge(temp_boundary_v, edge, edge_index, edge_count)
        integer, intent(inout) :: edge_index, edge_count
        real(kind(0d0)), intent(in), dimension(1:2, 1:2) :: edge
        real(kind(0d0)), dimension(1:edge_count, 1:2, 1:2) :: temp_boundary_v
    
        ! Increment edge index and store the edge
        edge_index = edge_index + 1
        temp_boundary_v(edge_index, 1, :) = edge(1, :)
        temp_boundary_v(edge_index, 2, :) = edge(2, :)
    
    end subroutine register_edge

    subroutine f_interpolate_2D(boundary_v, boundary_edge_count, spacing, interpolated_boundary_v, total_vertices)
        real(kind(0d0)), intent(in), dimension(:, :, :) :: boundary_v
        t_vec3, intent(in) :: spacing
        real(kind(0d0)), allocatable, intent(out), dimension(:, :) :: interpolated_boundary_v
        
        integer :: i, j, num_segments, total_vertices, boundary_edge_count
        real(kind(0d0)) :: x1, y1, x2, y2, edge_length, del_x, del_y, cell_width
        real(kind(0d0)), allocatable :: temp_boundary_v(:, :)
      
        ! Get the number of boundary edges

        cell_width = minval(spacing(1:2))
        num_segments = 0
      
        ! First pass: Calculate the total number of vertices including interpolated ones
        total_vertices = 1
        do i = 1, boundary_edge_count
           ! Get the coordinates of the two ends of the current edge
           x1 = boundary_v(i, 1, 1)
           y1 = boundary_v(i, 1, 2)
           x2 = boundary_v(i, 2, 1)
           y2 = boundary_v(i, 2, 2)

        !    print*, 'edge:', i, 'x1', x1 , 'y1', y1, 'x2', x2 , 'y2', y2
      
           ! Compute the length of the edge
           edge_length = dsqrt((x2 - x1)**2 + (y2 - y1)**2)
      
           ! Determine the number of segments
           if (edge_length > cell_width) then
              num_segments = 1000*ceiling(edge_length/cell_width)
           else
              num_segments = 1
           end if
      
           ! Each edge contributes num_segments vertices
           total_vertices = total_vertices + num_segments
        end do
      
        ! Allocate memory for the new boundary vertices array
        allocate(interpolated_boundary_v(1:total_vertices, 1:3))
      
        ! Fill the new boundary vertices array with original and interpolated vertices
        total_vertices = 1
        do i = 1, boundary_edge_count
           ! Get the coordinates of the two ends of the current edge
           x1 = boundary_v(i, 1, 1)
           y1 = boundary_v(i, 1, 2)
           x2 = boundary_v(i, 2, 1)
           y2 = boundary_v(i, 2, 2)
      
           ! Compute the length of the edge
           edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
      
           ! Determine the number of segments and interpolation step
           if (edge_length > cell_width) then
              num_segments = 1000*ceiling(edge_length / cell_width)
              del_x = (x2 - x1) / num_segments
              del_y = (y2 - y1) / num_segments
           else
              num_segments = 1
              del_x = 0d0
              del_y = 0d0
           end if
           
           interpolated_boundary_v(1, 1) = x1
           interpolated_boundary_v(1, 2) = y1
           interpolated_boundary_v(1, 3) = 0d0

           ! Add original and interpolated vertices to the output array
           do j = 1, num_segments - 1
              total_vertices = total_vertices + 1
              interpolated_boundary_v(total_vertices, 1) = x1 + j * del_x
              interpolated_boundary_v(total_vertices, 2) = y1 + j * del_y

            !   print*, i, j, interpolated_boundary_v(total_vertices, 1:2)

           end do
      
           ! Add the last vertex of the edge
           if (num_segments > 0) then
              total_vertices = total_vertices + 1
              interpolated_boundary_v(total_vertices, 1) = x2
              interpolated_boundary_v(total_vertices, 2) = y2
           end if
        end do
      
    end subroutine f_interpolate_2D

    subroutine f_interpolate_3D(model, spacing, interpolated_boundary_v, total_vertices)
        t_vec3, intent(in) :: spacing
        type(t_model), intent(in) :: model
        real(kind(0d0)), allocatable, intent(out), dimension(:, :) :: interpolated_boundary_v
        integer, intent(out) :: total_vertices

        integer :: i, j, k, num_triangles, num_segments
        real(kind(0d0)) :: x1, y1, z1, x2, y2, z2, edge_length, del_x, del_y, del_z, cell_width
        real(kind(0d0)), allocatable :: temp_boundary_v(:, :)

         ! Number of triangles in the model
        num_triangles = size(model%trs)

        cell_width = minval(spacing)

        ! Calculate the total number of vertices including interpolated ones
        total_vertices = 0
        do i = 1, num_triangles
            do j = 1, 3
                ! Get the coordinates of the two vertices of the current edge
                x1 = model%trs(i)%v(j, 1)
                y1 = model%trs(i)%v(j, 2)
                z1 = model%trs(i)%v(j, 3)
                ! Next vertex in the triangle (cyclic)
                x2 = model%trs(i)%v(mod(j, 3) + 1, 1)
                y2 = model%trs(i)%v(mod(j, 3) + 1, 2)
                z2 = model%trs(i)%v(mod(j, 3) + 1, 3)

                ! Compute the length of the edge
                edge_length = dsqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

                ! Determine the number of segments
                if (edge_length > cell_width) then
                    num_segments = 5*ceiling(edge_length / cell_width)
                else
                    num_segments = 1
                end if

                ! Each edge contributes num_segments vertices
                total_vertices = total_vertices + num_segments
            end do
        end do

        ! Allocate memory for the new boundary vertices array
        allocate(interpolated_boundary_v(1:total_vertices, 1:3))

        ! Second pass: Fill the new boundary vertices array with original and interpolated vertices
        total_vertices = 0
        do i = 1, num_triangles
            ! Loop through the 3 edges of each triangle
            do j = 1, 3
                ! Get the coordinates of the two vertices of the current edge
                x1 = model%trs(i)%v(j, 1)
                y1 = model%trs(i)%v(j, 2)
                z1 = model%trs(i)%v(j, 3)
                x2 = model%trs(i)%v(mod(j, 3) + 1, 1)  ! Next vertex in the triangle (cyclic)
                y2 = model%trs(i)%v(mod(j, 3) + 1, 2)
                z2 = model%trs(i)%v(mod(j, 3) + 1, 3)

                ! Compute the length of the edge
                edge_length = dsqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

                ! Determine the number of segments and interpolation step
                if (edge_length > cell_width) then
                    num_segments = 5*ceiling(edge_length / cell_width)
                    del_x = (x2 - x1) / num_segments
                    del_y = (y2 - y1) / num_segments
                    del_z = (z2 - z1) / num_segments
                else
                    num_segments = 1
                    del_x = 0d0
                    del_y = 0d0
                    del_z = 0d0
                end if

                ! Add original and interpolated vertices to the output array
                do k = 0, num_segments - 1
                    total_vertices = total_vertices + 1
                    interpolated_boundary_v(total_vertices, 1) = x1 + k * del_x
                    interpolated_boundary_v(total_vertices, 2) = y1 + k * del_y
                    interpolated_boundary_v(total_vertices, 3) = z1 + k * del_z
                end do

                ! Add the last vertex of the edge
                total_vertices = total_vertices + 1
                interpolated_boundary_v(total_vertices, 1) = x2
                interpolated_boundary_v(total_vertices, 2) = y2
                interpolated_boundary_v(total_vertices, 3) = z2
            end do
        end do

    end subroutine f_interpolate_3D

    function f_interpolated_distance(interpolated_boundary_v, total_vertices, point, spacing) result(distance)
        integer, intent(in) :: total_vertices
        real(kind(0d0)), intent(in), dimension(1:total_vertices, 1:3) :: interpolated_boundary_v 
        t_vec3, intent(in) :: point
        t_vec3, intent(in) :: spacing
        integer :: i
        real(kind(0d0)) :: xcc, ycc, zcc, &
                         & v1_x, v1_y, v1_z
        real(kind(0d0)) :: dist_buffer, min_dist
        real(kind(0d0)) :: distance

        xcc = point(1); ycc = point(2); zcc = point(3)
        distance = 1d12
        dist_buffer = 1d12
        min_dist = 1d12

        do i = 1, total_vertices
            v1_x = interpolated_boundary_v(i, 1)
            v1_y = interpolated_boundary_v(i, 2)
            v1_z = interpolated_boundary_v(i, 3)

            dist_buffer = dsqrt((xcc-v1_x)**2 + &
                                (ycc-v1_y)**2 + &
                                (zcc-v1_z)**2)

            if (min_dist > dist_buffer) then
                min_dist = dist_buffer
            end if

        end do

        distance = min_dist

        ! if (distance < 9d-3) then
        !     distance = 0d0
        ! end if

    end function f_interpolated_distance

    subroutine f_check_interpolation_2D(boundary_v, boundary_edge_count, spacing, interpolate)
        logical, intent(out) :: interpolate
        integer, INTENT(IN) ::  boundary_edge_count
        real(kind(0d0)), optional, intent(in), dimension(1:boundary_edge_count, 1:3, 1:2) :: boundary_v
        t_vec3, intent(in) :: spacing

        real(kind(0d0)) :: l1, cell_width
        real(kind(0d0)) :: v1_x, v1_y, v2_x, v2_y
        integer :: j
        
        cell_width = minval(spacing(1:2))
        interpolate = .false.

        do j = 1, boundary_edge_count
            v1_x = boundary_v(j, 1, 1)
            v2_x = boundary_v(j, 2, 1)
            v1_y = boundary_v(j, 1, 2)
            v2_y = boundary_v(j, 2, 2)

            l1 = dsqrt((v2_x - v1_x)**2 + (v2_y - v1_y)**2)

            if ((l1 > cell_width)) then
                interpolate = .true.
            end if
        end do

    end subroutine f_check_interpolation_2D

    subroutine f_check_interpolation_3D(model, spacing, interpolate)
        logical, intent(out) :: interpolate
        type(t_model), INTENT(IN) :: model
        t_vec3, intent(in) :: spacing

        real(kind(0d0)) :: l1, l2, l3, cell_width
        real(kind(0d0)) :: v1_x, v1_y, v1_z, &
                        & v2_x, v2_y, v2_z, &
                        & v3_x, v3_y, v3_z
        integer :: i
        
        cell_width = minval(spacing)
        interpolate = .false.

        do i = 1, model%ntrs
            v1_x = model%trs(i)%v(1, 1)
            v2_x = model%trs(i)%v(2, 1)
            v3_x = model%trs(i)%v(3, 1)

            v1_y = model%trs(i)%v(1, 2)
            v2_y = model%trs(i)%v(2, 2)
            v3_y = model%trs(i)%v(3, 2)

            v1_z = model%trs(i)%v(1, 3)
            v2_z = model%trs(i)%v(2, 3)
            v3_z = model%trs(i)%v(3, 3)

            l1 = dsqrt((v2_x - v1_x)**2 + (v2_y - v1_y)**2 + (v2_z - v1_z)**2)
            l2 = dsqrt((v3_x - v2_x)**2 + (v3_y - v2_y)**2 + (v3_z - v2_z)**2)
            l3 = dsqrt((v2_x - v1_x)**2 + (v2_y - v1_y)**2 + (v2_z - v1_z)**2)

            if ((l1 > cell_width) .or. (l2 > cell_width) .or. (l3 > cell_width)) then
                interpolate = .true.
            end if
        end do

    end subroutine f_check_interpolation_3D


      

end module m_model
