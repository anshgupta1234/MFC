!>
!! @file m_compute_levelset.fpp
!! @brief Contains module m_compute_levelset

#:include 'macros.fpp'

!> @brief This module is used to handle all operations related to immersed
!!              boundary methods (IBMs)
module m_compute_levelset

    ! Dependencies =============================================================

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

    private; public :: s_compute_cylinder_levelset, s_compute_circle_levelset, &
 s_compute_airfoil_levelset, &
 s_compute_3D_airfoil_levelset, &
 s_compute_rectangle_levelset, &
 s_compute_sphere_levelset, &
 s_compute_2D_STL_levelset, &
 s_compute_3D_STL_levelset

    real(kind(0d0)) :: x_centroid, y_centroid, z_centroid
    real(kind(0d0)) :: length_x, length_y, length_z
    real(kind(0d0)) :: radius

    type(bounds_info) :: x_boundary, y_boundary, z_boundary  !<
        !! These variables combine the centroid and length parameters associated with
        !! a particular patch to yield the locations of the patch boundaries in the
        !! x-, y- and z-coordinate directions. They are used as a means to concisely
        !! perform the actions necessary to lay out a particular patch on the grid.

contains

    !>  Initialize IBM module
    subroutine s_compute_2D_STL_levelset(levelset, levelset_norm, ib_patch_id, ghost_points, num_gps, ib_markers)
        integer, intent(INOUT) :: num_gps
        type(ghost_point), dimension(num_gps), intent(INOUT) :: ghost_points
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        real(kind(0d0)), dimension(3) :: dist_vec, dist_vec_buffer1, dist_vec_buffer2, dist_vec_buffer3
        integer, intent(IN) :: ib_patch_id
        real(kind(0d0)), dimension(0:num_gps) :: distance
        integer :: i, j, ii, jj, q, ii_buffer, jj_buffer, a !< Loop index variables
        integer :: ii_buffer_avg1, jj_buffer_avg1, ii_buffer_avg2, jj_buffer_avg2
        type(ghost_point) :: gp
        real(kind(0d0)) :: distance_q, distance_buffer
        type(integer_field), intent(INOUT) :: ib_markers

        do i = 0, m
            do j = 0, n
                distance_buffer = 1d12
                ii_buffer = 0
                jj_buffer = 0
                ii_buffer_avg1 = -1000
                jj_buffer_avg1 = -1000
                ii_buffer_avg2 = -1000
                jj_buffer_avg2 = -1000
                distance = 0d0

                do q = 1, num_gps
                    gp = ghost_points(q)
                    ii = gp%loc(1)
                    jj = gp%loc(2)

                    if (gp%IBB > 0) then
                        distance(q) = dsqrt((x_cc(i) - x_cc(ii))**2 &
                                    + (y_cc(j) - y_cc(jj))**2)
                    else
                        distance(q) = 1d12
                    end if

                    if (abs(distance_buffer - distance(q)) < 1d-12) then
                        if (ii_buffer_avg1 == -1000) then
                            ii_buffer_avg1 = ii
                            jj_buffer_avg1 = jj
                        else if (ii_buffer_avg1 /= -1000) then
                            ii_buffer_avg2 = ii
                            jj_buffer_avg2 = jj
                        end if
                    else if (distance_buffer > distance(q)) then
                            distance_buffer = distance(q)
                            ii_buffer = ii
                            jj_buffer = jj
                    end if
                end do

                dist_vec_buffer1(1) = x_cc(i) - x_cc(ii_buffer)
                dist_vec_buffer1(2) = y_cc(j) - y_cc(jj_buffer)
                dist_vec_buffer1(3) = 0

                dist_vec_buffer2(1) = x_cc(i) - x_cc(ii_buffer_avg1)
                dist_vec_buffer2(2) = y_cc(j) - y_cc(jj_buffer_avg1)
                dist_vec_buffer2(3) = 0

                dist_vec_buffer3(1) = x_cc(i) - x_cc(ii_buffer_avg2)
                dist_vec_buffer3(2) = y_cc(j) - y_cc(jj_buffer_avg2)
                dist_vec_buffer3(3) = 0


                if (ii_buffer_avg2 /= -1000) then
                    dist_vec(1) = (dist_vec_buffer1(1) + dist_vec_buffer2(1) + dist_vec_buffer3(1))
                    dist_vec(2) = (dist_vec_buffer1(2) + dist_vec_buffer2(2) + dist_vec_buffer3(2))
                    dist_vec(3) = 0
                    distance_q = dsqrt(dist_vec(1)**2 + dist_vec(2)**2 + dist_vec(3)**2)
                else if (ii_buffer_avg1 /= -1000) then
                    dist_vec(1) = (dist_vec_buffer1(1) + dist_vec_buffer2(1))
                    dist_vec(2) = (dist_vec_buffer1(2) + dist_vec_buffer2(2))
                    dist_vec(3) = 0
                    distance_q = dsqrt(dist_vec(1)**2 + dist_vec(2)**2 + dist_vec(3)**2)
                else
                    dist_vec(1) = x_cc(i) - x_cc(ii_buffer)
                    dist_vec(2) = y_cc(j) - y_cc(jj_buffer)
                    dist_vec(3) = 0
                    distance_q = distance_buffer
                end if

                if (ib_markers%sf(i, j, 0) /= 0) then
                    distance_q = -distance_q
                else
                    distance_q = abs(distance_q)
                end if

                levelset(i, j, 0, ib_patch_id) = distance_q

                if (distance_q == 0) then
                    levelset_norm(i, j, 0, ib_patch_id, :) = 0
                else
                    levelset_norm(i, j, 0, ib_patch_id, :) = &
                        dist_vec(:)/abs(distance_q)
                end if
            end do
        end do

        ! do i = 2, m-2
        !     do j = 0, n
        !         do a = 1, 3
        !             levelset_norm(i, j, 0, ib_patch_id, a) = levelset_norm(i-2, j, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i-1, j, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i, j, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i+1, j, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i+2, j, 0, ib_patch_id, a)
        !         end do
        !     end do
        ! end do

        ! do j = 2, n-2
        !     do i = 0, m
        !         do a = 1, 3
        !             levelset_norm(i, j, 0, ib_patch_id, a) = levelset_norm(i, j-2, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i, j-1, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i, j, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i, j+1, 0, ib_patch_id, a) &
        !                                                         + levelset_norm(i, j+2, 0, ib_patch_id, a)
        !         end do
        !         distance_buffer = dsqrt(levelset_norm(i, j, 0, ib_patch_id, 1)**2 + levelset_norm(i, j, 0, ib_patch_id, 2)**2)
        !         levelset_norm(i, j, 0, ib_patch_id, :) = levelset_norm(i, j, 0, ib_patch_id, :)/distance_buffer
        !     end do
        ! end do


    end subroutine s_compute_2D_STL_levelset

    subroutine s_compute_3D_STL_levelset(levelset, levelset_norm, ib_patch_id, ghost_points, num_gps, ib_markers)
        integer, intent(INOUT) :: num_gps
        type(ghost_point), dimension(num_gps), intent(INOUT) :: ghost_points
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        real(kind(0d0)), dimension(3) :: dist_vec
        integer, intent(IN) :: ib_patch_id
        real(kind(0d0)), dimension(0:num_gps) :: distance
        integer :: i, j, k, ii, jj, kk, q, ii_buffer, jj_buffer, kk_buffer !< Loop index variables
        type(ghost_point) :: gp
        real(kind(0d0)) :: distance_buffer
        type(integer_field), intent(INOUT) :: ib_markers

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    distance_buffer = 1d12
                    ii_buffer = 0
                    jj_buffer = 0
                    kk_buffer = 0
                    distance = 0d0

                    do q = 1, num_gps
                        gp = ghost_points(q)
                        ii = gp%loc(1)
                        jj = gp%loc(2)
                        kk = gp%loc(3)

                        if (gp%IBB > 0) then
                            distance(q) = dsqrt((x_cc(i) - x_cc(ii))**2 &
                                        + (y_cc(j) - y_cc(jj))**2 & 
                                        + (z_cc(k) - z_cc(kk))**2)
                        else
                            distance(q) = 1d12
                        end if

                        if (distance_buffer > distance(q)) then
                            distance_buffer = distance(q)
                            ii_buffer = ii
                            jj_buffer = jj
                            kk_buffer = kk
                        end if
                    end do

                    if (ib_markers%sf(i, j, k) /= 0) then
                        distance_buffer = -distance_buffer
                    else
                        distance_buffer = abs(distance_buffer)
                    end if

                    dist_vec(1) = x_cc(i) - x_cc(ii_buffer)
                    dist_vec(2) = y_cc(j) - y_cc(jj_buffer)
                    dist_vec(3) = z_cc(k) - z_cc(kk_buffer)

                    levelset(i, j, k, ib_patch_id) = distance_buffer

                    if (distance_buffer == 0) then
                        levelset_norm(i, j, k, ib_patch_id, :) = 0
                    else
                        levelset_norm(i, j, k, ib_patch_id, :) = &
                            dist_vec(:)/abs(distance_buffer)
                    end if
                end do
            end do
        end do

    end subroutine s_compute_3D_STL_levelset

    subroutine s_compute_circle_levelset(levelset, levelset_norm, ib_patch_id)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(kind(0d0)) :: radius, dist
        real(kind(0d0)) :: x_centroid, y_centroid
        real(kind(0d0)), dimension(3) :: dist_vec

        integer :: i, j !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid

        do i = 0, m
            do j = 0, n

                dist_vec(1) = x_cc(i) - x_centroid
                dist_vec(2) = y_cc(j) - y_centroid
                dist_vec(3) = 0
                dist = dsqrt(sum(dist_vec**2))
                levelset(i, j, 0, ib_patch_id) = dist - radius
                if (dist == 0) then
                    levelset_norm(i, j, 0, ib_patch_id, :) = 0
                else
                    levelset_norm(i, j, 0, ib_patch_id, :) = &
                        dist_vec(:)/dist
                end if

            end do
        end do

    end subroutine s_compute_circle_levelset

    subroutine s_compute_airfoil_levelset(levelset, levelset_norm, ib_patch_id)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(kind(0d0)) :: radius, dist, global_dist
        integer :: global_id
        real(kind(0d0)) :: x_centroid, y_centroid, x_act, y_act, theta
        real(kind(0d0)), dimension(3) :: dist_vec

        integer :: i, j, k !< Loop index variables

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        theta = pi*patch_ib(ib_patch_id)%theta/180d0

        do i = 0, m
            do j = 0, n

                if (patch_ib(ib_patch_id)%theta /= dflt_real) then
                    x_act = (x_cc(i) - x_centroid)*cos(theta) - (y_cc(j) - y_centroid)*sin(theta) + x_centroid
                    y_act = (x_cc(i) - x_centroid)*sin(theta) + (y_cc(j) - y_centroid)*cos(theta) + y_centroid
                else
                    x_act = x_cc(i)
                    y_act = y_cc(j)
                end if

                if (y_act >= y_centroid) then
                    do k = 1, Np
                        dist_vec(1) = x_cc(i) - airfoil_grid_u(k)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_u(k)%y
                        dist_vec(3) = 0
                        dist = dsqrt(sum(dist_vec**2))
                        if (k == 1) then
                            global_dist = dist
                            global_id = k
                        else
                            if (dist < global_dist) then
                                global_dist = dist
                                global_id = k
                            end if
                        end if
                    end do
                    dist_vec(1) = x_cc(i) - airfoil_grid_u(global_id)%x
                    dist_vec(2) = y_cc(j) - airfoil_grid_u(global_id)%y
                    dist_vec(3) = 0
                    dist = global_dist
                else
                    do k = 1, Np
                        dist_vec(1) = x_cc(i) - airfoil_grid_l(k)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_l(k)%y
                        dist_vec(3) = 0
                        dist = dsqrt(sum(dist_vec**2))
                        if (k == 1) then
                            global_dist = dist
                            global_id = k
                        else
                            if (dist < global_dist) then
                                global_dist = dist
                                global_id = k
                            end if
                        end if
                    end do
                    dist_vec(1) = x_cc(i) - airfoil_grid_l(global_id)%x
                    dist_vec(2) = y_cc(j) - airfoil_grid_l(global_id)%y
                    dist_vec(3) = 0
                    dist = global_dist
                end if

                levelset(i, j, 0, ib_patch_id) = dist
                if (dist == 0) then
                    levelset_norm(i, j, 0, ib_patch_id, :) = 0
                else
                    levelset_norm(i, j, 0, ib_patch_id, :) = &
                        dist_vec(:)/dist
                end if

            end do
        end do

    end subroutine s_compute_airfoil_levelset

    subroutine s_compute_3D_airfoil_levelset(levelset, levelset_norm, ib_patch_id)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(kind(0d0)) :: radius, dist, dist_surf, dist_side, global_dist
        integer :: global_id
        real(kind(0d0)) :: x_centroid, y_centroid, z_centroid, lz, z_max, z_min, x_act, y_act, theta
        real(kind(0d0)), dimension(3) :: dist_vec

        integer :: i, j, k, l !< Loop index variables

        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid
        lz = patch_ib(ib_patch_id)%length_z
        theta = pi*patch_ib(ib_patch_id)%theta/180d0

        z_max = z_centroid + lz/2
        z_min = z_centroid - lz/2

        do l = 0, p
            do j = 0, n
                do i = 0, m

                    if (patch_ib(ib_patch_id)%theta /= dflt_real) then
                        x_act = (x_cc(i) - x_centroid)*cos(theta) - (y_cc(j) - y_centroid)*sin(theta) + x_centroid
                        y_act = (x_cc(i) - x_centroid)*sin(theta) + (y_cc(j) - y_centroid)*cos(theta) + y_centroid
                    else
                        x_act = x_cc(i)
                        y_act = y_cc(j)
                    end if

                    if (y_act >= y_centroid) then
                        do k = 1, Np
                            dist_vec(1) = x_cc(i) - airfoil_grid_u(k)%x
                            dist_vec(2) = y_cc(j) - airfoil_grid_u(k)%y
                            dist_vec(3) = 0
                            dist_surf = dsqrt(sum(dist_vec**2))
                            if (k == 1) then
                                global_dist = dist_surf
                                global_id = k
                            else
                                if (dist_surf < global_dist) then
                                    global_dist = dist_surf
                                    global_id = k
                                end if
                            end if
                        end do
                        dist_vec(1) = x_cc(i) - airfoil_grid_u(global_id)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_u(global_id)%y
                        dist_vec(3) = 0
                        dist_surf = global_dist
                    else
                        do k = 1, Np
                            dist_vec(1) = x_cc(i) - airfoil_grid_l(k)%x
                            dist_vec(2) = y_cc(j) - airfoil_grid_l(k)%y
                            dist_vec(3) = 0
                            dist_surf = dsqrt(sum(dist_vec**2))
                            if (k == 1) then
                                global_dist = dist_surf
                                global_id = k
                            else
                                if (dist_surf < global_dist) then
                                    global_dist = dist_surf
                                    global_id = k
                                end if
                            end if
                        end do
                        dist_vec(1) = x_cc(i) - airfoil_grid_l(global_id)%x
                        dist_vec(2) = y_cc(j) - airfoil_grid_l(global_id)%y
                        dist_vec(3) = 0
                        dist_surf = global_dist
                    end if

                    dist_side = min(abs(z_cc(l) - z_min), abs(z_max - z_cc(l)))

                    if (dist_side < dist_surf) then
                        levelset(i, j, l, ib_patch_id) = dist_side
                        if (dist_side == abs(z_cc(l) - z_min)) then
                            levelset_norm(i, j, l, ib_patch_id, :) = (/0, 0, -1/)
                        else
                            levelset_norm(i, j, l, ib_patch_id, :) = (/0, 0, 1/)
                        end if
                    else
                        levelset(i, j, l, ib_patch_id) = dist_surf
                        if (dist == 0) then
                            levelset_norm(i, j, l, ib_patch_id, :) = 0
                        else
                            levelset_norm(i, j, l, ib_patch_id, :) = &
                                dist_vec(:)/dist_surf
                        end if
                    end if

                end do
            end do
        end do

    end subroutine s_compute_3D_airfoil_levelset

    !>  Initialize IBM module
    subroutine s_compute_rectangle_levelset(levelset, levelset_norm, ib_patch_id)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm

        integer :: ib_patch_id
        real(kind(0d0)) :: top_right(2), bottom_left(2)
        real(kind(0d0)) :: x, y, min_dist
        real(kind(0d0)) :: side_dists(4)

        integer :: i, j, k !< Loop index variables

        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid

        top_right(1) = x_centroid + length_x/2
        top_right(2) = y_centroid + length_y/2

        bottom_left(1) = x_centroid - length_x/2
        bottom_left(2) = y_centroid - length_y/2

        do i = 0, m
            do j = 0, n

                x = x_cc(i)
                y = y_cc(j)

                if ((x > bottom_left(1) .and. x < top_right(1)) .or. &
                    (y > bottom_left(2) .and. y < top_right(2))) then

                    side_dists(1) = bottom_left(1) - x
                    side_dists(2) = x - top_right(1)
                    side_dists(3) = bottom_left(2) - y
                    side_dists(4) = y - top_right(2)

                    min_dist = minval(abs(side_dists))

                    if (i == 50 .and. j == 0) then
                    end if

                    if (min_dist == abs(side_dists(1))) then
                        levelset(i, j, 0, ib_patch_id) = side_dists(1)
                        if (side_dists(1) == 0) then
                            levelset_norm(i, j, 0, ib_patch_id, 1) = 0d0
                        else
                            levelset_norm(i, j, 0, ib_patch_id, 1) = side_dists(1)/ &
                                                                     abs(side_dists(1))
                        end if

                    else if (min_dist == abs(side_dists(2))) then
                        levelset(i, j, 0, ib_patch_id) = side_dists(2)
                        if (side_dists(2) == 0) then
                            levelset_norm(i, j, 0, ib_patch_id, 1) = 0d0
                        else
                            levelset_norm(i, j, 0, ib_patch_id, 1) = side_dists(2)/ &
                                                                     abs(side_dists(2))
                        end if

                    else if (min_dist == abs(side_dists(3))) then
                        if (side_dists(3) == 0) then
                            levelset_norm(i, j, 0, ib_patch_id, 1) = 0d0
                        else
                            levelset_norm(i, j, 0, ib_patch_id, 1) = side_dists(3)/ &
                                                                     abs(side_dists(3))
                        end if

                    else if (min_dist == abs(side_dists(4))) then
                        if (side_dists(4) == 0) then
                            levelset_norm(i, j, 0, ib_patch_id, 1) = 0d0
                        else
                            levelset_norm(i, j, 0, ib_patch_id, 1) = side_dists(4)/ &
                                                                     abs(side_dists(4))
                        end if

                    end if

                end if

            end do
        end do

    end subroutine s_compute_rectangle_levelset

    subroutine s_compute_sphere_levelset(levelset, levelset_norm, ib_patch_id)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(kind(0d0)) :: radius, dist
        real(kind(0d0)) :: x_centroid, y_centroid, z_centroid
        real(kind(0d0)), dimension(3) :: dist_vec

        integer :: i, j, k !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    dist_vec(1) = x_cc(i) - x_centroid
                    dist_vec(2) = y_cc(j) - y_centroid
                    dist_vec(3) = z_cc(k) - z_centroid
                    dist = dsqrt(sum(dist_vec**2))
                    levelset(i, j, k, ib_patch_id) = dist - radius
                    if (dist == 0) then
                        levelset_norm(i, j, k, ib_patch_id, :) = (/1, 0, 0/)
                    else
                        levelset_norm(i, j, k, ib_patch_id, :) = &
                            dist_vec(:)/dist
                    end if
                end do
            end do
        end do

    end subroutine s_compute_sphere_levelset

    subroutine s_compute_cylinder_levelset(levelset, levelset_norm, ib_patch_id)

        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs), intent(INOUT) :: levelset
        real(kind(0d0)), dimension(0:m, 0:n, 0:p, num_ibs, 3), intent(INOUT) :: levelset_norm
        integer, intent(IN) :: ib_patch_id

        real(kind(0d0)) :: radius, dist
        real(kind(0d0)) :: x_centroid, y_centroid, z_centroid
        real(kind(0d0)) :: length_x, length_y, length_z
        real(kind(0d0)), dimension(3) :: pos_vec, centroid_vec, dist_vec, dist_sides_vec, dist_surface_vec
        real(kind(0d0)) :: dist_side, dist_surface, side_pos
        type(bounds_info) :: boundary
        integer :: i, j, k !< Loop index variables

        radius = patch_ib(ib_patch_id)%radius
        x_centroid = patch_ib(ib_patch_id)%x_centroid
        y_centroid = patch_ib(ib_patch_id)%y_centroid
        z_centroid = patch_ib(ib_patch_id)%z_centroid
        length_x = patch_ib(ib_patch_id)%length_x
        length_y = patch_ib(ib_patch_id)%length_y
        length_z = patch_ib(ib_patch_id)%length_z

        if (length_x /= 0d0) then
            boundary%beg = x_centroid - 0.5*length_x
            boundary%end = x_centroid + 0.5*length_x
            dist_sides_vec = (/1, 0, 0/)
            dist_surface_vec = (/0, 1, 1/)
        else if (length_y /= 0d0) then
            boundary%beg = y_centroid - 0.5*length_y
            boundary%end = y_centroid + 0.5*length_y
            dist_sides_vec = (/0, 1, 0/)
            dist_surface_vec = (/1, 0, 1/)
        else if (length_z /= 0d0) then
            boundary%beg = z_centroid - 0.5*length_z
            boundary%end = z_centroid + 0.5*length_z
            dist_sides_vec = (/0, 0, 1/)
            dist_surface_vec = (/1, 1, 0/)
        end if

        do i = 0, m
            do j = 0, n
                do k = 0, p
                    pos_vec = [x_cc(i), y_cc(j), z_cc(k)]
                    centroid_vec = [x_centroid, y_centroid, z_centroid]
                    side_pos = dot_product(pos_vec, dist_sides_vec)
                    dist_side = min(abs(side_pos - boundary%beg), &
                                    abs(boundary%end - side_pos))

                    dist_surface = norm2((pos_vec - centroid_vec)*dist_surface_vec) &
                                   - radius

                    if (dist_side < abs(dist_surface)) then
                        levelset(i, j, k, ib_patch_id) = -dist_side
                        if (dist_side == abs(side_pos - boundary%beg)) then
                            levelset_norm(i, j, k, ib_patch_id, :) = -dist_sides_vec
                        else
                            levelset_norm(i, j, k, ib_patch_id, :) = dist_sides_vec
                        end if
                    else
                        levelset(i, j, k, ib_patch_id) = dist_surface

                        levelset_norm(i, j, k, ib_patch_id, :) = &
                            (pos_vec - centroid_vec)*dist_surface_vec
                        levelset_norm(i, j, k, ib_patch_id, :) = &
                            levelset_norm(i, j, k, ib_patch_id, :)/ &
                            norm2(levelset_norm(i, j, k, ib_patch_id, :))
                    end if
                end do
            end do
        end do

    end subroutine s_compute_cylinder_levelset

end module m_compute_levelset
