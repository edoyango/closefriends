! Module to perform direct search for pairs of points within a fixed cutoff.
! cellList subroutine
!   dim:      An integer defining dimension of the coordinates.
!   npoints:  An integer with number of points.
!   maxpercell: An integer with max number of points in a cell.
!   x:        An real/double precision array of dim x points dimension containing the list of points to find pairs of.
!   cutoff:   The distance (same type as x) which defines a pair
!   maxnpair: An integer of the maximum number of pairs expected.
!   npairs:   An integer with the number of pairs found.
!   pair_i:   An integer array of maxnpairs length with the "left-sided" point in a pair.
!   pair_j:   An integer array of maxnpairs length with the "right-sided" point in a pair.

module closeFriends
 
    public
    ! private:: sort_hashes, coordsToHash, hashIncrement, rearrange, find_starts, update_pair_list
 
    interface sort_hashes
 
       subroutine sort_hashes(n, hashes, idx) bind(C)
 
          use iso_c_binding
 
          integer(c_int), value, intent(in):: n
          integer(c_int), intent(in):: hashes(*)
          integer(c_int), intent(out):: idx(*)
 
       end subroutine sort_hashes
 
    end interface sort_hashes
 
 contains
 
    !---------------------------------------------------------------------------
    function coordsToHash(n, dim, x, minx, ngridx, cutoff) result(gridhash)
       ! converts raw 2/3D positions to gridhashes based on the grid dimensions
       ! described by minx (lowest grid coordinate), and ngridx (number of grid
       ! cells in xyz direction).
 
       implicit none
       integer, intent(in):: n, dim, ngridx(3)
       double precision, intent(in):: x(dim, n), minx(dim), cutoff
       integer:: i, icell(3), gridhash(n)
 
       if (dim == 2) icell(3) = 1
 
       do i = 1, n
          icell(1:dim) = int((x(1:dim, i) - minx(1:dim))/cutoff) + 1
          gridhash(i) = ngridx(1)*ngridx(2)*(icell(3) - 1) + ngridx(1)*(icell(2) - 1) + icell(1)
       end do
 
    end function coordsToHash
 
    !---------------------------------------------------------------------------
    function hashIncrement(ngridx, dcellx, dcelly, dcellz) result(dhash)
       ! calculates the change in hash, given a change in xyz cell indices
 
       implicit none
       integer, intent(in):: ngridx(3), dcellx, dcelly, dcellz
       integer:: dhash
 
       dhash = ngridx(1)*ngridx(2)*dcellx + ngridx(1)*dcelly + dcellz
 
    end function hashIncrement
 
    !---------------------------------------------------------------------------
    subroutine rearrange(n, dim, idx, gridhash, x, undo)
       ! rearranges gridhash and x, using the indices obtained from a sort by
       ! key.
 
       implicit none
       logical, intent(in):: undo
       integer, intent(in):: n, dim, idx(n)
       integer, intent(inout):: gridhash(n)
       double precision, intent(inout):: x(dim, n)
       integer, allocatable:: tmpgridhash(:)
       double precision, allocatable:: tmpx(:, :)
       integer:: i
 
       allocate (tmpgridhash(n), tmpx(dim, n))
 
       if (undo) then
 
          do i = 1, n
             tmpgridhash(idx(i)) = gridhash(i)
             tmpx(:, idx(i)) = x(:, i)
          end do
 
       else
 
          do i = 1, n
             tmpgridhash(i) = gridhash(idx(i))
             tmpx(:, i) = x(:, idx(i))
          end do
       
       end if
 
       gridhash(:) = tmpgridhash(:)
       x(:, :) = tmpx(:, :)
 
    end subroutine rearrange
 
    !---------------------------------------------------------------------------
    subroutine find_starts(n, ngridx, gridhash, starts)
       ! finds the starting index for each grid cell. See "Building the grid
       ! using Sorting" section in https://developer.download.nvidia.com/assets/cuda/files/particles.pdf
 
       implicit none
       integer, intent(in):: n, ngridx(3)
       integer, intent(in):: gridhash(n)
       integer, intent(out):: starts(:)
       integer:: i
 
       starts(1:gridhash(1)) = 1
 
       do i = 2, n
          if (gridhash(i) /= gridhash(i - 1)) starts(gridhash(i - 1) + 1:gridhash(i)) = i
       end do
 
       starts(gridhash(n) + 1:ngridx(1)*ngridx(2)*ngridx(3)) = n + 1
 
    end subroutine find_starts
 
    !---------------------------------------------------------------------------
    subroutine update_pair_list(dim, n, maxnpair, cutoff, jstart, jend, i, x, npairs, pair_i, pair_j)
       ! subroutine to loop over sequential grid cells and check if j point is
       ! within a given i point
 
       implicit none
       integer, intent(in):: jstart, jend, i, n, maxnpair, dim
       double precision, intent(in):: cutoff, x(dim, n)
       integer, intent(inout):: npairs, pair_i(maxnpair), pair_j(maxnpair)
    !    integer, intent(inout):: npairs, pair_i(0:), pair_j(0:)
       integer:: j
       double precision:: r2
 
       do j = jstart, jend
          r2 = sum((x(:, i) - x(:, j))**2)
          if (r2 <= cutoff*cutoff) then
             npairs = npairs + 1
             pair_i(npairs) = i
             pair_j(npairs) = j
            !  npairs = npairs + merge(1, 0, r2 <= cutoff*cutoff)
          end if
       end do
 
    end subroutine update_pair_list
 
    !---------------------------------------------------------------------------
    subroutine cellList(dim, npoints, x, cutoff, maxnpair, npairs, pair_i, pair_j)
       ! main subroutine to perform pair search
 
       implicit none
       integer, intent(in):: dim, npoints, maxnpair
       double precision, intent(in):: cutoff
       double precision, intent(inout):: x(dim, npoints)
       integer, intent(out)::  npairs, pair_i(maxnpair), pair_j(maxnpair)
       integer:: i, hashi, hashj, ngridx(3)
       double precision:: minx(3), maxx(3)
       integer, allocatable:: gridhash(:), idx(:), starts(:)
 
       ! determine grid min-max coordinates, with concessions made for dim=2 case
       minx(1:dim) = minval(x, dim=2)
       if (dim == 2) minx(3) = 0.d0
       maxx(1:dim) = maxval(x, dim=2)
       if (dim == 2) maxx(3) = 0.d0
       
       ! creating buffer layers
       minx(:) = minx(:) - 2.d0*cutoff
       maxx(:) = maxx(:) + 2.d0*cutoff
 
       ! determining no. of grid cells and adjusting maximum extent
       ngridx(:) = int((maxx(:) - minx(:))/cutoff) + 1
       maxx(:) = maxx(:) + ngridx(:)*cutoff
 
       ! convert coordinates to cell grid hashes
       allocate (gridhash(npoints))
       gridhash = coordsToHash(npoints, dim, x, minx, ngridx, cutoff)
 
       ! get sorting indices to reorder points by ascending grid hash
       allocate (idx(npoints))
       call sort_hashes(npoints, gridhash, idx)
 
       ! reorder points (x, gridhash) by ascending grid hash
       call rearrange(npoints, dim, idx, gridhash, x, undo=.false.)
 
       ! find starting points of each grid cell
       allocate (starts(product(ngridx)))
       call find_starts(npoints, ngridx, gridhash, starts)
 
       ! perform pair search
       npairs = 0
       do hashi = gridhash(1), gridhash(npoints) ! for all relevant grid cell points
          do i = starts(hashi), starts(hashi + 1) - 1 ! loop over all points in the cell
             ! loop over "other" points in cell + next cell (1st adjacent tell: top of hashi cell)
             call update_pair_list(dim, npoints, maxnpair, cutoff, i + 1, starts(hashi + 2) - 1, i, x, npairs, pair_i, &
                pair_j)
             ! loop over 2nd-4rd adjacent cells (bottom-south-west to top-south-west)
             hashj = hashi + hashIncrement(ngridx, -1, -1, -1)
             call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x, npairs, &
                pair_i, pair_j)
             ! loop over 5th-7th adjacent cells (bottom-centre-west to top-centre-west)
             hashj = hashi + hashIncrement(ngridx, -1, 0, -1)
             call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x, npairs, &
                pair_i, pair_j)
             ! loop over 8th-10th adjacent cells (bottom-north-west to top-north-west)
             hashj = hashi + hashIncrement(ngridx, -1, 1, -1)
             call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x, npairs, &
                pair_i, pair_j)
             ! loop over 11th-13th adjacent cells (bottom-south-centre to top-south-centre)
             hashj = hashi + hashIncrement(ngridx, 0, -1, -1)
             call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x, npairs, &
                pair_i, pair_j)
          end do
       end do
 
       ! call rearrange(npoints, dim, idx, gridhash, x, undo=.true.)
 
       ! do i = 1, npairs
       !    pair_i(i) = idx(pair_i(i))
       !    pair_j(i) = idx(pair_j(i))
       ! end do
 
    end subroutine cellList
 
 end module closeFriends