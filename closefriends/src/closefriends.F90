module closefriends

   implicit none

   !----------------------------------------------------------------------------
   interface sort_hashes
      subroutine sort_hashes(n, hashes, idx) bind(C)
         ! performs argsort using C++ sort function
         use iso_c_binding
         integer(c_int), value, intent(in):: n
         integer(c_int), intent(in):: hashes(*)
         integer(c_int), intent(out):: idx(*)
      end subroutine sort_hashes
   end interface sort_hashes

contains

   !----------------------------------------------------------------------------
   subroutine rearrange(n, dim, idx, gridhash, x)
      ! rearranges gridhash and x, using the indices obtained from a sort by
      ! key.

      implicit none
      integer, intent(in):: n, dim, idx(n)
      integer, intent(inout):: gridhash(n)
      double precision, intent(inout):: x(dim, n)
      integer, allocatable:: tmpgridhash(:)
      double precision, allocatable:: tmpx(:, :)
      integer:: i

      allocate (tmpgridhash(n), tmpx(dim, n))

      do i = 1, n
         tmpgridhash(i) = gridhash(idx(i))
         tmpx(:, i) = x(:, idx(i))
      end do

      gridhash(:) = tmpgridhash(:)
      x(:, :) = tmpx(:, :)

   end subroutine rearrange

   !----------------------------------------------------------------------------
   subroutine rearrange_with_xtmp(n, dim, idx, gridhash, x, x_tmp)
      ! same as above, except doesn't modify x, and returns x_tmp instead.

      implicit none
      integer, intent(in):: n, dim, idx(n)
      integer, intent(inout):: gridhash(n)
      double precision, intent(in):: x(dim, n)
      integer, allocatable:: tmpgridhash(:)
      double precision, intent(out):: x_tmp(dim, n)
      integer:: i

      allocate (tmpgridhash(n))

      do i = 1, n
         tmpgridhash(i) = gridhash(idx(i))
         x_tmp(:, i) = x(:, idx(i))
      end do

      gridhash(:) = tmpgridhash(:)

   end subroutine rearrange_with_xtmp

   !---------------------------------------------------------------------------
   subroutine find_starts(dim, n, ngridx, gridhash, starts)
      ! finds the starting index for each grid cell. See "Building the grid
      ! using Sorting" section in https://developer.download.nvidia.com/assets/cuda/files/particles.pdf

      implicit none
      integer, intent(in):: dim, n, ngridx(dim)
      integer, intent(in):: gridhash(n)
      integer, intent(out):: starts(:)
      integer:: i

      starts(1:gridhash(1)) = 1

      do i = 2, n
         if (gridhash(i) /= gridhash(i - 1)) starts(gridhash(i - 1) + 1:gridhash(i)) = i
      end do

      starts(gridhash(n) + 1:product(ngridx)) = n + 1

   end subroutine find_starts

   !---------------------------------------------------------------------------
   subroutine update_pair_list(dim, n, maxnpair, cutoff, jstart, jend, i, x, npairs, pairs)
      ! subroutine to loop over sequential grid cells and check if j point is
      ! within a given i point

      implicit none
      integer, intent(in):: jstart, jend, i, n, maxnpair, dim
      double precision, intent(in):: cutoff, x(dim, n)
      integer, intent(inout):: npairs, pairs(2, maxnpair)
      !    integer, intent(inout):: npairs, pair_i(0:), pair_j(0:)
      integer:: j
      double precision:: r2

      do j = jstart, jend
         r2 = sum((x(:, i) - x(:, j))**2)
         if (r2 <= cutoff*cutoff) then
            npairs = npairs + 1
            pairs(1, npairs) = i
            pairs(2, npairs) = j
            !  npairs = npairs + merge(1, 0, r2 <= cutoff*cutoff)
         end if
      end do

   end subroutine update_pair_list

   !----------------------------------------------------------------------------
   function idx_to_hash(dim, idx, ngridx) result(hash)
      ! function to convert an arbitrary-dimensioned index array to a single
      ! hash integer

      implicit none
      integer, intent(in):: dim, idx(dim), ngridx(dim)
      integer:: hash, i, cum_prod

      hash = 1
      cum_prod = 1
      do i = 1, dim
         hash = hash + (idx(i) - 1)*cum_prod
         cum_prod = cum_prod*ngridx(i)
      end do

   end function idx_to_hash

   !----------------------------------------------------------------------------
   function didx_to_dhash(dim, idx, ngridx) result(hash)
      ! function to convert incremental, arbitrary-dimensioned, index array to
      ! a single incremental hash integer

      implicit none
      integer, intent(in):: dim, idx(dim), ngridx(dim)
      integer:: hash, i, cum_prod

      hash = 0
      cum_prod = 1
      do i = 1, dim
         hash = hash + idx(i)*cum_prod
         cum_prod = cum_prod*ngridx(i)
      end do

   end function didx_to_dhash

   !----------------------------------------------------------------------------
   function hash_to_idx(dim, hash, ngridx) result(idx)
      ! function to convert a single hash integer to an arbitrary-dimensioned
      ! index array

      implicit none
      integer, intent(in):: dim, hash, ngridx(dim)
      integer:: idx(dim), i, cum_prod, rem

      rem = hash - 1
      cum_prod = product(ngridx)
      do i = dim, 1, -1
         cum_prod = cum_prod/ngridx(i)
         idx(i) = int(real(rem)/cum_prod) + 1
         rem = mod(hash - 1, cum_prod)
      end do

   end function hash_to_idx

   !----------------------------------------------------------------------------
   subroutine get_adjacent_cells_hash_increment(dim, ngridx, adjHash)
      ! function to determine which adjacent cells to search

      implicit none
      integer, intent(in):: dim, ngridx(dim)
      integer, intent(out):: adjHash(:)
      integer:: nAdj, i, idx(dim), quasi_hash(dim)

      nAdj = (3**(dim - 1) - 1)/2
      quasi_hash(:) = 3

      do i = 0, nAdj - 1
         idx = hash_to_idx(dim, 3*i + 1, quasi_hash) - 2
         adjHash(i + 1) = didx_to_dhash(dim, idx, ngridx)
      end do

   end subroutine get_adjacent_cells_hash_increment

   !----------------------------------------------------------------------------
   function coords_to_hash(n, dim, x, minx, ngridx, cutoff) result(gridhash)
      ! function to convert raw double coordinates to integer hashes

      implicit none
      integer, intent(in):: n, dim, ngridx(dim)
      double precision, intent(in):: x(dim, n), minx(dim), cutoff
      integer:: i, icell(dim), gridhash(n)

      do i = 1, n
         icell(:) = int((x(:, i) - minx(:))/cutoff) + 1
         gridhash(i) = idx_to_hash(dim, icell, ngridx)
      end do

   end function coords_to_hash

   !----------------------------------------------------------------------------
   subroutine query_pairs(dim, npoints, x, cutoff, maxnpair, npairs, pairs)
      ! performs the pair search, and reorders x

      implicit none
      integer, intent(in):: dim, npoints, maxnpair
      double precision, intent(in):: cutoff
      double precision, intent(inout):: x(dim, npoints)
      integer, intent(out):: npairs, pairs(2, maxnpair)
      integer:: i, i_adj, hashi, hashj, ngridx(dim), n_adj
      double precision:: minx(dim), maxx(dim)
      integer, allocatable:: gridhash(:), idx(:), starts(:), adj_cell_hash_increment(:)

      ! determine bounding box from raw coordinates
      minx(:) = minval(x, dim=2)
      maxx(:) = maxval(x, dim=2)

      ! extend bounding box by two cells in each dimensions
      minx(:) = minx(:) - 2.d0*cutoff
      maxx(:) = maxx(:) + 2.d0*cutoff

      ! calculate number of grid cells and adjust maximum extent
      ngridx(:) = int((maxx(:) - minx(:))/cutoff) + 1
      maxx(:) = maxx(:) + ngridx(:)*cutoff

      ! convert coordinates to integer hashes
      allocate (gridhash(npoints))
      gridhash = coords_to_hash(npoints, dim, x, minx, ngridx, cutoff)

      ! perform argsort with idx as keys
      allocate (idx(npoints))
      call sort_hashes(npoints, gridhash, idx)

      ! rearrange gridhash and x using keys
      call rearrange(npoints, dim, idx, gridhash, x)

      ! find where each cell starts in the point array
      allocate (starts(product(ngridx)))
      call find_starts(dim, npoints, ngridx, gridhash, starts)

      ! calculate which adjacent cells to search, based on the dimension of the
      ! problem
      n_adj = (3**(dim - 1) - 1)/2
      allocate (adj_cell_hash_increment(n_adj))
      call get_adjacent_cells_hash_increment(dim, ngridx, adj_cell_hash_increment)

      ! perform sweep to find pairs
      npairs = 0
      do hashi = gridhash(1), gridhash(npoints)
         do i = starts(hashi), starts(hashi + 1) - 1
            call update_pair_list(dim, npoints, maxnpair, cutoff, i + 1, starts(hashi + 2) - 1, i, x, npairs, pairs)
            do i_adj = 1, n_adj
               hashj = hashi + adj_cell_hash_increment(i_adj)
               call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x, &
                                     npairs, pairs)
            end do
         end do
      end do

      pairs(:, 1:npairs) = pairs(:, 1:npairs) - 1  ! -1 for Python 0-indexing

   end subroutine query_pairs

   subroutine query_pairs_noreorder(dim, npoints, x, cutoff, maxnpair, npairs, pairs)
      ! performs the pair search, without reordering x.
      ! steps are the same, except the rearrange_withtmpx subroutine is used
      ! instead of rearrange, and subsequently, x_tmp is used instead of x.

      implicit none
      integer, intent(in):: dim, npoints, maxnpair
      double precision, intent(in):: cutoff, x(dim, npoints)
      integer, intent(out):: npairs, pairs(2, maxnpair)
      integer:: i, i_adj, hashi, hashj, ngridx(dim), n_adj
      double precision:: minx(dim), maxx(dim)
      integer, allocatable:: gridhash(:), idx(:), starts(:), adj_cell_hash_increment(:)
      double precision, allocatable:: x_tmp(:, :)

      minx(:) = minval(x, dim=2)
      maxx(:) = maxval(x, dim=2)

      minx(:) = minx(:) - 2.d0*cutoff
      maxx(:) = maxx(:) + 2.d0*cutoff

      ngridx(:) = int((maxx(:) - minx(:))/cutoff) + 1
      maxx(:) = maxx(:) + ngridx(:)*cutoff

      allocate (gridhash(npoints))
      gridhash = coords_to_hash(npoints, dim, x, minx, ngridx, cutoff)

      allocate (idx(npoints), x_tmp(dim, npoints))
      call sort_hashes(npoints, gridhash, idx)

      call rearrange_with_xtmp(npoints, dim, idx, gridhash, x, x_tmp)

      allocate (starts(product(ngridx)))
      call find_starts(dim, npoints, ngridx, gridhash, starts)

      n_adj = (3**(dim - 1) - 1)/2
      allocate (adj_cell_hash_increment(n_adj))
      call get_adjacent_cells_hash_increment(dim, ngridx, adj_cell_hash_increment)

      npairs = 0
      do hashi = gridhash(1), gridhash(npoints)
         do i = starts(hashi), starts(hashi + 1) - 1
            call update_pair_list(dim, npoints, maxnpair, cutoff, i + 1, starts(hashi + 2) - 1, i, x_tmp, npairs, &
                                  pairs)
            do i_adj = 1, n_adj
               hashj = hashi + adj_cell_hash_increment(i_adj)
               call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x_tmp, &
                                     npairs, pairs)
            end do
         end do
      end do

      do i = 1, npairs
         pairs(1, i) = idx(pairs(1, i)) - 1 ! -1 for Python 0-indexing
         pairs(2, i) = idx(pairs(2, i)) - 1
      end do

   end subroutine query_pairs_noreorder

end module closefriends
