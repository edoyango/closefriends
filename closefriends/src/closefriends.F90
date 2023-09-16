module closefriends

   implicit none

   interface sort_hashes
      subroutine sort_hashes(n, hashes, idx) bind(C)
         use iso_c_binding
         integer(c_int), value, intent(in):: n
         integer(c_int), intent(in):: hashes(*)
         integer(c_int), intent(out):: idx(*)
      end subroutine sort_hashes
   end interface sort_hashes

contains

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

   function idxToHash(dim, idx, ngridx) result(hash)

      implicit none
      integer, intent(in):: dim, idx(dim), ngridx(dim)
      integer:: hash, i, cum_prod

      hash = 1
      cum_prod = 1
      do i = 1, dim
         hash = hash + (idx(i) - 1)*cum_prod
         cum_prod = cum_prod*ngridx(i)
      end do

   end function idxToHash

   function didxTodHash(dim, idx, ngridx) result(hash)

      implicit none
      integer, intent(in):: dim, idx(dim), ngridx(dim)
      integer:: hash, i, cum_prod

      hash = 0
      cum_prod = 1
      do i = 1, dim
         hash = hash + idx(i)*cum_prod
         cum_prod = cum_prod*ngridx(i)
      end do

   end function didxTodHash

   function hashToIdx(dim, hash, ngridx) result(idx)

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

   end function hashToIdx

   subroutine getAdjacentCellsHashIncrement(dim, ngridx, adjHash)

      implicit none
      integer, intent(in):: dim, ngridx(dim)
      integer, intent(out):: adjHash(:)
      integer:: nAdj, i, idx(dim), quasi_hash(dim)

      nAdj = (3**(dim - 1) - 1)/2
      quasi_hash(:) = 3

      do i = 0, nAdj - 1
         idx = hashToIdx(dim, 3*i + 1, quasi_hash) - 2
         adjHash(i + 1) = didxTodHash(dim, idx, ngridx)
      end do

   end subroutine getAdjacentCellsHashIncrement

   function coordsToHash(n, dim, x, minx, ngridx, cutoff) result(gridhash)

      implicit none
      integer, intent(in):: n, dim, ngridx(dim)
      double precision, intent(in):: x(dim, n), minx(dim), cutoff
      integer:: i, icell(dim), gridhash(n)

      do i = 1, n
         icell(:) = int((x(:, i) - minx(:))/cutoff) + 1
         gridhash(i) = idxToHash(dim, icell, ngridx)
      end do

   end function coordsToHash

   subroutine cellList(dim, npoints, x, cutoff, maxnpair, npairs, pair_i, pair_j)

      implicit none
      integer, intent(in):: dim, npoints, maxnpair
      double precision, intent(in):: cutoff
      double precision, intent(inout):: x(dim, npoints)
      integer, intent(out):: npairs, pair_i(maxnpair), pair_j(maxnpair)
      integer:: i, i_adj, hashi, hashj, ngridx(dim), n_adj
      double precision:: minx(dim), maxx(dim)
      integer, allocatable:: gridhash(:), idx(:), starts(:), adj_cell_hash_increment(:)

      minx(:) = minval(x, dim=2)
      maxx(:) = maxval(x, dim=2)

      minx(:) = minx(:) - 2.d0*cutoff
      maxx(:) = maxx(:) + 2.d0*cutoff

      ngridx(:) = int((maxx(:) - minx(:))/cutoff) + 1
      maxx(:) = maxx(:) + ngridx(:)*cutoff

      allocate (gridhash(npoints))
      gridhash = coordsToHash(npoints, dim, x, minx, ngridx, cutoff)

      allocate (idx(npoints))
      call sort_hashes(npoints, gridhash, idx)

      call rearrange(npoints, dim, idx, gridhash, x, undo=.false.)

      allocate (starts(product(ngridx)))
      call find_starts(dim, npoints, ngridx, gridhash, starts)

      n_adj = (3**(dim - 1) - 1)/2
      allocate (adj_cell_hash_increment(n_adj))
      call getAdjacentCellsHashIncrement(dim, ngridx, adj_cell_hash_increment)

      npairs = 0
      do hashi = gridhash(1), gridhash(npoints)
         do i = starts(hashi), starts(hashi + 1) - 1
            call update_pair_list(dim, npoints, maxnpair, cutoff, i + 1, starts(hashi + 2) - 1, i, x, npairs, pair_i, &
                                  pair_j)
            do i_adj = 1, n_adj
               hashj = hashi + adj_cell_hash_increment(i_adj)
               call update_pair_list(dim, npoints, maxnpair, cutoff, starts(hashj), starts(hashj + 3) - 1, i, x, &
                                     npairs, pair_i, pair_j)
            end do
         end do
      end do

   end subroutine cellList

end module closefriends
