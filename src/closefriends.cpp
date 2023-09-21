#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>

namespace py = pybind11;

// return a pointer for 2D indexing (might be unnecessary) -------------------------------------------------------------
template <typename T>
std::vector<T*> indexArray(const py::array_t<T, py::array::c_style> &input_array) {
    py::buffer_info buf_info = input_array.request();
    // if (buf_info.ndim != 2)
        // throw std::runtime_error("Expected a 2D array");

    T* ptr = static_cast<T *>(buf_info.ptr);
    const int npoints = buf_info.shape[0];
    const int ndims = buf_info.shape[1];

    std::vector<T*> arr(npoints);

    for (int i = 0; i < npoints; ++i, ptr += ndims)
        arr[i] = ptr;

    return arr;
}

// convert an arbitrary-dimensioned index array to a single hash integer -----------------------------------------------
int idxToHash(const int ndims, const std::vector<int> &idx, const std::vector<int> &ngridx) {
    int hash = 0, cum_prod = 1;
    for (int d = 0; d < ndims; ++d) {
        hash += idx[d]*cum_prod;
        cum_prod *= ngridx[d];
    }
    return hash;
}

// convert a single hash integer to an arbitrary-dimensioned index array -----------------------------------------------
std::vector<int> hashToIdx(const int ndims, const int hash, const std::vector<int> &ngridx) {
    std::vector<int> idx(ndims);
    int rem = hash, cum_prod = 1;
    for (int d = 0; d < ndims; ++d) 
        cum_prod *= ngridx[d];
    for (int i = ndims-1; i >= 0; --i) {
        cum_prod /= ngridx[i];
        idx[i] = static_cast<int>(static_cast<double>(rem)/cum_prod);
        rem = hash % cum_prod;
    }
    return idx;
}

// performs argsort ----------------------------------------------------------------------------------------------------
std::vector<int> sortHashes(int n, std::vector<int> hashes) {

    std::vector<int> idx(n);

    // filling idx with numbers 0:n-1
    for (int i = 0; i < n; ++i) idx[i] = i;

    // sorting idx based on hashes
    std::sort(idx.begin(), idx.end(), 
        [&](const int& a, const int& b) {
            return (hashes[a] < hashes[b]);
        }
    );

    return idx;
}

// same as above, except doesn't modify x, and returns tmp_x instead ---------------------------------------------------
void rearrangeWithTmpx(const int ndims, 
                       const int npoints, 
                       const std::vector<int> &idx, 
                       std::vector<int> &gridhash, 
                       const std::vector<double*> &x, 
                       std::vector<double*> &tmp_x) {

    std::vector<int> tmp_gridhash(npoints);
    
    for (int i = 0; i < npoints; ++i) {
        tmp_gridhash[i] = gridhash[idx[i]];
        for (int d = 0; d < ndims; ++d)
            tmp_x[i][d] = x[idx[i]][d];
    }

    for (int i = 0; i < npoints; ++i)
        gridhash[i] = tmp_gridhash[i];

}

// convert raw double coordinates to integer hashes --------------------------------------------------------------------
std::vector<int> coordsToHash(const int ndims, 
                              const int npoints, 
                              const std::vector<double*> &x, 
                              const std::vector<double> &minx, 
                              const std::vector<int> &ngridx,
                              const double cutoff) {
    std::vector<int> gridhash(npoints), icell(ndims);
    for (int i = 0; i < npoints; ++i) {
        for (int d = 0; d < ndims; ++d)
            icell[d] = static_cast<int>((x[i][d]-minx[d])/cutoff);
        gridhash[i] = idxToHash(ndims, icell, ngridx);
    }
    return gridhash;
}

// finds the starting index for each grid cell. See "Building the grid using Sorting" section in -----------------------
// https://developer.download.nvidia.com/assets/cuda/files/particles.pdf
std::vector<int> findStarts(const int ndims, 
                            const int npoints, 
                            const std::vector<int> &ngridx, 
                            const std::vector<int> &gridhash) {
    int ngridx_prod = 1;
    for (int d = 0; d < ndims; ++d) 
        ngridx_prod *= ngridx[d];
    std::vector<int> starts(ngridx_prod);
    for (int i = 0; i <= gridhash[0]; ++i) 
        starts[i] = 0;
    for (int i = 1; i < npoints; ++i)
        if (gridhash[i] != gridhash[i-1])
            for (int j = gridhash[i-1]+1; j <= gridhash[i]; ++j) 
                starts[j] = i;
    for (int i = gridhash[npoints-1]; i < ngridx_prod; ++i) starts[i] = npoints;

    return starts;
}

// determine which adjacent cells to search ----------------------------------------------------------------------------
std::vector<int> getAdjacentCellsHashIncrement(const int ndims, const std::vector<int> &ngridx) {
    const int nAdj = (std::pow(3, ndims-1)-1)/2;
    std::vector<int> adjHashInc(nAdj);
    const std::vector<int> quasi_hash(ndims, 3);
    for (int i = 0; i < nAdj; ++i) {
        std::vector<int> idx = hashToIdx(ndims, 3*i, quasi_hash);
        for (int d = 0; d < ndims; ++d)
            idx[d] -= 1;
        adjHashInc[i] = idxToHash(ndims, idx, ngridx);
    }
    return adjHashInc;
}

// loop over sequential grid cells and check if j point is within a given i point --------------------------------------
void updatePairList(const int ndims, 
                    const int npoints,
                    const int maxnpair, 
                    const double cutoff,
                    const int jstart,
                    const int jend, 
                    const int i,
                    const std::vector<double*> &x,
                    int &npairs,
                    py::array_t<int, py::array::c_style> pairs) {
    for (int j = jstart; j < jend; ++j) {
        double r2 = 0.;
        for (int d = 0; d < ndims; ++d) 
            r2 += (x[i][d]-x[j][d])*(x[i][d]-x[j][d]);
        if (r2 <= cutoff*cutoff) {
            *pairs.mutable_data(npairs, 0) = i;
            *pairs.mutable_data(npairs, 1) = j;
            npairs++;
        }
    }
}

// performs the pair search, without reordering x ----------------------------------------------------------------------
py::array_t<int, py::array::c_style> 
query_pairs(py::array_t<double, py::array::c_style> &input_array, double cutoff, int maxnpair) {

    // create wrapper pointer for 2D indexing
    const std::vector<double*> x = indexArray(input_array);
    const int npoints = input_array.request().shape[0]; // get npoints from input array
    const int ndims = input_array.request().shape[1]; // get dims from input array

    // determine bounding box from raw coordinates
    std::vector<double> mingridx(ndims, std::numeric_limits<double>::max()), 
        maxgridx(ndims, -std::numeric_limits<double>::max());

    for (int i = 0; i < npoints; ++i)
        for (int d = 0; d < ndims; ++d) {
            mingridx[d] = std::min(mingridx[d], x[i][d]);
            maxgridx[d] = std::max(maxgridx[d], x[i][d]);
        }

    // extend bounding box by two cells in each dimensions, calculate number of grid cells and adjust maximum extent
    std::vector<int> ngridx(ndims);
    for (int d = 0; d < ndims; ++d) {
        mingridx[d] -= 2.*cutoff;
        maxgridx[d] += 2.*cutoff;
        ngridx[d] = static_cast<int>((maxgridx[d]-mingridx[d])/cutoff) + 1;
        maxgridx[d] = mingridx[d] + ngridx[d]*cutoff;
    }

    // convert coordinates to integer hashes
    std::vector<int> gridhash = coordsToHash(ndims, npoints, x, mingridx, ngridx, cutoff);

    // perform argsort with idx as keys
    std::vector<int> idx = sortHashes(npoints, gridhash);

    // copy input data into tmp_input_array and create wrapper pointer
    py::array_t<double, py::array::c_style> tmp_input_array = input_array.attr("copy")();
    std::vector<double*> tmp_x = indexArray(tmp_input_array);

    // rearrange gridhash and x using keys
    rearrangeWithTmpx(ndims, npoints, idx, gridhash, x, tmp_x);

    // find where each cell starts in the point array
    std::vector<int> starts = findStarts(ndims, npoints, ngridx, gridhash);

    // calculate which adjacent cells to search, based on the dimension of the problem
    std::vector<int> adjCellHashIncrement = getAdjacentCellsHashIncrement(ndims, ngridx);

    // perform sweep to find pairs
    py::array_t<int, py::array::c_style> pairs({maxnpair, 2});
    int npairs = 0;
    for (int hashi = gridhash[0]; hashi <= gridhash[npoints-1]; ++hashi)
        for (int i = starts[hashi]; i < starts[hashi+1]; ++i) {
            updatePairList(ndims, npoints, maxnpair, cutoff, i+1, starts[hashi+2], i, tmp_x, npairs, pairs);
            for (uint iAdj = 0; iAdj < adjCellHashIncrement.size(); ++iAdj) {
                const int hashj = hashi + adjCellHashIncrement[iAdj];
                updatePairList(ndims, npoints, maxnpair, cutoff, starts[hashj], starts[hashj+3], i, tmp_x, npairs, pairs);
            }
        }

    // modify pair indicies to refer to original positions
    for (int i = 0; i < npairs; ++i) {
        *pairs.mutable_data(i, 0) = idx[*pairs.data(i, 0)];
        *pairs.mutable_data(i, 1) = idx[*pairs.data(i, 1)];
    }

    // truncate array to exclude empty entries of pairs
    py::slice slice_row(0, npairs, 1);
    py::array rarray = pairs[slice_row];
    return rarray;
}

// register module
PYBIND11_MODULE(closefriends, m) {
    m.def("query_pairs", &query_pairs, "A function which sums up all the elements of a 2D array", 
        py::arg("data"), 
        py::arg("r"), 
        py::arg("maxnpair"));
}