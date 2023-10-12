#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <set>
#include <tuple>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

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
py::array_t<int, py::array::c_style> sortHashes(int n, std::vector<int> hashes) {

    py::array_t<int, py::array::c_style> idx(n);
    int* idx_ptr = idx.mutable_data();

    // filling idx with numbers 0:n-1
    for (int i = 0; i < n; ++i) idx_ptr[i] = i;

    // sorting idx based on hashes
    std::sort(idx_ptr, idx_ptr+n, 
        [&](const int& a, const int& b) {
            return (hashes[a] < hashes[b]);
        }
    );

    return idx;
}

// same as above, except doesn't modify x, and returns tmp_x instead ---------------------------------------------------
void rearrangeWithTmpx(const int ndims, 
                       const int npoints, 
                       const py::detail::unchecked_reference<int, 1> idx, 
                       std::vector<int> &gridhash, 
                       const py::detail::unchecked_reference<double, 2> &x, 
                       py::detail::unchecked_mutable_reference<double, 2> &tmp_x) {

    std::vector<int> tmp_gridhash(npoints);
    
    for (int i = 0; i < npoints; ++i) {
        tmp_gridhash[i] = gridhash[idx(i)];
        for (int d = 0; d < ndims; ++d)
            tmp_x(i, d) = x(idx(i), d);
    }

    for (int i = 0; i < npoints; ++i)
        gridhash[i] = tmp_gridhash[i];

}

// convert raw double coordinates to integer hashes --------------------------------------------------------------------
std::vector<int> coordsToHash(const int ndims, 
                              const int npoints, 
                              const py::detail::unchecked_reference<double, 2> &x, 
                              const std::vector<double> &minx, 
                              const std::vector<int> &ngridx,
                              const double cutoff) {
    std::vector<int> gridhash(npoints), icell(ndims);
    for (int i = 0; i < npoints; ++i) {
        for (int d = 0; d < ndims; ++d)
            icell[d] = static_cast<int>((x(i, d)-minx[d])/cutoff);
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
    for (int i = 0; i <= gridhash[0]; ++i) starts[i] = 0;
    for (int i = 1; i < npoints; ++i)
        if (gridhash[i] != gridhash[i-1])
            for (int j = gridhash[i-1]+1; j <= gridhash[i]; ++j) 
                starts[j] = i;
    for (int i = gridhash[npoints-1]+1; i < ngridx_prod; ++i) starts[i] = npoints;

    return starts;
}

// determine which adjacent cells to search ----------------------------------------------------------------------------
std::vector<int> getAdjacentCellsHashIncrement(const int ndims, const std::vector<int> &ngridx) {
    const int nAdj = (std::pow(3, ndims-1)-1)/2;
    std::vector<int> adjHashInc(nAdj);
    const std::vector<int> quasi_ngridx(ndims, 3);
    for (int i = 0; i < nAdj; ++i) {
        std::vector<int> idx = hashToIdx(ndims, 3*i, quasi_ngridx);
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
                    const py::detail::unchecked_mutable_reference<double, 2>  &x,
                    int &npairs,
                    py::detail::unchecked_mutable_reference<int, 1> pair_i_data,
                    py::detail::unchecked_mutable_reference<int, 1> pair_j_data) {
    for (int j = jstart; j < jend; ++j) {
        double r2 = 0.;
        for (int d = 0; d < ndims; ++d) 
            r2 += (x(i, d)-x(j, d))*(x(i, d)-x(j, d));
        if (r2 <= cutoff*cutoff) {
            pair_i_data(npairs) = i;
            pair_j_data(npairs) = j;
            npairs = std::min(maxnpair, npairs + 1);
        }
    }
}

// overwrite input array
void overwrite_input_x(py::array_t<double, py::array::c_style> target, py::array_t<double, py::array::c_style> source) {
    py::buffer_info target_info = target.request();
    py::buffer_info source_info = source.request();
    double* target_ptr = static_cast<double*>(target_info.ptr);
    double* source_ptr = static_cast<double*>(source_info.ptr);
    std::copy(source_ptr, source_ptr+target_info.size, target_ptr);
}

// performs the pair search, without reordering x ----------------------------------------------------------------------
py::object
query_pairs(py::array_t<double, py::array::c_style> &input_array, 
            double cutoff, 
            int maxnpair, 
            const std::string &output_type = "seperate-ndarrays",
            const bool retain_order = false) {

    // create wrapper pointer for 2D indexing
    py::detail::unchecked_reference<double, 2> x = input_array.unchecked<2>();
    const int npoints = input_array.request().shape[0]; // get npoints from input array
    const int ndims = input_array.request().shape[1]; // get dims from input array

    // determine bounding box from raw coordinates
    std::vector<double> mingridx(ndims, std::numeric_limits<double>::max()), 
        maxgridx(ndims, -std::numeric_limits<double>::max());

    for (int i = 0; i < npoints; ++i)
        for (int d = 0; d < ndims; ++d) {
            mingridx[d] = std::min(mingridx[d], x(i, d));
            maxgridx[d] = std::max(maxgridx[d], x(i, d));
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
    py::array_t<int, py::array::c_style> idx = sortHashes(npoints, gridhash);
    py::detail::unchecked_reference<int, 1> idx_data = idx.unchecked<1>();

    // copy input data into tmp_input_array and create wrapper pointer
    py::array_t<double, py::array::c_style> tmp_input_array = input_array.attr("copy")();
    py::detail::unchecked_mutable_reference<double, 2> tmp_x = tmp_input_array.mutable_unchecked<2>();

    // rearrange gridhash and x using keys
    rearrangeWithTmpx(ndims, npoints, idx_data, gridhash, x, tmp_x);

    // find where each cell starts in the point array
    std::vector<int> starts = findStarts(ndims, npoints, ngridx, gridhash);

    // calculate which adjacent cells to search, based on the dimension of the problem
    std::vector<int> adjCellHashIncrement = getAdjacentCellsHashIncrement(ndims, ngridx);

    // perform sweep to find pairs
    // py::array_t<int, py::array::c_style> pairs({maxnpair, 2});
    // py::detail::unchecked_mutable_reference<int, 2> pairs_data = pairs.mutable_unchecked<2>();
    py::array_t<int, py::array::c_style> pair_i(maxnpair), pair_j(maxnpair);
    py::detail::unchecked_mutable_reference<int, 1> pair_i_data = pair_i.mutable_unchecked<1>(), pair_j_data = pair_j.mutable_unchecked<1>();
    int npairs = 0;
    for (int hashi = gridhash[0]; hashi <= gridhash[npoints-1]; ++hashi)
        for (int i = starts[hashi]; i < starts[hashi+1]; ++i) {
            updatePairList(ndims, npoints, maxnpair, cutoff, i+1, starts[hashi+2], i, tmp_x, npairs, pair_i_data, pair_j_data);
            for (uint iAdj = 0; iAdj < adjCellHashIncrement.size(); ++iAdj) {
                const int hashj = hashi + adjCellHashIncrement[iAdj];
                updatePairList(ndims, npoints, maxnpair, cutoff, starts[hashj], starts[hashj+3], i, tmp_x, npairs, pair_i_data, pair_j_data);
            }
        }

    if (npairs==maxnpair) throw py::buffer_error("Max number of pairs, as indicated by maxnpair, isn't large enough!");

    // modify pair indices to refer to original positions
    // also ensures pairs[k, 0] < pairs[k, 1]
    if (retain_order) {
        for (int k = 0; k < npairs; ++k) {
            int i = idx_data(*pair_i.data(k));
            int j = idx_data(*pair_j.data(k));
            pair_i_data(k) = std::min(i, j);
            pair_j_data(k) = std::max(i, j);
        }
    } else {
        overwrite_input_x(input_array, tmp_input_array);
    }

    // truncate array to exclude empty entries of pairs
    py::slice row_slice(0, npairs, 1);
    py::array_t<int> pair_i_trunc(pair_i[row_slice]), pair_j_trunc(pair_j[row_slice]);

    if (output_type == "seperate-ndarrays") {

        if (retain_order) {
            return py::make_tuple(pair_i_trunc, pair_j_trunc);
        } else {
            return py::make_tuple(pair_i_trunc, pair_j_trunc, idx);
        }

    } else if (output_type == "ndarray") {

        py::array_t<int, py::array::c_style> pairs({npairs, 2});
        py::detail::unchecked_mutable_reference<int, 2> pairs_data = pairs.mutable_unchecked<2>();

        for (int k = 0; k < npairs; ++k) {
            pairs_data(k, 0) = (*pair_i_trunc.data(k));
            pairs_data(k, 1) = (*pair_j_trunc.data(k));
        }

        if (retain_order) {
            return pairs;
        } else {
            return py::make_tuple(pairs, idx);
        }

    } else if (output_type == "set") {

        // perform sweep to find pairs
        std::set<std::tuple<int, int>> pairs_set;
        // py::detail::unchecked_reference<int, 2> rarray_data = rarray.unchecked<2>();
        for (int i = 0; i < npairs; ++i) {
            std::tuple<int, int> pair = std::make_tuple((*pair_i_trunc.data(i)), *(pair_j_trunc.data(i)));
            pairs_set.insert(pair);
        }

        if (retain_order) {
            return py::cast(pairs_set);
        } else {
            return py::make_tuple(py::cast(pairs_set), idx);
        }
    } else {
        throw std::invalid_argument("name '" + output_type + "' is not defined");
    }
}

// register module
PYBIND11_MODULE(closefriends, m) {
    m.def("query_pairs", &query_pairs, "A function which sums up all the elements of a 2D array", 
        py::arg("data"), 
        py::arg("r"), 
        py::arg("maxnpair"),
        py::arg("output_type") = "seperate-ndarrays",
        py::arg("retain_order") = false);
}