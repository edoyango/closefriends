#include <algorithm>
#include <numeric>

extern "C" {
    
    void sort_hashes(int n, int* hashes, int* idx) {

        // filling idx with numbers 0:n-1
        for (int i = 0; i < n; ++i) {
            idx[i] = i;
        }

        // sorting idx based on hashes
        std::sort(idx, idx + n, 
            [&](const int& a, const int& b) {
                return (hashes[a] < hashes[b]);
            }
        );

        // adapting idx for 1-based indexing
        for (int i = 0; i < n; ++i) {
            idx[i]++;
        }
    }
}