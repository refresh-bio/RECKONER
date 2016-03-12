#include <kmc_api/kmer_api.h>
#include <kmc_api/kmc_file.h>
#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------
// Main function.
//----------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "ERROR: Not enough arguments." << endl;
        exit(EXIT_FAILURE);
    }

    CKMCFile kmc_file;
    if (!kmc_file.OpenForListing(argv[1])) {
        cerr << "ERROR: cannot open KMC files." << endl;
        exit(EXIT_FAILURE);
    }

    uint32 k, mode, counter_size, prefix_length, signature_length, min_count;
    uint64 max_count, total;
    kmc_file.Info(k, mode, counter_size, prefix_length, signature_length, min_count, max_count, total);

    CKmerAPI kmer(k);
    float count;
    int iCount;

    const int histoSize = 256;
	const int maxCutoff = 5;
    vector<uint64> histogram;
    histogram.resize(histoSize, 0);

    while (kmc_file.ReadNextKmer(kmer, count)) {
        iCount = static_cast<int> (round(count));
        if (iCount >= histoSize) {
            iCount = histoSize - 1;
        }

        ++histogram[iCount];
    }

    uint64 firstMin = histogram[2];

	int cutoff = 2;
    for (int i = 3; i < histoSize; ++i) {
        if (firstMin >= histogram[i]) {
            firstMin = histogram[i];
        }
        else {
            cutoff = i - 1;
            break;
        }
    }
	cout << (cutoff > maxCutoff ? maxCutoff : cutoff) << endl;
}
