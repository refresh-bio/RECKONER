/*
* RECKONER - Read Error Corrector Based on KMC
*
* This software is distributed under GNU GPL 3 license.
*
* Authors: Yun Heo, Maciej Dlugosz
* Version: 2.0
*
*/

#ifndef QUERYKMCDB_H
#define QUERYKMCDB_H


#include <cassert>
#include <kmc_api/kmc_file.h>
#include <kmc_api/kmer_api.h>



class QueryKMCDb {
private:
    CKMCFile& kmc_file;
    CKMCFile* kmc_long_file;
    CKmerAPI kmer_api;
    CKmerAPI kmer_long_api;

    uint32 kmer_length;
    uint32 long_kmer_length;

    template<typename KMER_TYPE>
    bool get_canonical_kmer(const KMER_TYPE& kmer, std::string& canonical, uint32 kmer_length);

public:
    QueryKMCDb(CKMCFile& _kmc_file, CKMCFile* _kmc_long_file, uint32 _kmer_length, uint32 _long_kmer_length) :
        kmc_file(_kmc_file),
        kmc_long_file(_kmc_long_file),
        kmer_api(_kmer_length),
        kmer_long_api(_long_kmer_length),
        kmer_length(_kmer_length),
        long_kmer_length(_long_kmer_length)
    {}

    bool use_long_kmer() { return kmc_long_file != nullptr; }

    template<typename KMER_TYPE>
    bool query_text(const KMER_TYPE& kmer, float& kmer_quality);

    template<typename KMER_TYPE>
    bool query_text(const KMER_TYPE& kmer, float& kmer_quality, std::size_t& num_forward_kmers, std::size_t& num_reverse_kmers);

    template<typename KMER_TYPE>
    bool query_text_long_kmer(const KMER_TYPE& kmer, float& kmer_quality);
};



//----------------------------------------------------------------------
// Determines the canonical k-mer sequence (if it is not the same, as given k-mer).
//----------------------------------------------------------------------

template<typename KMER_TYPE>
bool QueryKMCDb::get_canonical_kmer(const KMER_TYPE& kmer, std::string& canonical, uint32 kmer_length) {
    std::string reverse_complement(kmer_length, '0');
    bool surely_reverse = false;
    for (int i = 0, i_comp = static_cast<int>(kmer_length) - 1; i < static_cast<int>(kmer_length); ++i, --i_comp) {
        char symbol = 'A';
        if (kmer[i_comp] == 'A') {
            symbol = 'T';
        }
        else if (kmer[i_comp] == 'C') {
            symbol = 'G';
        }
        else if (kmer[i_comp] == 'G') {
            symbol = 'C';
        }

        if (kmer[i] == symbol || surely_reverse) {
            reverse_complement[i] = symbol;
        }
        else if (kmer[i] < symbol) {
            return false;
        }
        else if (kmer[i] > symbol) {
            reverse_complement[i] = symbol;
            surely_reverse = true;
        }
    }
    if (surely_reverse) {
        canonical.swap(reverse_complement);
        return true;
    }
    return false;
}



//----------------------------------------------------------------------
// Checks in the KMC database if k-mer exists, if yes it returns kmer counter.
//----------------------------------------------------------------------

template<typename KMER_TYPE>
inline bool QueryKMCDb::query_text(const KMER_TYPE& kmer, float& kmer_quality) {
    std::string canonical;
    if (get_canonical_kmer(kmer, canonical, kmer_length)) {
        const bool res = kmer_api.from_string_low_level(canonical);
        assert(res);
    }
    else {
        const bool res = kmer_api.from_string_low_level(kmer);
        assert(res);
    }

    return kmc_file.CheckKmer(kmer_api, kmer_quality);
}



//----------------------------------------------------------------------
// Checks in the KMC database if k-mer exists, if yes it returns kmer counter. Determines a number of reverses.
//----------------------------------------------------------------------

template<typename KMER_TYPE>
inline bool QueryKMCDb::query_text(const KMER_TYPE& kmer, float& kmer_quality, std::size_t& num_forward_kmers, std::size_t& num_reverse_kmers) {
    assert(kmer_api.from_string_low_level(kmer));

    if (kmc_file.CheckKmer(kmer_api, kmer_quality)) {
        ++num_forward_kmers;
        return true;
    }

    kmer_api.reverse();

    if (kmc_file.CheckKmer(kmer_api, kmer_quality)) {
        ++num_reverse_kmers;
        return true;
    }

    return false;
}



//----------------------------------------------------------------------
// Checks in the KMC database if a long k-mer exists, if yes it returns kmer counter.
//----------------------------------------------------------------------

template<typename KMER_TYPE>
inline bool QueryKMCDb::query_text_long_kmer(const KMER_TYPE& kmer, float& kmer_quality) {
    std::string canonical;
    if (get_canonical_kmer(kmer, canonical, long_kmer_length)) {
        const bool res = kmer_long_api.from_string_low_level(canonical);
        assert(res);
    }
    else {
        const bool res = kmer_long_api.from_string_low_level(kmer);
        assert(res);
    }

    return kmc_long_file->CheckKmer(kmer_long_api, kmer_quality);
}

#endif
