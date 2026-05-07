#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include "rshash.hpp"
#include "contig_table.hpp"
#include "projected_hits.hpp"

class ReferenceIndex {
public:
    ReferenceIndex() = default;

    void load(const std::string& basename) {
        std::cout << "loading RSHash dictionary...\n";
        m_dict.load(basename + ".rshash");

        std::cout << "loading contig table...\n";
        {
            std::ifstream in(basename + ".ctab", std::ios::binary);
            cereal::BinaryInputArchive ar(in);
            ar(m_ctab);
        }

        std::cout << "loading reference info...\n";
        {
            std::ifstream in(basename + ".refinfo", std::ios::binary);
            cereal::BinaryInputArchive ar(in);
            ar(m_ref_names, m_ref_lens);
        }

        std::cout << "loaded reference index: " << m_ref_names.size()
                  << " references, " << m_dict.num_real_unitigs() << " unitigs\n";
    }

    projected_hits query_kmer(uint64_t kmer_fw, uint64_t kmer_rc) const {
        constexpr uint32_t invalid_u32 = std::numeric_limits<uint32_t>::max();

        auto loc = m_dict.locate_kmer(kmer_fw, kmer_rc);
        if (!loc.has_value()) {
            return {invalid_u32, invalid_u32, false, invalid_u32,
                    static_cast<uint32_t>(m_dict.getk()), {}};
        }

        auto& r = loc.value();
        auto span = m_ctab.entries(r.unitig_id);

        return projected_hits{
            static_cast<uint32_t>(r.unitig_id),
            r.contig_pos,
            r.is_forward,
            r.contig_len,
            static_cast<uint32_t>(m_dict.getk()),
            span
        };
    }

    void streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4>& query,
                         std::vector<std::optional<LocateResult>>& results,
                         uint64_t& extensions) const {
        m_dict.streaming_query_ref(query, results, extensions);
    }

    void streaming_query_cb(const seqan3::bitpacked_sequence<seqan3::dna4>& query,
                            uint64_t& extensions, uint64_t* offsets_buf, uint64_t* kmer_buf,
                            void* ctx, streaming_hit_fn_t on_hit) const {
        m_dict.streaming_query_ref_cb(query, extensions, offsets_buf, kmer_buf, ctx, on_hit);
    }

    std::pair<uint64_t, uint64_t> query_buffer_sizes() const { return m_dict.query_buffer_sizes(); }

    uint64_t k() const { return m_dict.getk(); }
    uint64_t num_refs() const { return m_ref_names.size(); }
    const std::string& ref_name(size_t i) const { return m_ref_names[i]; }
    uint64_t ref_len(size_t i) const { return m_ref_lens[i]; }
    ContigTable& contig_table() const { return m_ctab; }
    const RSHash& dict() const { return m_dict; }

private:
    mutable RSHash m_dict;
    mutable ContigTable m_ctab;
    std::vector<std::string> m_ref_names;
    std::vector<uint64_t> m_ref_lens;
};
