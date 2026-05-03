#pragma once

#include <vector>
#include <cmath>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "compact_vector.hpp"
#include "EliasFano.hpp"

struct contig_span {
    ::bits::compact_vector::iterator start_;
    ::bits::compact_vector::iterator stop_;
    size_t len_;

    bool empty() const { return len_ == 0; }
    size_t size() const { return len_; }
    auto begin() const { return start_; }
    auto end() const { return stop_; }
};

class ContigTable {
public:
    uint64_t m_ref_len_bits{0};
    uint64_t m_ref_shift{0};
    uint64_t m_pos_mask{0};
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> m_ctg_offsets;
    ::bits::compact_vector m_ctg_entries;

    ContigTable()
        : m_ctg_offsets(std::vector<uint64_t>{0}, 1)
    {}

    void set_encoding(uint64_t ref_len_bits) {
        m_ref_len_bits = ref_len_bits;
        m_ref_shift = ref_len_bits + 1;
        m_pos_mask = (ref_len_bits >= 64) ? ~uint64_t(0)
                     : ((uint64_t(1) << ref_len_bits) - 1);
    }

    contig_span entries(uint64_t unitig_id) {
        auto start_pos = m_ctg_offsets.select(unitig_id);
        auto end_pos = m_ctg_offsets.select(unitig_id + 1);
        size_t len = end_pos - start_pos;
        return {m_ctg_entries.at(start_pos),
                m_ctg_entries.at(start_pos + len), len};
    }

    uint64_t encode_entry(uint64_t ref_id, uint64_t pos, bool is_fw) const {
        uint64_t e = ref_id;
        e <<= m_ref_shift;
        e |= (pos << 1);
        e |= is_fw ? 1 : 0;
        return e;
    }

    uint32_t decode_ref_id(uint64_t entry) const {
        return static_cast<uint32_t>(entry >> m_ref_shift);
    }

    uint32_t decode_pos(uint64_t entry) const {
        return static_cast<uint32_t>((entry >> 1) & m_pos_mask);
    }

    bool decode_orientation(uint64_t entry) const {
        return (entry & 0x1);
    }

    template <class Archive>
    void save(Archive& ar) const {
        ar(m_ref_len_bits);
        ar(m_ctg_offsets);
        ar(m_ctg_entries);
    }

    template <class Archive>
    void load(Archive& ar) {
        ar(m_ref_len_bits);
        ar(m_ctg_offsets);
        ar(m_ctg_entries);
        set_encoding(m_ref_len_bits);
    }
};
