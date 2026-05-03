#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <gtl/phmap.hpp>
#include "rshash.hpp"
#include "contig_table.hpp"

struct rank_count {
    uint64_t rank;
    uint32_t len;
    uint64_t count;
};

void build_rshash_from_cf_seg(
    const std::string& cf_basename,
    uint64_t k, uint8_t level,
    uint8_t m1, uint8_t m2, uint8_t m3,
    uint8_t t1, uint8_t t2, uint16_t t3,
    RSHash& index,
    gtl::flat_hash_map<uint64_t, rank_count>& id_to_rank,
    std::vector<uint64_t>& segment_order)
{
    std::ifstream seg_file(cf_basename + ".cf_seg");
    if (!seg_file.is_open()) {
        std::cerr << "error: cannot open " << cf_basename << ".cf_seg\n";
        std::exit(1);
    }

    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> sequences;
    uint64_t idx = 0;
    uint64_t seg_id;
    std::string seg_seq;
    while (seg_file >> seg_id >> seg_seq) {
        segment_order.push_back(seg_id);
        id_to_rank[seg_id] = {idx, static_cast<uint32_t>(seg_seq.length()), 0};

        seqan3::bitpacked_sequence<seqan3::dna4> bps;
        bps.resize(seg_seq.size());
        for (size_t i = 0; i < seg_seq.size(); ++i) {
            switch (seg_seq[i]) {
                case 'A': case 'a': bps[i] = seqan3::dna4{}.assign_char('A'); break;
                case 'C': case 'c': bps[i] = seqan3::dna4{}.assign_char('C'); break;
                case 'G': case 'g': bps[i] = seqan3::dna4{}.assign_char('G'); break;
                case 'T': case 't': bps[i] = seqan3::dna4{}.assign_char('T'); break;
                default: bps[i] = seqan3::dna4{}.assign_char('A'); break;
            }
        }
        sequences.push_back(std::move(bps));
        ++idx;
    }

    std::cout << "parsed " << idx << " segments from .cf_seg\n";

    index = RSHash(k, level, m1, m2, m3, t1, t2, t3);
    index.build(sequences);
}

void build_contig_table_from_cf_seq(
    const std::string& cf_basename,
    uint64_t k,
    gtl::flat_hash_map<uint64_t, rank_count>& id_to_rank,
    const std::vector<uint64_t>& segment_order,
    ContigTable& ctab,
    std::vector<std::string>& ref_names,
    std::vector<uint64_t>& ref_lens)
{
    const std::string refstr = "Reference";
    const auto hlen = refstr.length();

    // First pass: count occurrences, compute reference lengths
    size_t num_refs = 0;
    size_t max_ref_len = 0;
    {
        std::ifstream ifile(cf_basename + ".cf_seq");
        if (!ifile.is_open()) {
            std::cerr << "error: cannot open " << cf_basename << ".cf_seq\n";
            std::exit(1);
        }

        uint64_t refctr = 0;
        bool first = true;
        uint64_t current_offset = 0;

        std::string tok;
        while (ifile >> tok) {
            if (tok.compare(0, hlen, refstr) == 0) {
                auto sp = tok.find("Sequence:");
                auto ep = sp + 9;

                if (!first) {
                    uint64_t len = current_offset + (k - 1);
                    ref_lens.push_back(len);
                    max_ref_len = std::max(static_cast<uint64_t>(max_ref_len), len);
                    ++refctr;
                }

                ref_names.push_back(tok.substr(ep));
                current_offset = 0;
                first = false;
            } else {
                bool is_n_tile = false;
                if (!((tok.back() == '-') or (tok.back() == '+'))) {
                    if (tok.front() == 'N') {
                        is_n_tile = true;
                    } else {
                        std::cerr << "error: unexpected tiling entry: " << tok << "\n";
                        std::exit(1);
                    }
                }

                if (is_n_tile) {
                    tok.erase(0, 1);
                    uint64_t num_ns = std::stoul(tok, nullptr, 10);
                    if (current_offset > 0) {
                        current_offset += (k - 1);
                    }
                    current_offset += num_ns;
                } else {
                    tok.pop_back();
                    uint64_t id = std::stoul(tok, nullptr, 10);
                    auto rit = id_to_rank.find(id);
                    if (rit == id_to_rank.end()) {
                        std::cerr << "error: segment " << id << " not found in .cf_seg\n";
                        std::exit(1);
                    }
                    rit->second.count += 1;
                    current_offset += rit->second.len - (k - 1);
                }
            }
        }

        if (!first) {
            uint64_t len = current_offset + (k - 1);
            ref_lens.push_back(len);
            max_ref_len = std::max(static_cast<uint64_t>(max_ref_len), len);
        }

        num_refs = ref_lens.size();
    }

    uint64_t ref_len_bits = static_cast<uint64_t>(std::ceil(std::log2(max_ref_len + 1)));
    uint64_t num_ref_bits = static_cast<uint64_t>(std::ceil(std::log2(num_refs + 1)));
    uint64_t total_ctg_bits = ref_len_bits + num_ref_bits + 1;

    ctab.set_encoding(ref_len_bits);

    std::cout << "references: " << num_refs << ", segments: " << id_to_rank.size() << "\n";
    std::cout << "max ref len: " << max_ref_len << " (" << ref_len_bits << " bits)\n";

    // Build cumulative offset vector
    uint64_t tot_seg_occ = 0;
    {
        std::vector<uint64_t> contig_offsets;
        contig_offsets.reserve(segment_order.size() + 1);
        contig_offsets.push_back(0);
        for (auto seg_id : segment_order) {
            tot_seg_occ += id_to_rank[seg_id].count;
            contig_offsets.push_back(tot_seg_occ);
        }

        // Convert counts to running offsets for second pass
        for (size_t i = 0; i < segment_order.size(); ++i) {
            auto seg_id = segment_order[i];
            id_to_rank[seg_id].count = contig_offsets[i];
        }

        uint64_t ef_universe = contig_offsets.back() + 1;
        ctab.m_ctg_offsets = sux::bits::EliasFano<sux::util::AllocType::MALLOC>(
            contig_offsets, ef_universe);
    }

    std::cout << "total segment occurrences: " << tot_seg_occ << "\n";

    // Second pass: fill contig entries
    {
        auto builder = ::bits::compact_vector::builder(tot_seg_occ, total_ctg_bits);

        std::ifstream ifile(cf_basename + ".cf_seq");
        uint64_t refctr = 0;
        bool first = true;
        uint64_t current_offset = 0;

        std::string tok;
        while (ifile >> tok) {
            if (tok.compare(0, hlen, refstr) == 0) {
                if (!first) { ++refctr; }
                first = false;
                current_offset = 0;
            } else {
                bool is_n_tile = false;
                bool is_fw = true;
                if (tok.back() == '-') {
                    is_fw = false;
                } else if (tok.back() == '+') {
                    is_fw = true;
                } else if (tok.front() == 'N') {
                    is_n_tile = true;
                } else {
                    std::cerr << "error: unexpected tiling entry: " << tok << "\n";
                    std::exit(1);
                }

                if (is_n_tile) {
                    tok.erase(0, 1);
                    uint64_t num_ns = std::stoul(tok, nullptr, 10);
                    if (current_offset > 0) {
                        current_offset += (k - 1);
                    }
                    current_offset += num_ns;
                } else {
                    tok.pop_back();
                    uint64_t id = std::stoul(tok, nullptr, 10);
                    auto& v = id_to_rank[id];
                    auto entry_idx = v.count;
                    uint64_t encoded = ctab.encode_entry(refctr, current_offset, is_fw);
                    builder.set(entry_idx, encoded);
                    v.count += 1;
                    current_offset += v.len - (k - 1);
                }
            }
        }

        builder.build(ctab.m_ctg_entries);
    }
}

void build_reference_index(
    const std::string& cf_basename,
    const std::string& output_basename,
    uint64_t k, uint8_t level,
    uint8_t m1, uint8_t m2, uint8_t m3,
    uint8_t t1, uint8_t t2, uint16_t t3)
{
    RSHash index;
    gtl::flat_hash_map<uint64_t, rank_count> id_to_rank;
    std::vector<uint64_t> segment_order;

    std::cout << "building RSHash dictionary from .cf_seg...\n";
    build_rshash_from_cf_seg(cf_basename, k, level, m1, m2, m3, t1, t2, t3,
                             index, id_to_rank, segment_order);

    std::cout << "saving RSHash dictionary...\n";
    index.save(output_basename + ".rshash");

    ContigTable ctab;
    std::vector<std::string> ref_names;
    std::vector<uint64_t> ref_lens;

    std::cout << "building contig table from .cf_seq...\n";
    build_contig_table_from_cf_seq(cf_basename, k, id_to_rank, segment_order,
                                   ctab, ref_names, ref_lens);

    std::cout << "saving contig table...\n";
    {
        std::ofstream out(output_basename + ".ctab", std::ios::binary);
        cereal::BinaryOutputArchive ar(out);
        ar(ctab);
    }

    std::cout << "saving reference info...\n";
    {
        std::ofstream out(output_basename + ".refinfo", std::ios::binary);
        cereal::BinaryOutputArchive ar(out);
        ar(ref_names, ref_lens);
    }

    std::cout << "reference index built successfully.\n";
    std::cout << "  segments: " << segment_order.size() << "\n";
    std::cout << "  references: " << ref_names.size() << "\n";
}
