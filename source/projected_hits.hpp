#pragma once

#include <cstdint>
#include <iostream>
#include "contig_table.hpp"

struct ref_pos {
    uint32_t pos;
    bool isFW;
};

struct projected_hits {
    uint32_t contigIdx_;
    uint32_t contigPos_;
    bool contigOrientation_;
    uint32_t contigLen_;
    uint32_t k_;
    contig_span refRange;

    bool empty() const { return refRange.empty(); }

    uint32_t contig_id() const { return contigIdx_; }
    uint32_t contig_pos() const { return contigPos_; }
    uint32_t contig_len() const { return contigLen_; }
    bool hit_fw_on_contig() const { return contigOrientation_; }

    ref_pos decode_hit(uint64_t v, const ContigTable& ctab) const {
        bool contigFW = ctab.decode_orientation(v);
        uint32_t rpos;
        bool rfw;
        if (contigFW and contigOrientation_) {
            rpos = ctab.decode_pos(v) + contigPos_;
            rfw = true;
        } else if (contigFW and !contigOrientation_) {
            rpos = ctab.decode_pos(v) + contigPos_;
            rfw = false;
        } else if (!contigFW and contigOrientation_) {
            rpos = ctab.decode_pos(v) + contigLen_ - (contigPos_ + k_);
            rfw = false;
        } else {
            rpos = ctab.decode_pos(v) + contigLen_ - (contigPos_ + k_);
            rfw = true;
        }
        return {rpos, rfw};
    }

    friend std::ostream& operator<<(std::ostream& os, const projected_hits& h) {
        os << "{ contig_id:" << h.contigIdx_
           << " pos:" << h.contigPos_
           << " ori:" << (h.contigOrientation_ ? "fw" : "rc")
           << " len:" << h.contigLen_
           << " refs:" << h.refRange.size() << " }";
        return os;
    }
};
