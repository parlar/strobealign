#include "doctest.h"
#include "sam.hpp"
#include "revcomp.hpp"

TEST_CASE("Formatting unmapped SAM record") {
    klibpp::KSeq kseq;
    kseq.name = "read1";
    kseq.seq = "ACGT";
    kseq.qual = ">#BB";
    std::string sam_string;
    References references;

    SUBCASE("without RG") {
        Sam sam(sam_string, references);
        sam.add_unmapped(kseq);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\n");
    }

    SUBCASE("with RG") {
        Sam sam(sam_string, references, CigarOps::EQX, "rg1");
        sam.add_unmapped(kseq);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t>#BB\tRG:Z:rg1\n");
    }

    SUBCASE("without qualities") {
        Sam sam(sam_string, references);
        klibpp::KSeq kseq_noqual;
        kseq_noqual.name = "read1";
        kseq_noqual.seq = "ACGT";
        kseq_noqual.qual = "";
        sam.add_unmapped(kseq_noqual);

        CHECK(sam_string == "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n");
    }

    SUBCASE("do not output single-end unmapped") {
        Sam sam{sam_string, references, CigarOps::EQX, "", false};
        sam.add_unmapped(kseq);
        CHECK(sam_string == "");
    }

    SUBCASE("do not output paired-end unmapped") {
        Sam sam{sam_string, references, CigarOps::EQX, "", false};
        sam.add_unmapped_pair(kseq, kseq);
        CHECK(sam_string == "");
    }
}

TEST_CASE("Sam::add") {
    References references;
    references.add("contig1", "AACCGGTT");

    klibpp::KSeq record;
    record.name = "readname";
    record.seq = "AACGT";
    record.qual = ">#BB";

    Alignment aln;
    aln.ref_id = 0;
    aln.is_unaligned = false;
    aln.is_revcomp = true;
    aln.ref_start = 2;
    aln.edit_distance = 3;
    aln.score = 9;
    aln.cigar = Cigar("2S2=1X3=3S");

    std::string read_rc = reverse_complement(record.seq);
    Details details;
    SUBCASE("Cigar =/X") {
        std::string sam_string;
        Sam sam(sam_string, references);
        sam.add(aln, record, read_rc, 55, PRIMARY, details);
        CHECK(sam_string ==
            "readname\t16\tcontig1\t3\t55\t2S2=1X3=3S\t*\t0\t0\tACGTT\tBB#>\tNM:i:3\tAS:i:9\tms:i:9\tde:f:0.5000\trl:i:5\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0\tXp:i:0\tXs:i:0\tXc:i:0\tXW:f:1.00\tMD:Z:2G3\tnn:i:0\tXG:f:0.67\n"
        );
    }
    SUBCASE("Cigar M") {
        std::string sam_string;
        Sam sam(sam_string, references, CigarOps::M);
        sam.add(aln, record, read_rc, 55, PRIMARY, details);
        CHECK(sam_string ==
            "readname\t16\tcontig1\t3\t55\t2S6M3S\t*\t0\t0\tACGTT\tBB#>\tNM:i:3\tAS:i:9\tms:i:9\tde:f:0.5000\trl:i:5\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0\tXp:i:0\tXs:i:0\tXc:i:0\tXW:f:1.00\tMD:Z:2G3\tnn:i:0\tXG:f:0.67\n"
        );
    }
    SUBCASE("secondary alignment has tp:A:S") {
        std::string sam_string;
        Sam sam(sam_string, references);
        sam.add(aln, record, read_rc, 0, SECONDARY_ALN, details);
        // Secondary: flag 0x110 = SECONDARY|REVERSE, SEQ/QUAL are *, tp:A:S
        CHECK(sam_string.find("\ttp:A:S\t") != std::string::npos);
        CHECK(sam_string.find("\t0x") == std::string::npos); // sanity: no hex in output
    }
    SUBCASE("supplementary alignment has tp:A:S") {
        std::string sam_string;
        Sam sam(sam_string, references);
        sam.add(aln, record, read_rc, 55, SUPPLEMENTARY_ALN, details);
        CHECK(sam_string.find("\ttp:A:S\t") != std::string::npos);
    }
}

TEST_CASE("Pair with one unmapped SAM record") {
    References references;
    references.add("contig1", "ACGT");
    std::string sam_string;
    Sam sam(sam_string, references);

    Alignment aln1;
    aln1.ref_id = 0;
    aln1.is_unaligned = false;
    aln1.ref_start = 2;
    aln1.is_revcomp = true;
    aln1.edit_distance = 17;
    aln1.score = 9;
    aln1.cigar = Cigar("2=");

    Alignment aln2;
    aln2.is_unaligned = true;

    klibpp::KSeq record1;
    klibpp::KSeq record2;
    record1.name = "readname";
    record1.seq = "AACC";
    record1.qual = "#!B<";
    record2.name = "readname";
    record2.seq = "GGTT";
    record2.qual = "IHB#";
    std::string read1_rc = reverse_complement(record1.seq);
    std::string read2_rc = reverse_complement(record2.seq);

    int mapq1 = 55;
    int mapq2 = 57;
    bool is_proper = false;
    std::array<Details, 2> details;

    sam.add_pair(
        aln1,
        aln2,
        record1,
        record2,
        read1_rc,
        read2_rc,
        mapq1,
        mapq2,
        is_proper,
        PRIMARY,
        details
    );
    // 89: PAIRED,MUNMAP,REVERSE,READ1
    // 165: PAIRED,UNMAP,MREVERSE,READ2
    CHECK(sam_string ==
      "readname\t89\tcontig1\t3\t55\t2=\t=\t3\t0\tGGTT\t<B!#\tNM:i:17\tAS:i:9\tms:i:9\tde:f:8.5000\trl:i:4\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0\tXp:i:0\tXs:i:0\tXc:i:0\tXW:f:1.00\tMD:Z:2\tnn:i:0\tXG:f:0.50\n"
      "readname\t165\tcontig1\t3\t0\t*\t=\t3\t0\tGGTT\tIHB#\n"
    );
}

TEST_CASE("TLEN zero when reads map to different contigs") {
    References references;
    references.add("contig1", "ACGT");
    references.add("contig2", "GGAATTCC");
    std::string sam_string;

    Alignment aln1;
    aln1.ref_id = 0;
    aln1.is_unaligned = false;
    aln1.ref_start = 2;
    aln1.is_revcomp = false;
    aln1.edit_distance = 17;
    aln1.score = 9;
    aln1.cigar = Cigar("2=");

    Alignment aln2;
    aln2.is_unaligned = false;
    aln2.ref_id = 1;
    aln2.ref_start = 3;
    aln2.is_revcomp = false;
    aln2.edit_distance = 2;
    aln2.score = 4;
    aln2.cigar = Cigar("3=");

    klibpp::KSeq record1;
    klibpp::KSeq record2;
    record1.name = "readname";
    record1.seq = "AACC";
    record1.qual = "#!B<";
    record2.name = "readname";
    record2.seq = "GGTT";
    record2.qual = "IHB#";
    std::string read1_rc = reverse_complement(record1.seq);
    std::string read2_rc = reverse_complement(record2.seq);

    int mapq1 = 55;
    int mapq2 = 57;
    bool is_proper = false;
    std::array<Details, 2> details;

    Sam sam(sam_string, references);

    sam.add_pair(
        aln1,
        aln2,
        record1,
        record2,
        read1_rc,
        read2_rc,
        mapq1,
        mapq2,
        is_proper,
        PRIMARY,
        details
    );
    // 65: PAIRED,READ1
    // 129: PAIRED,READ2
    CHECK(sam_string ==
    "readname\t65\tcontig1\t3\t55\t2=\tcontig2\t4\t0\tAACC\t#!B<\tNM:i:17\tAS:i:9\tms:i:9\tde:f:8.5000\trl:i:4\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0\tXp:i:0\tXs:i:0\tXc:i:0\tXW:f:1.00\tMD:Z:2\tnn:i:0\tXG:f:0.50\tMC:Z:3=\tMQ:i:57\n"
    "readname\t129\tcontig2\t4\t57\t3=\tcontig1\t3\t0\tGGTT\tIHB#\tNM:i:2\tAS:i:4\tms:i:4\tde:f:0.6667\trl:i:4\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0\tXp:i:0\tXs:i:0\tXc:i:0\tXW:f:1.00\tMD:Z:3\tnn:i:0\tXG:f:0.00\tMC:Z:2=\tMQ:i:55\n"
    );
}

TEST_CASE("MD tag generation") {
    References references;
    references.add("ref", "AACCGGTTAA");

    std::string sam_string;
    Sam sam(sam_string, references);

    SUBCASE("all matches") {
        Cigar cigar("5=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:5");
    }

    SUBCASE("single mismatch") {
        // ref[2]='C'
        Cigar cigar("2=1X2=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:2C2");
    }

    SUBCASE("consecutive mismatches") {
        // ref[1]='A', ref[2]='C'
        Cigar cigar("1=2X2=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:1A0C2");
    }

    SUBCASE("deletion") {
        // ref[2]='C', ref[3]='C'
        Cigar cigar("2=2D2=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:2^CC2");
    }

    SUBCASE("insertion") {
        // Insertions consume query but not reference, invisible in MD
        Cigar cigar("2=2I3=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:5");
    }

    SUBCASE("soft clip") {
        // Soft clips don't consume reference
        Cigar cigar("2S3=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:3");
    }

    SUBCASE("mismatch at start") {
        // ref[0]='A'
        Cigar cigar("1X4=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:0A4");
    }

    SUBCASE("mismatch at end") {
        // ref[4]='G'
        Cigar cigar("4=1X");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:4G0");
    }

    SUBCASE("with ref_start offset") {
        // Starting at ref pos 3: ref[3]='C', ref[4]='G', ref[5]='G'
        Cigar cigar("1X1=1X");
        CHECK(sam.compute_md_tag(cigar, 0, 3) == "MD:Z:0C1G0");
    }

    SUBCASE("deletion after mismatch") {
        // ref[0]='A' (mismatch), then ref[1..2]='AC' (deletion), then ref[3..4]='CG' (match)
        // Zero between mismatch and deletion is required by SAM spec
        Cigar cigar("1X2D2=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:0A0^AC2");
    }

    SUBCASE("complex: soft clip + match + mismatch + deletion + match") {
        // 1S: skip | 3= at ref[0..2]='AAC' | 1X at ref[3]='C' → wrong, ref[3]='G'
        // Actually ref = "AACCGGTTAA"
        // 3= at ref[0..2] → match_count=3
        // 1X at ref[3] → emit 3, emit 'C', ref_pos=4 — wait, ref[3]='C'
        // Actually: ref = A(0) A(1) C(2) C(3) G(4) G(5) T(6) T(7) A(8) A(9)
        // 1S: skip | 3= at ref[0..2]: match=3, ref_pos=3
        // 1X at ref[3]='C': emit "3", emit 'C', ref_pos=4
        // 1D at ref[4]='G': emit "0", emit "^G", ref_pos=5
        // 2= at ref[5..6]: match=2, ref_pos=7
        Cigar cigar("1S3=1X1D2=");
        CHECK(sam.compute_md_tag(cigar, 0, 0) == "MD:Z:3C0^G2");
    }
}

TEST_CASE("MC tag generation") {
    References references;
    references.add("ref", "ACGT");
    std::string sam_string;

    SUBCASE("EQX mode") {
        Sam sam(sam_string, references);
        Cigar cigar("2=1X1=");
        CHECK(sam.compute_mc_tag(cigar) == "MC:Z:2=1X1=");
    }

    SUBCASE("M mode") {
        Sam sam(sam_string, references, CigarOps::M);
        Cigar cigar("2=1X1=");
        CHECK(sam.compute_mc_tag(cigar) == "MC:Z:4M");
    }
}
