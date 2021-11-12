# coding: utf-8

"""Extract variable regions from an aligned segment, in a flexible
and SNP-tolerant manner

Modified and heavily trimmed down version of read_extractor.py (v0.1.0)
from the variant_extractor project

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import numpy as np
import pandas as pd

from collections import defaultdict

from scripts.util import translate, reverse_complement


class ReadExtractor:
    """Extract variable regions from a pysam AlignedSegment
    """

    RefSeq = ""

    def __init__(self, read):
        """Build the extactor object for a read (pysam.AlignedSegment)
        or a pair of reads if using paired-end sequencing

        Parameters
        ----------
        read: pysam.AlignedSegment
        """

        self.read = read

        # Build our own mutation string to store mutational information
        # Since both the CIGAR and MD string don't fit our needs
        # Format: Position:Ref:Alt;...
        # Where position is relative to the reference (0-indexed)
        # For insertions, the position is the position on the reference
        # after the insertion
        # For deletions, the position is the position on the reference
        # that was deleted
        # Store it as a list of tuples, (Position, Ref, Alt) for now.
        # Mutations will be individually serialized then joined by ';' later
        # to serialize into one big string
        self.mutation_str = []

        # Any invalidation errors that flag this variant as successfully extracted,
        # but not passing filters, will be stored in this array
        # Later when writing to disk we'll serialize this array as a semicolon-delimited string
        self.invalid_errors = []

        # Store SNPs
        self.dna_snps = []

        # Read data from the pysam.AlignedSegment object into python variables
        self.load_read()

    def load_read(self):
        """Load data in from the pysam.AlignedSegment object into Python
        """

        # Nucleotide sequence of the read
        self.read_seq = self.read.get_forward_sequence()

        # Don't try to do anything else if this read has no sequence
        if not self.read_seq:
            return

        # If reverse complement, flip the sequence and the quality scores
        if self.read.is_reverse:
            self.read_seq = reverse_complement(self.read_seq)

        # Don't try to do anything else if this read is unmapped
        if self.read.is_unmapped:
            return

        # Get the reference sequence
        self.reference_seq = ReadExtractor.RefSeq[self.read.reference_name]

        """Expand CIGAR tuples to a list of CIGAR operations on the read (query)

        https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        https://drive5.com/usearch/manual/cigar.html
        https://samtools.github.io/hts-specs/SAMv1.pdf

        Op                  Code    Description
        -----------------------------------------------------------------------------------------
        M	BAM_CMATCH      0       Match (alignment column containing two letters). This could
                                    contain two different letters (mismatch) or two identical
                                    letters. USEARCH generates CIGAR strings containing Ms rather
                                    than X's and ='s (see below).
        I	BAM_CINS        1       Insertion (gap in the query sequence).
        D	BAM_CDEL        2       Deletion (gap in the target sequence).
        N	BAM_CREF_SKIP   3       skipped region from the reference
        S	BAM_CSOFT_CLIP  4       Segment of the query sequence that does not appear in the
                                    alignment. This is used with soft clipping, where the
                                    full-length query sequence is given (field 10 in the SAM record)
                                    . In this case, S operations specify segments at the start and/
                                    or end of the query that do not appear in a local alignment.
        H	BAM_CHARD_CLIP  5       Segment of the query sequence that does not appear in the
                                    alignment. This is used with hard clipping, where only the
                                    aligned segment of the query sequences is given (field 10 in
                                    the SAM record). In this case, H operations specify segments at
                                    the start and/or end of the query that do not appear in the SAM
                                    record.
        P	BAM_CPAD        6       padding (silent deletion from padded reference)
        =	BAM_CEQUAL      7       Alignment column containing two identical letters. USEARCH can
                                    read CIGAR strings using this operation, but does not generate
                                    them.
        X	BAM_CDIFF       8       Alignment column containing a mismatch, i.e. two different
                                    letters. USEARCH can read CIGAR strings using this operation,
                                    but does not generate them.
        B	BAM_CBACK       9
        """

        self.cigar_ops = []
        for op_group in self.read.cigartuples:
            # First element of the tuple is the operation code
            # Second element of the tuple is the number of operations

            # Create a new list [# of operations] long and add it to the
            # master operations list
            self.cigar_ops.extend([op_group[0],] * op_group[1])

        # Reset the cigar index
        self.cigar_i = 0

        # Start the reference at the position it is mapped onto the read
        # using read.reference_start
        self.ref_i = self.read.reference_start

        # Start the read at the position it is mapped onto the reference
        # using read.query_alignment_start
        self.read_i = 0

    def crawl_to(self, destination):
        """Iterate (consume bases) through both the read and the reference
        Use the CIGAR operations and other stats to stay on the same
        "aligned" base (as if we did a multiple sequence alignment on the read and ref)

        Parameters
        ----------
        destination: int
            - Index on the reference of where we want to crawl to
        """

        while self.ref_i < destination:

            # If we've reached the end of the CIGAR string, break out
            if self.cigar_i >= len(self.cigar_ops) or self.read_i >= len(self.read_seq):
                return

            # Grab the current CIGAR operation
            op = self.cigar_ops[self.cigar_i]

            """
            https://samtools.github.io/hts-specs/SAMv1.pdf
            ---------------------------------------------------
            | Op  | Code  | Consume Read  | Consume Reference |
            ---------------------------------------------------
            | M   | 0     | Yes           | Yes               |
            | I   | 1     | Yes           | No                |
            | D   | 2     | No            | Yes               |
            | N   | 3     | No            | Yes               |
            | S   | 4     | Yes           | No                |
            | H   | 5     | No            | No                |
            | P   | 6     | No            | No                |
            | =   | 7     | Yes           | Yes               |
            | X   | 8     | Yes           | Yes               |
            | B   | 9     | ?             | ?                 |
            ---------------------------------------------------
            """

            # MATCH - can be match or mismatch (SNP)
            if op == 0 or op == 7 or op == 8:

                # Check for SNPs
                # If the OP code is 0, then we have to check both the read
                # and the reference to see if there's a mismatch
                # If bowtie2 gave us the OP code of 8, then we know there's a mismatch
                if (
                    # Check for a mismatch OP code or a base mismatch for a
                    # generic 0 OP code
                    (
                        (op == 8)
                        or (
                            op == 0
                            and self.read_seq[self.read_i]
                            != self.reference_seq[self.ref_i]
                        )
                    )
                    and
                    # If the reference has an X as the base, then
                    # ignore any SNPs at this position
                    (self.reference_seq[self.ref_i] != "X")
                ):
                    # Add substitution information to mutation string
                    self.mutation_str.append(
                        (
                            self.read.query_name,
                            self.ref_i,
                            self.reference_seq[self.ref_i],
                            self.read_seq[self.read_i],
                        )
                    )

                self.read_i += 1
                self.ref_i += 1

            # Insertion or Soft Clip
            elif op == 1 or op == 4:

                # Add insertion information to mutation string
                self.mutation_str.append(
                    (self.read.query_name, self.ref_i, "", self.read_seq[self.read_i])
                )

                self.read_i += 1

            # Deletion or Skip
            elif op == 2 or op == 3:

                # Add deletion information to mutation string
                self.mutation_str.append(
                    (
                        self.read.query_name,
                        self.ref_i,
                        self.reference_seq[self.ref_i],
                        "",
                    )
                )

                self.ref_i += 1

            # Hard Clip, Padding
            else:
                # Do nothing
                pass

            # Always iterate the CIGAR index
            self.cigar_i += 1

        # END WHILE

    def get_dna_snps(self):
        """Store list of NT SNPs/indels"""

        # Join adjacent indels
        self.dna_snps = []
        i = 0
        while i < len(self.mutation_str):
            (query_name, pos, ref, alt) = self.mutation_str[i]
            # mut is a tuple: (Position, Ref, Alt)

            # Offset the position back to 1-indexed, starting at the genome start
            pos = pos + 1

            # If it's a SNP, then add and continue
            if ref and alt:
                i += 1

                # Actually, skip adding it if either the ref or the alt
                # is an ambiguous base (N)
                # This is useless data bloat and should be removed as
                # early as possible
                if alt not in ["A", "C", "G", "T"]:
                    continue

                self.dna_snps.append((query_name, pos, ref, alt, self.read.reference_name))
                continue

            # Check ahead for adjacent positions and the same indel type
            j = i
            while j < len(self.mutation_str) and (
                # Both insertions
                (
                    (not self.mutation_str[j][2] and not ref)
                    # Both deletions
                    or (not self.mutation_str[j][3] and not alt)
                )
                # New position must be adjacent to the previous one
                and self.mutation_str[j][1] == int(pos - 1 + (j - i))
            ):
                j += 1

            # Get adjacent indels
            adj_muts = self.mutation_str[i:j]
            # Combine bases, but keep first position and type
            self.dna_snps.append(
                (
                    query_name,
                    pos,
                    "".join([m[2] for m in adj_muts]),
                    "".join([m[3] for m in adj_muts]),
                    self.read.reference_name
                )
            )
            # Skip ahead to the end of the adjacent mutations
            i = j

    def process_all(self):
        """Do everything, return everything"""
        # Skip read without read_seq
        if not self.read_seq:
            return None

        # Travel to the end of the read
        # so that we can collect additional mutations (if they exist)
        # Don't throw an error once we reach the end
        self.crawl_to(len(self.reference_seq))
        self.get_dna_snps()
        return self.dna_snps
