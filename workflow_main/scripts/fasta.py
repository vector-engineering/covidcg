# coding: utf-8

"""FASTA helper functions

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""


def read_fasta_file(lines):
    """Read a reference FASTA file. This function is built for many entries in
    the FASTA file but we should be loading in one reference at a time for 
    simplicity's sake. If we need to extract from multiple references then we
    should just run the program multiple times instead.

    Parameters
    ----------
    lines: list of str
        - Output from file.readlines()
    
    Returns
    -------
    entries: dict
        key: Name of the FASTA entry
        value: DNA sequence of the FASTA entry
    """

    entries = {}

    # Read sequences
    cur_entry = ""
    cur_seq = ""
    for i, line in enumerate(lines):
        # Strip whitespace
        line = line.strip()

        # If not the name of an entry, add this line to the current sequence
        # (some FASTA files will have multiple lines per sequence)
        if ">" not in line:

            # Skip if empty
            if not line:
                pass
            else:
                cur_seq = cur_seq + line

                # Force to uppercase
                _cur_seq = list(cur_seq.upper())

                # Throw an error for non canonical bases
                # https://www.bioinformatics.org/sms/iupac.html
                # for char in _cur_seq:
                #     if char not in 'ATCGURYSWKMBDHVN':
                #         error_msg = 'Non canonical base: \"{}\" in sequence {} on line {}.'.format(char, line, i)
                #         raise Exception(error_msg)

                # IUPAC also defines gaps as '-' or '.',
                # but the reference shouldn't have gaps.
                # Maybe I can add this in later...

                # Replace the sequence with the edited one
                cur_seq = "".join(_cur_seq)

        # Start of another entry = end of the previous entry
        if ">" in line or i == (len(lines) - 1):
            # Avoid capturing the first one and pushing an empty sequence
            if cur_entry:
                entry_name = cur_entry.split()[0]
                entry_desc = ""
                if len(cur_entry.split()) > 1:
                    entry_desc = " ".join(cur_entry.split()[1:])

                entries[entry_name] = {}
                entries[entry_name]["name"] = entry_name
                entries[entry_name]["sequence"] = cur_seq
                entries[entry_name]["description"] = entry_desc

            # Clear the entry and sequence
            cur_entry = line[1:]
            cur_seq = ""

    return entries
