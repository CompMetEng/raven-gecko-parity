"""Write proteomics_edge.tsv from the full_ecModel tutorial's proteomics.tsv.

Two conditions, three replicates then two, built from 400 rows of the real data and edited
so that the loaders meet the cases the real file does not contain: a missing value, a zero
(counted as missing), a row with a large spread between replicates, a protein that the
model does not have, an empty identifier and an identifier that occurs twice. Two files are
written: proteomics_edge.tsv marks a missing value as NA, and proteomics_empty.tsv leaves the
cell empty (the last cells of a line empty included).

    python make_proteomics_edge.py <GECKO>/tutorials/full_ecModel/data/proteomics.tsv proteomics_edge.tsv proteomics_empty.tsv
"""

import math
import sys

rows = [line.rstrip("\n").split("\t") for line in open(sys.argv[1], encoding="utf-8")]
header, body = rows[0], rows[1:401]
out = [["uniprotID_S288c", "A_1", "A_2", "A_3", "B_1", "B_2"]]
def scaled(text, factor):
    value = float(text) * factor
    return "NaN" if math.isnan(value) else f"{value:.6E}"


for ident, a, b, c in body:
    out.append([ident, a, b, c, scaled(a, 1.5), scaled(c, 0.5)])

MISSING = "<missing>"
out[1][2] = MISSING                              # one replicate not measured
out[2][1] = "0"                                  # a zero counts as missing
out[3][1:4] = ["1.0E-03", "5.0E+00", "9.0E-02"]  # large spread in condition A
out[4][4:6] = [MISSING, MISSING]                # condition B not measured at all
out[5][1:] = [MISSING, MISSING, MISSING, "2.0E-01", "3.0E-01"]  # condition A not measured at all
out.append(["P00000", "1.0", "1.0", "1.0", "1.0", "1.0"])   # not in the model
out.append(["", "1.0", "1.0", "1.0", "1.0", "1.0"])         # empty identifier
out.append(list(out[10]))                                    # identifier twice

for path, marker in ((sys.argv[2], "NA"), (sys.argv[3], "")):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        for row in out:
            handle.write("\t".join(marker if cell == MISSING else cell for cell in row) + "\n")
