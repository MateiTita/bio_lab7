import math

import matplotlib.pyplot as plt

# ---------------------------------------------------------------- enzymes
# name, recognition sequence (5'->3'), cleavage site = number of bases of the
# recognition sequence that stay on the left fragment (top strand).
# All three recognise a 6-bp palindrome and leave a 4-nt 5' overhang.
ENZYMES = [
    {"name": "EcoRI",   "site": "GAATTC", "cut": 1},   # G^AATTC
    {"name": "BamHI",   "site": "GGATCC", "cut": 1},   # G^GATCC
    {"name": "HindIII", "site": "AAGCTT", "cut": 1},   # A^AGCTT
]

# ---------------------------------------------------------------- inputs
SINGLE_SEQUENCE = ("CXCL8 mRNA", "data/cxcl8/cxcl8.fasta")          # task 1-2
VIRAL_GENOMES = [                                                    # task 3
    ("phiX174", "data/phix174/phix174.fasta"),
    ("Ebola",   "data/ebola/ebola.fna"),
]

LADDER = [100, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000,
          5000, 7000, 10000, 20000]

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def read_fasta(path):
    with open(path) as f:
        lines = f.read().splitlines()
    header = lines[0][1:]
    sequence = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))
    return header, sequence.upper()


def cut_positions(sequence, enzyme):
    """Positions where the top strand is cut (number of bases left of the cut)."""
    site = enzyme["site"]
    positions = []
    start = sequence.find(site)
    while start != -1:
        positions.append(start + enzyme["cut"])
        start = sequence.find(site, start + 1)
    return positions


def digest(sequence, enzymes):
    """Cut a linear sequence with one or more enzymes -> (cuts, fragment lengths)."""
    cuts = sorted(set(p for e in enzymes for p in cut_positions(sequence, e)))
    edges = [0] + cuts + [len(sequence)]
    fragments = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
    return cuts, fragments


def print_enzyme_table():
    print("Restriction enzymes")
    print(f"  {'name':<8} {'recognition sequence':<22} cleavage site")
    for e in ENZYMES:
        site, cut = e["site"], e["cut"]
        bottom = "".join(COMPLEMENT[b] for b in site)
        bottom_cut = len(site) - cut          # palindrome: mirrored on the other strand
        print(f"  {e['name']:<8} 5'-{site}-3'" + " " * 11
              + f"5'-{site[:cut]}^{site[cut:]}-3'")
        print(f"  {'':<8} 3'-{bottom}-5'" + " " * 11
              + f"3'-{bottom[:bottom_cut]}^{bottom[bottom_cut:]}-5'")
    print()


def build_lanes(label, sequence):
    """One lane per enzyme plus one with all three together."""
    lanes = []
    for enzymes in [[e] for e in ENZYMES] + [ENZYMES]:
        name = enzymes[0]["name"] if len(enzymes) == 1 else "All 3"
        cuts, fragments = digest(sequence, enzymes)
        lanes.append({"name": name, "group": label, "fragments": fragments})
        print(f"  {name:<8} {len(cuts)} cut(s) at {cuts if cuts else '-'}")
        print(f"  {'':<8} fragments (bp): {sorted(fragments, reverse=True)}")
    return lanes


def draw_gel(lanes, title, filename):
    sizes = [s for lane in lanes for s in lane["fragments"]]
    low, high = min(sizes) * 0.7, max(sizes) * 1.6
    ladder = [s for s in LADDER if low <= s <= max(sizes) * 1.25]
    lanes = [{"name": "Ladder", "group": "", "fragments": ladder}] + lanes

    fig, ax = plt.subplots(figsize=(1.6 + 1.25 * len(lanes), 9))
    ax.set_facecolor("black")
    ax.set_yscale("log")            # migration distance ~ -log(size): big on top
    ax.set_ylim(low, high)
    min_gap = (math.log10(high) - math.log10(low)) / 60   # label spacing

    for x, lane in enumerate(lanes, 1):
        is_ladder = lane["name"] == "Ladder"
        ax.hlines(high / 1.08, x - 0.3, x + 0.3, colors="dimgray", linewidth=5)  # well
        last_label = None
        for size in sorted(set(lane["fragments"]), reverse=True):
            ax.hlines(size, x - 0.3, x + 0.3, linewidth=4,
                      colors="white" if is_ladder else "lime")
            # skip a label when it would sit on top of the previous one
            if last_label is None or math.log10(last_label) - math.log10(size) > min_gap:
                ax.text(x + 0.33, size, str(size), color="white", va="center", fontsize=7)
                last_label = size

    names = [lane["name"] if not lane["group"] else f"{lane['group']}\n{lane['name']}"
             for lane in lanes]
    ax.set_xticks(range(1, len(lanes) + 1), names, fontsize=9)
    ax.set_xlim(0.4, len(lanes) + 0.9)
    ax.set_yticks(ladder, [str(s) for s in ladder])
    ax.minorticks_off()
    ax.set_ylabel("fragment length (bp), log scale  -  top: long, bottom: short")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(filename, dpi=150)
    print(f"  saved {filename}\n")


print_enzyme_table()

# ---- tasks 1-2: one sequence of 1000-3000 nt
label, path = SINGLE_SEQUENCE
header, sequence = read_fasta(path)
print(f"{header}\n  length: {len(sequence)} bp")
draw_gel(build_lanes(label, sequence),
         f"Restriction digest of {label} ({len(sequence)} bp)", "gel_single.png")

# ---- task 3: two viral genomes side by side
lanes = []
for label, path in VIRAL_GENOMES:
    header, sequence = read_fasta(path)
    print(f"{header}\n  length: {len(sequence)} bp")
    lanes += build_lanes(f"{label}\n{len(sequence)} bp", sequence)
draw_gel(lanes, "Restriction digest: " + " vs ".join(name for name, _ in VIRAL_GENOMES),
         "gel_viruses.png")

plt.show()
