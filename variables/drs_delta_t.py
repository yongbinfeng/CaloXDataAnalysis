"""
DRS same-row Delta-T variables for both 6 mm and nominal 3 mm channels.

Defines event-by-event:

    DeltaT = t_right - t_left

for same-row channel combinations in:

    Quartz  6 mm
    Plastic 6 mm
    Quartz  3 mm nominal
    Plastic 3 mm nominal

The timing input is the already-corrected DRS column:

    {channel}_Time_cfd_mcp

which is CFD + reference-channel corrected + MCP corrected + group-MPV corrected
and converted to ns by variables/drs.py.
"""

from itertools import combinations


# ============================================================
# Full detector grids
# ============================================================

QUARTZ_GRID = [
    [None,  "002", None,  None ],
    ["006", "004", "206", "204"],
    ["016", "014", "216", "214"],
    ["026", "024", "226", "224"],
    [None,  "030", None,  None ],
    [None,  "034", None,  None ],
    ["106", "104", "306", "304"],
    ["116", "114", "316", "314"],
    ["126", "124", "326", "324"],
    [None,  "134", None,  "334"],
]

PLASTIC_GRID = [
    [None,  "000", "202", "200"],
    ["012", "010", "212", "210"],
    ["022", "020", "222", "220"],
    ["032", None,  "232", "230"],
    ["102", "100", "302", "300"],
    ["112", "110", "312", "310"],
    ["122", "120", "322", "320"],
    ["132", "130", "332", "330"],
]


# ============================================================
# Split into 6 mm and nominal 3 mm sections
# ============================================================
# Keep original row numbers so the HTML output matches the physical grid.

DETECTOR_BLOCKS = [
    {
        "family": "Quartz6mm",
        "family_display": "6 mm FSHA Quartz",
        "material": "Quartz",
        "size": "6mm",
        "grid": QUARTZ_GRID,
        "rows": range(0, 6),
    },
    {
        "family": "Plastic6mm",
        "family_display": "6 mm Toray PJR-FB750 Plastic",
        "material": "Plastic",
        "size": "6mm",
        "grid": PLASTIC_GRID,
        "rows": range(0, 4),
    },
    {
        "family": "Quartz3mm",
        "family_display": "3 mm FSHA Quartz",
        "material": "Quartz",
        "size": "3mm",
        "grid": QUARTZ_GRID,
        "rows": range(6, 10),
    },
    {
        "family": "Plastic3mm",
        "family_display": "3 mm Toray PJR-FB750 Plastic",
        "material": "Plastic",
        "size": "3mm",
        "grid": PLASTIC_GRID,
        "rows": range(4, 8),
    },
]


def code_to_channel(code):
    """
    Convert compact code like '106' into full DRS branch name:

        '106' -> 'DRS_Board1_Group0_Channel6'
    """
    return f"DRS_Board{int(code[0])}_Group{int(code[1])}_Channel{int(code[2])}"


def build_block_pairs(block, require_complete_rows=False):
    """
    Build all same-row C(N,2) combinations.

    For complete 4-channel rows:
        C(4,2) = 6 pairs

    For partial rows:
        C(3,2) = 3 pairs
        C(2,2) = 1 pair

    We keep partial rows because some real 6 mm / 3 mm rows are not fully populated.
    """
    pairs = []

    grid = block["grid"]

    for irow in block["rows"]:
        row = grid[irow]
        row4 = list(row) + [None] * (4 - len(row))
        valid = [ch for ch in row4 if ch is not None]

        if require_complete_rows and len(valid) != 4:
            continue

        if len(valid) < 2:
            continue

        for left, right in combinations(valid, 2):
            pairs.append({
                "family": block["family"],
                "family_display": block["family_display"],
                "material": block["material"],
                "size": block["size"],
                "row": irow,
                "row_channels": valid,
                "left": left,
                "right": right,
            })

    return pairs


DELTA_T_PAIRS = []
for block in DETECTOR_BLOCKS:
    DELTA_T_PAIRS.extend(build_block_pairs(block, require_complete_rows=False))


def dt_col_name(family, irow, left, right):
    return f"DeltaT_{family}_row{irow}_{left}_{right}"


def good_col_name(family, irow, left, right):
    return f"GoodDeltaT_{family}_row{irow}_{left}_{right}"


def define_drs_delta_t_pairs(ctx, amp_threshold=100.0):
    """
    Define corrected Delta-T columns on ctx.rdf.

    Uses:
        {channel}_Time_cfd_mcp

    Also defines a boolean good-column per pair requiring:
        finite left/right times
        left/right peak_value > amp_threshold
    """
    rdf = ctx.rdf
    available = set(str(c) for c in rdf.GetColumnNames())

    defined_pairs = []

    for p in DELTA_T_PAIRS:
        family = p["family"]
        irow = p["row"]
        left = p["left"]
        right = p["right"]

        ch_left = code_to_channel(left)
        ch_right = code_to_channel(right)

        t_left = f"{ch_left}_Time_cfd_mcp"
        t_right = f"{ch_right}_Time_cfd_mcp"

        amp_left = f"{ch_left}_peak_value"
        amp_right = f"{ch_right}_peak_value"

        needed = [t_left, t_right, amp_left, amp_right]
        missing = [x for x in needed if x not in available]

        if missing:
            print(
                f"[DeltaT] skipping {family} row {irow} {left}-{right}: "
                f"missing {missing}",
                flush=True,
            )
            continue

        dt_col = dt_col_name(family, irow, left, right)
        good_col = good_col_name(family, irow, left, right)

        # Delta-T convention:
        #   Δt = t_right - t_left
        rdf = rdf.Define(dt_col, f"{t_right} - {t_left}")

        good_expr = (
            f"std::isfinite({t_left}) && "
            f"std::isfinite({t_right}) && "
            f"{amp_left} > {amp_threshold} && "
            f"{amp_right} > {amp_threshold}"
        )

        rdf = rdf.Define(good_col, good_expr)

        available.add(dt_col)
        available.add(good_col)

        out = dict(p)
        out.update({
            "left_channel": ch_left,
            "right_channel": ch_right,
            "dt_col": dt_col,
            "good_col": good_col,
        })

        defined_pairs.append(out)

    ctx.delta_t_pair_defs = defined_pairs

    print(
        f"[DeltaT] defined {len(defined_pairs)} Delta-T pair columns "
        f"including 6 mm and 3 mm nominal channels",
        flush=True,
    )

    return rdf
