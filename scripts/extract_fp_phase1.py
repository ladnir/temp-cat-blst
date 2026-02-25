import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "avxbls" / "src" / "fp12_avx.c"

TARGETS = [
    "add_fp_8x1w",
    "sub_fp_8x1w",
    "asx4_fp_8x1w",
    "asax6_fp_8x1w",
    "assx6_fp_8x1w",
    "add_fp_4x2w",
    "sub_fp_4x2w",
    "asx2_fp_4x2w",
    "add_fpx2_8x1w",
    "sub_fpx2_8x1w",
    "aaasx2_fpx2_8x1w",
    "aassx2_fpx2_8x1w",
    "asx4_fpx2_8x1w",
    "ssaax2_fpx2_8x1w",
    "ssasx2_fpx2_8x1w",
    "add_fpx2_4x2w",
    "sub_fpx2_4x2w",
    "asx2_fpx2_4x2w",
    "aass_fpx2_4x2w",
    "mul_fpx2_8x1w_v1",
    "mul_fpx2_8x1w_v2",
    "mul_fpx2_8x1w_v3",
    "mul_fpx2_8x1w_v4",
    "mul_fpx2_4x2w_v1",
    "redc_fpx2_8x1w",
    "redc_fpx2_4x2w",
]


def find_block(text: str, name: str) -> tuple[int, int]:
    # Find function signature start.
    m = re.search(rf"\b{name}\s*\(", text)
    if not m:
        return (-1, -1)
    # Walk back to line start.
    start = text.rfind("\n", 0, m.start()) + 1
    # Find opening brace for function.
    brace = text.find("{", m.end())
    if brace < 0:
        return (-1, -1)
    depth = 0
    i = brace
    while i < len(text):
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                return (start, end)
        i += 1
    return (-1, -1)


def line_no(text: str, idx: int) -> int:
    return text.count("\n", 0, idx) + 1


def main() -> int:
    text = SRC.read_text(encoding="utf-8", errors="ignore")
    for name in TARGETS:
        s, e = find_block(text, name)
        if s < 0:
            print(f"MISSING {name}")
            continue
        print(f"\n=== {name} [{line_no(text,s)}-{line_no(text,e)}] ===")
        print(text[s:e].strip())
    return 0


if __name__ == "__main__":
    sys.exit(main())
