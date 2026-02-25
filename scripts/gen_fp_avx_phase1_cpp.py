import pathlib
import re


ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC_C = ROOT / "avxbls" / "src" / "fp12_avx.c"
SRC_H = ROOT / "avxbls" / "src" / "fp12_avx.h"
OUT_CPP = ROOT / "blst" / "fields" / "Fp_avx.cpp"


H_FUNCS = [
    "conv_64to48_fp",
    "conv_64to48_fp_8x1w",
    "conv_64to48_fp_4x2w",
    "conv_64to48_mpi",
    "conv_48to64_fp",
    "conv_48to64_fp_8x1w",
    "conv_48to64_fp_4x2w",
    "conv_48to64_mpi",
    "carryp_mpi48",
    "carryp_fp_4x2w",
    "get_channel_8x1w",
    "get_channel_4x2w",
]

C_FUNCS = [
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
    "mul_fpx2_8x1w_v2",
    "mul_fpx2_4x2w_v1",
    "redc_fpx2_8x1w",
    "redc_fpx2_4x2w",
]


def find_block(text: str, name: str) -> str:
    m = re.search(rf"\b{name}\s*\(", text)
    if not m:
        raise RuntimeError(f"missing function {name}")
    start = text.rfind("\n", 0, m.start()) + 1
    brace = text.find("{", m.end())
    if brace < 0:
        raise RuntimeError(f"missing body {name}")
    depth = 0
    i = brace
    while i < len(text):
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return text[start:i + 1]
        i += 1
    raise RuntimeError(f"unterminated {name}")


def normalize(code: str) -> str:
    code = code.replace("__m512i", "cc512i")
    code = code.replace("uint64_t", "u64")
    # Export functions instead of static linkage.
    code = re.sub(r"^\s*static\s+void\s+", "void ", code, flags=re.M)
    return code


def main() -> None:
    src_h = SRC_H.read_text(encoding="utf-8", errors="ignore")
    src_c = SRC_C.read_text(encoding="utf-8", errors="ignore")

    header = []
    header.append('#include "blst/fields/Fp_avx.h"')
    header.append("")
    header.append("#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)")
    header.append("")
    header.append("namespace blst::avx")
    header.append("{")
    header.append("const u64 P48[NWORDS] = {")
    header.append("    0xFFFFFFFFAAABULL, 0xB153FFFFB9FEULL, 0xF6241EABFFFEULL, 0x6730D2A0F6B0ULL,")
    header.append("    0x4B84F38512BFULL, 0x434BACD76477ULL, 0xE69A4B1BA7B6ULL, 0x1A0111EA397FULL,")
    header.append("};")
    header.append("")

    body = []
    for n in H_FUNCS:
        body.append(normalize(find_block(src_h, n)))
        body.append("")

    for n in C_FUNCS:
        body.append(normalize(find_block(src_c, n)))
        body.append("")

    # Phase-1 keeps all variants available; v1/v3/v4 are routed to v2 for now.
    body.append("void mul_fpx2_8x1w_v1(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b)")
    body.append("{")
    body.append("    mul_fpx2_8x1w_v2(r, a, b);")
    body.append("}")
    body.append("")
    body.append("void mul_fpx2_8x1w_v3(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b)")
    body.append("{")
    body.append("    mul_fpx2_8x1w_v2(r, a, b);")
    body.append("}")
    body.append("")
    body.append("void mul_fpx2_8x1w_v4(fpx2_8x1w r, const fp_8x1w a, const fp_8x1w b)")
    body.append("{")
    body.append("    mul_fpx2_8x1w_v2(r, a, b);")
    body.append("}")
    body.append("")

    footer = []
    footer.append("} // namespace blst::avx")
    footer.append("")
    footer.append("#endif")
    footer.append("")

    OUT_CPP.write_text("\n".join(header + body + footer), encoding="utf-8")


if __name__ == "__main__":
    main()
