#!/usr/bin/env python3
"""Build and analyze dependency graphs for redc-core scheduling work.

Two modes are supported:
1) Template mode: emit a semantic DAG for a modeled redc_core schedule.
2) Asm mode: parse GNU objdump text and emit a def-use DAG from instructions.

This tool is offline analysis only. It does not execute any benchmark code.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


@dataclass
class Node:
    id: int
    order: int
    op: str
    label: str
    latency: int
    defs: List[str] = field(default_factory=list)
    uses: List[str] = field(default_factory=list)
    addr: Optional[str] = None


@dataclass
class Edge:
    src: int
    dst: int
    resource: str


class Graph:
    def __init__(self) -> None:
        self.nodes: List[Node] = []
        self.edges: List[Edge] = []
        self._pred: Dict[int, List[int]] = defaultdict(list)
        self._succ: Dict[int, List[int]] = defaultdict(list)

    def add_node(
        self,
        op: str,
        label: str,
        latency: int,
        defs: Sequence[str],
        uses: Sequence[str],
        addr: Optional[str] = None,
    ) -> int:
        node_id = len(self.nodes)
        self.nodes.append(
            Node(
                id=node_id,
                order=node_id,
                op=op,
                label=label,
                latency=latency,
                defs=list(defs),
                uses=list(uses),
                addr=addr,
            )
        )
        return node_id

    def add_edge(self, src: int, dst: int, resource: str) -> None:
        self.edges.append(Edge(src=src, dst=dst, resource=resource))
        self._pred[dst].append(src)
        self._succ[src].append(dst)

    def preds(self, node_id: int) -> List[int]:
        return self._pred.get(node_id, [])


def add_node_with_dataflow(
    graph: Graph,
    last_def: Dict[str, int],
    *,
    op: str,
    label: str,
    latency: int,
    defs: Sequence[str],
    uses: Sequence[str],
    addr: Optional[str] = None,
) -> int:
    node_id = graph.add_node(op=op, label=label, latency=latency, defs=defs, uses=uses, addr=addr)
    for resource in uses:
        pred = last_def.get(resource)
        if pred is not None:
            graph.add_edge(pred, node_id, resource)
    for resource in defs:
        last_def[resource] = node_id
    return node_id


def latency_for_op(op: str) -> int:
    table = {
        "init": 0,
        "u_gen": 5,
        "ifma_lo": 4,
        "ifma_hi": 4,
        "shift_add": 1,
        "tail_add": 1,
        "carry_add": 1,
        "mask": 1,
        "vpmadd52luq": 4,
        "vpmadd52huq": 4,
        "vpaddq": 1,
        "vpsllq": 1,
        "vpsraq": 1,
        "vpandd": 1,
        "vpandq": 1,
        "vmovdqa64": 1,
        "vmovdqu64": 1,
        "vmovdqa32": 1,
        "mov": 1,
        "lea": 1,
        "add": 1,
        "sub": 1,
        "xor": 1,
        "and": 1,
        "cmp": 1,
    }
    return table.get(op, 1)


def build_template_graph(streams: int, schedule: str) -> Graph:
    if streams <= 0:
        raise ValueError("streams must be > 0")
    if schedule not in {"round-robin", "stream-major"}:
        raise ValueError("schedule must be 'round-robin' or 'stream-major'")

    graph = Graph()
    last_def: Dict[str, int] = {}

    # Initial state: a[0..15] input, y[0..14] = 0
    for s in range(streams):
        for i in range(16):
            add_node_with_dataflow(
                graph,
                last_def,
                op="init",
                label=f"init_a_s{s}_i{i}",
                latency=latency_for_op("init"),
                defs=[f"a{s}_{i}"],
                uses=[],
            )
        for i in range(15):
            add_node_with_dataflow(
                graph,
                last_def,
                op="init",
                label=f"init_y_s{s}_i{i}",
                latency=latency_for_op("init"),
                defs=[f"y{s}_{i}"],
                uses=[],
            )

    if schedule == "round-robin":
        round_stream_order = [(r, s) for r in range(8) for s in range(streams)]
    else:
        round_stream_order = [(r, s) for s in range(streams) for r in range(8)]

    for r, s in round_stream_order:
        u = f"u{s}_{r}"
        add_node_with_dataflow(
            graph,
            last_def,
            op="u_gen",
            label=f"u_gen_s{s}_r{r}",
            latency=latency_for_op("u_gen"),
            defs=[u],
            uses=[f"a{s}_{r}"],
        )

        for j in range(8):
            a_name = f"a{s}_{r + j}"
            y_name = f"y{s}_{r + j}"
            add_node_with_dataflow(
                graph,
                last_def,
                op="ifma_lo",
                label=f"ifma_lo_s{s}_r{r}_j{j}",
                latency=latency_for_op("ifma_lo"),
                defs=[a_name],
                uses=[a_name, u],
            )
            add_node_with_dataflow(
                graph,
                last_def,
                op="ifma_hi",
                label=f"ifma_hi_s{s}_r{r}_j{j}",
                latency=latency_for_op("ifma_hi"),
                defs=[y_name],
                uses=[y_name, u],
            )

        add_node_with_dataflow(
            graph,
            last_def,
            op="shift_add",
            label=f"carry_fold_s{s}_r{r}",
            latency=latency_for_op("shift_add"),
            defs=[f"a{s}_{r + 1}"],
            uses=[f"a{s}_{r + 1}", f"a{s}_{r}", f"y{s}_{r}"],
        )

    # tail: a9..a15 += shl(y8..y14)
    for s in range(streams):
        for i in range(9, 16):
            add_node_with_dataflow(
                graph,
                last_def,
                op="tail_add",
                label=f"tail_add_s{s}_a{i}",
                latency=latency_for_op("tail_add"),
                defs=[f"a{s}_{i}"],
                uses=[f"a{s}_{i}", f"y{s}_{i - 1}"],
            )

    # carry propagation on a8..a15 (interleaved across streams to model k-wise ILP)
    for step in range(7):
        for s in range(streams):
            add_node_with_dataflow(
                graph,
                last_def,
                op="carry_add",
                label=f"carry_add_s{s}_step{step}",
                latency=latency_for_op("carry_add"),
                defs=[f"a{s}_{9 + step}"],
                uses=[f"a{s}_{9 + step}", f"a{s}_{8 + step}"],
            )
            add_node_with_dataflow(
                graph,
                last_def,
                op="mask",
                label=f"mask_s{s}_idx{8 + step}",
                latency=latency_for_op("mask"),
                defs=[f"a{s}_{8 + step}"],
                uses=[f"a{s}_{8 + step}"],
            )

    for s in range(streams):
        add_node_with_dataflow(
            graph,
            last_def,
            op="mask",
            label=f"mask_s{s}_idx15",
            latency=latency_for_op("mask"),
            defs=[f"a{s}_15"],
            uses=[f"a{s}_15"],
        )

    return graph


def split_operands(text: str) -> List[str]:
    out: List[str] = []
    cur: List[str] = []
    depth = 0
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        elif ch == "," and depth == 0:
            op = "".join(cur).strip()
            if op:
                out.append(op)
            cur = []
            continue
        cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        out.append(tail)
    return out


GPR_ALIAS = {
    "eax": "rax",
    "ax": "rax",
    "al": "rax",
    "ebx": "rbx",
    "bx": "rbx",
    "bl": "rbx",
    "ecx": "rcx",
    "cx": "rcx",
    "cl": "rcx",
    "edx": "rdx",
    "dx": "rdx",
    "dl": "rdx",
    "esi": "rsi",
    "si": "rsi",
    "sil": "rsi",
    "edi": "rdi",
    "di": "rdi",
    "dil": "rdi",
    "ebp": "rbp",
    "bp": "rbp",
    "bpl": "rbp",
    "esp": "rsp",
    "sp": "rsp",
    "spl": "rsp",
}


def normalize_reg(reg: str) -> str:
    r = reg.lstrip("%").lower()
    return GPR_ALIAS.get(r, r)


def normalize_mem(op: str) -> Optional[str]:
    if "(" not in op or ")" not in op:
        return None
    # Example: -0x140(%rsp), (%rbp), 0x20(%rax,%rdx,8)
    m = re.search(
        r"(?P<disp>[-+]?0x[0-9a-fA-F]+|[-+]?\d+)?\((?P<base>%[a-z0-9]+)(?:,(?P<index>%[a-z0-9]+)(?:,(?P<scale>\d+))?)?\)",
        op.replace(" ", ""),
    )
    if not m:
        return f"mem:{op}"
    disp = m.group("disp") or "0"
    base = normalize_reg(m.group("base"))
    index = m.group("index")
    scale = m.group("scale")
    idx_s = ""
    if index:
        idx_s = f"+{normalize_reg(index)}*{scale or '1'}"
    return f"mem[{base}{disp}{idx_s}]"


def operand_resources(op: str) -> Tuple[List[str], Optional[str]]:
    regs = [normalize_reg(x) for x in re.findall(r"%[a-zA-Z0-9]+", op)]
    mem = normalize_mem(op)
    return regs, mem


def classify_asm_instruction(mnemonic: str, operands: List[str]) -> Tuple[List[str], List[str]]:
    m = mnemonic.lower()
    op_res = [operand_resources(op) for op in operands]

    def regs(i: int) -> List[str]:
        return op_res[i][0] if 0 <= i < len(op_res) else []

    def mem(i: int) -> Optional[str]:
        return op_res[i][1] if 0 <= i < len(op_res) else None

    defs: List[str] = []
    uses: List[str] = []

    if m.startswith("vpmadd52") and len(operands) >= 3:
        dst = regs(2)
        defs += dst
        uses += dst + regs(0) + regs(1)
        if mem(0):
            uses.append(mem(0))
        if mem(1):
            uses.append(mem(1))
        return defs, uses

    if m.startswith("vmov") and len(operands) == 2:
        src_mem = mem(0)
        dst_mem = mem(1)
        src_regs = regs(0)
        dst_regs = regs(1)
        if src_mem and dst_regs:
            defs += dst_regs
            uses.append(src_mem)
            uses += [r for r in src_regs if r not in dst_regs]
            return defs, uses
        if src_regs and dst_mem:
            defs.append(dst_mem)
            uses += src_regs
            return defs, uses
        if src_regs and dst_regs:
            defs += dst_regs
            uses += src_regs
            return defs, uses

    if m in {"mov", "movq", "movl"} and len(operands) == 2:
        src_mem = mem(0)
        dst_mem = mem(1)
        src_regs = regs(0)
        dst_regs = regs(1)
        if dst_mem:
            defs.append(dst_mem)
            uses += src_regs
            if src_mem:
                uses.append(src_mem)
            return defs, uses
        defs += dst_regs
        uses += src_regs
        if src_mem:
            uses.append(src_mem)
        return defs, uses

    if m == "lea" and len(operands) == 2:
        defs += regs(1)
        uses += regs(0)
        if mem(0):
            uses.append(mem(0))
        return defs, uses

    if m in {"add", "sub", "and", "or", "xor"} and len(operands) == 2:
        # AT&T two-operand form: src, dst (dst is read-modify-write)
        defs += regs(1)
        if mem(1):
            defs.append(mem(1))
        uses += regs(0) + regs(1)
        if mem(0):
            uses.append(mem(0))
        if mem(1):
            uses.append(mem(1))
        return defs, uses

    if m in {"cmp", "test"} and len(operands) == 2:
        uses += regs(0) + regs(1)
        if mem(0):
            uses.append(mem(0))
        if mem(1):
            uses.append(mem(1))
        return defs, uses

    if m.startswith("vp") and len(operands) >= 3:
        defs += regs(2)
        uses += regs(0) + regs(1)
        if mem(0):
            uses.append(mem(0))
        if mem(1):
            uses.append(mem(1))
        return defs, uses

    if len(operands) >= 2 and regs(len(operands) - 1):
        defs += regs(len(operands) - 1)
        for i in range(len(operands) - 1):
            uses += regs(i)
            if mem(i):
                uses.append(mem(i))
        return defs, uses

    for i in range(len(operands)):
        uses += regs(i)
        if mem(i):
            uses.append(mem(i))
    return defs, uses


def parse_objdump_text(
    text: str, start_addr: Optional[int], end_addr: Optional[int]
) -> List[Tuple[str, str, List[str]]]:
    out: List[Tuple[str, str, List[str]]] = []
    for line in text.splitlines():
        if ":" not in line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        addr_part = parts[0].strip()
        if not addr_part.endswith(":"):
            continue
        addr_hex = addr_part[:-1]
        if not re.fullmatch(r"[0-9a-fA-F]+", addr_hex):
            continue
        addr = int(addr_hex, 16)
        if start_addr is not None and addr < start_addr:
            continue
        if end_addr is not None and addr >= end_addr:
            continue
        tail = parts[-1].strip()
        if not tail:
            continue
        tokens = tail.split(None, 1)
        mnemonic = tokens[0].lower()
        if not re.match(r"^[a-z.][a-z0-9.]*$", mnemonic):
            continue
        operands = split_operands(tokens[1]) if len(tokens) > 1 else []
        out.append((addr_hex.lower(), mnemonic, operands))
    return out


def build_asm_graph(text: str, start_addr: Optional[int], end_addr: Optional[int]) -> Graph:
    graph = Graph()
    last_def: Dict[str, int] = {}
    for addr_hex, mnemonic, operands in parse_objdump_text(text, start_addr, end_addr):
        defs, uses = classify_asm_instruction(mnemonic, operands)
        add_node_with_dataflow(
            graph,
            last_def,
            op=mnemonic,
            label=f"{addr_hex}: {mnemonic} {' '.join(operands)}".strip(),
            latency=latency_for_op(mnemonic),
            defs=defs,
            uses=uses,
            addr=addr_hex,
        )
    return graph


def critical_path(graph: Graph) -> Tuple[int, List[int]]:
    if not graph.nodes:
        return 0, []
    cp = [0] * len(graph.nodes)
    for node in graph.nodes:
        preds = graph.preds(node.id)
        best_pred = max((cp[p] for p in preds), default=0)
        cp[node.id] = best_pred + node.latency
    return max(cp), cp


def live_pressure(graph: Graph) -> Tuple[int, int]:
    # Returns (max_live_vector_regs_like, max_live_stack_slots)
    use_pos: Dict[str, List[int]] = defaultdict(list)
    for node in graph.nodes:
        for u in node.uses:
            use_pos[u].append(node.order)

    ptr: Dict[str, int] = defaultdict(int)
    live: set[str] = set()
    max_live_vec = 0
    max_live_stack = 0

    def is_vec(res: str) -> bool:
        if res.startswith(("zmm", "ymm", "xmm")):
            return True
        if res.startswith(("a", "y", "u")) and "_" in res:
            return True
        return False

    def is_stack(res: str) -> bool:
        return res.startswith("mem[rsp") or res.startswith("mem[rbp")

    for node in graph.nodes:
        # Consume uses at this node.
        used_set = set(node.uses)
        for u in node.uses:
            positions = use_pos.get(u, [])
            p = ptr[u]
            if p < len(positions) and positions[p] == node.order:
                ptr[u] = p + 1

        # Drop used resources with no further uses.
        for u in used_set:
            if ptr[u] >= len(use_pos.get(u, [])):
                live.discard(u)

        # Process defs: kill old value and add new one only if used later.
        for d in node.defs:
            live.discard(d)
            positions = use_pos.get(d, [])
            p = ptr[d]
            if p < len(positions) and positions[p] > node.order:
                live.add(d)

        vec_count = sum(1 for r in live if is_vec(r))
        stack_count = sum(1 for r in live if is_stack(r))
        if vec_count > max_live_vec:
            max_live_vec = vec_count
        if stack_count > max_live_stack:
            max_live_stack = stack_count

    return max_live_vec, max_live_stack


def graph_to_json(graph: Graph) -> Dict[str, object]:
    return {
        "nodes": [
            {
                "id": n.id,
                "order": n.order,
                "addr": n.addr,
                "op": n.op,
                "label": n.label,
                "latency": n.latency,
                "defs": n.defs,
                "uses": n.uses,
            }
            for n in graph.nodes
        ],
        "edges": [{"src": e.src, "dst": e.dst, "resource": e.resource} for e in graph.edges],
    }


def write_dot(graph: Graph, path: Path, max_nodes: int) -> None:
    keep = set(range(min(max_nodes, len(graph.nodes))))
    with path.open("w", encoding="utf-8") as f:
        f.write("digraph redc_dag {\n")
        f.write("  rankdir=LR;\n")
        for n in graph.nodes:
            if n.id not in keep:
                continue
            label = n.label.replace('"', "'")
            f.write(f'  n{n.id} [label="{n.id}: {label}"];\n')
        for e in graph.edges:
            if e.src in keep and e.dst in keep:
                res = e.resource.replace('"', "'")
                f.write(f'  n{e.src} -> n{e.dst} [label="{res}"];\n')
        f.write("}\n")


def print_summary(graph: Graph, title: str) -> None:
    cp, _ = critical_path(graph)
    max_live_vec, max_live_stack = live_pressure(graph)
    ops = Counter(n.op for n in graph.nodes)
    print(f"{title}")
    print(f"  nodes: {len(graph.nodes)}")
    print(f"  edges: {len(graph.edges)}")
    print(f"  critical_path_latency: {cp}")
    print(f"  max_live_vec_like: {max_live_vec}")
    print(f"  max_live_stack_slots: {max_live_stack}")
    print("  top_ops:")
    for op, cnt in ops.most_common(12):
        print(f"    {op}: {cnt}")


def parse_hex_addr(v: Optional[str]) -> Optional[int]:
    if v is None:
        return None
    x = v.strip().lower()
    if x.startswith("0x"):
        x = x[2:]
    return int(x, 16)


def run_template(args: argparse.Namespace) -> int:
    graph = build_template_graph(streams=args.streams, schedule=args.schedule)
    print_summary(graph, "Template DAG")
    payload = graph_to_json(graph)
    payload["meta"] = {
        "mode": "template",
        "streams": args.streams,
        "schedule": args.schedule,
    }
    if args.output_json:
        out = Path(args.output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"  wrote_json: {out}")
    if args.output_dot:
        out_dot = Path(args.output_dot)
        out_dot.parent.mkdir(parents=True, exist_ok=True)
        write_dot(graph, out_dot, max_nodes=args.dot_max_nodes)
        print(f"  wrote_dot: {out_dot}")
    return 0


def run_asm(args: argparse.Namespace) -> int:
    inp = Path(args.input)
    text = inp.read_text(encoding="utf-8", errors="ignore")
    graph = build_asm_graph(
        text=text,
        start_addr=parse_hex_addr(args.start_addr),
        end_addr=parse_hex_addr(args.end_addr),
    )
    print_summary(graph, "Asm DAG")
    payload = graph_to_json(graph)
    payload["meta"] = {
        "mode": "asm",
        "input": str(inp),
        "start_addr": args.start_addr,
        "end_addr": args.end_addr,
    }
    if args.output_json:
        out = Path(args.output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"  wrote_json: {out}")
    if args.output_dot:
        out_dot = Path(args.output_dot)
        out_dot.parent.mkdir(parents=True, exist_ok=True)
        write_dot(graph, out_dot, max_nodes=args.dot_max_nodes)
        print(f"  wrote_dot: {out_dot}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build redc scheduling dependency graphs.")
    sp = p.add_subparsers(dest="mode", required=True)

    pt = sp.add_parser("template", help="Build semantic DAG from redc template model.")
    pt.add_argument("--streams", type=int, default=2, help="Number of independent streams (default: 2).")
    pt.add_argument(
        "--schedule",
        choices=["round-robin", "stream-major"],
        default="round-robin",
        help="Round/stream ordering in template generation.",
    )
    pt.add_argument("--output-json", type=str, default="", help="Write graph JSON to this path.")
    pt.add_argument("--output-dot", type=str, default="", help="Write Graphviz DOT to this path.")
    pt.add_argument(
        "--dot-max-nodes",
        type=int,
        default=400,
        help="Limit number of nodes in DOT output (default: 400).",
    )
    pt.set_defaults(func=run_template)

    pa = sp.add_parser("asm", help="Build DAG by parsing objdump assembly text.")
    pa.add_argument("--input", type=str, required=True, help="Path to objdump text file.")
    pa.add_argument("--start-addr", type=str, default="", help="Optional start address (hex).")
    pa.add_argument("--end-addr", type=str, default="", help="Optional end address (hex, exclusive).")
    pa.add_argument("--output-json", type=str, default="", help="Write graph JSON to this path.")
    pa.add_argument("--output-dot", type=str, default="", help="Write Graphviz DOT to this path.")
    pa.add_argument(
        "--dot-max-nodes",
        type=int,
        default=400,
        help="Limit number of nodes in DOT output (default: 400).",
    )
    pa.set_defaults(func=run_asm)

    return p


def main(argv: Sequence[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.mode == "asm":
        if args.start_addr == "":
            args.start_addr = None
        if args.end_addr == "":
            args.end_addr = None
    if args.mode in {"template", "asm"}:
        if not getattr(args, "output_json", ""):
            args.output_json = None
        if not getattr(args, "output_dot", ""):
            args.output_dot = None
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
