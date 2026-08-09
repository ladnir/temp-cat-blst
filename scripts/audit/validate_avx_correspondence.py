#!/usr/bin/env python3
# Copyright 2026 Category Labs, Inc.
# SPDX-License-Identifier: Apache-2.0

"""Validate the source hashes and required fields in the AVX proof ledger."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


VALID_STATUSES = {"in_review", "accepted"}
HASH_FIELDS = (
    ("cpp_path", "cpp_sha256"),
    ("constant_path", "constant_sha256"),
)
ACCEPTED_REQUIRED_FIELDS = (
    "cpp_path",
    "cpp_commit",
    "cpp_sha256",
    "lean_declaration",
    "theorem",
    "portable_evidence",
    "native_evidence",
    "residual_assumptions",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def split_field(value: str) -> list[str]:
    return value.split("|") if value else []


def resolve_repo_path(repo: Path, value: str, row_id: str) -> Path:
    relative = Path(value)
    require(not relative.is_absolute(), f"{row_id}: absolute path is not allowed: {value}")
    resolved = (repo / relative).resolve()
    try:
        resolved.relative_to(repo)
    except ValueError as error:
        raise ValueError(f"{row_id}: path escapes the repository: {value}") from error
    require(resolved.is_file(), f"{row_id}: mapped file does not exist: {value}")
    return resolved


def canonical_source_bytes(path: Path, row_id: str) -> bytes:
    """Return the platform-independent byte stream recorded by the ledger.

    Every mapped correspondence artifact is text. Git stores these files with
    LF line endings, while a Windows worktree may materialize CRLF. Hashing the
    raw worktree bytes would make an accepted proof row platform-dependent.
    Reject legacy bare CR and normalize only the ordinary CRLF spelling.
    """
    contents = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in contents, f"{row_id}: mapped text contains bare CR: {path}")
    return contents


def validate_hashes(repo: Path, row: dict[str, str], path_field: str, hash_field: str) -> int:
    row_id = row["id"]
    paths = split_field(row[path_field])
    hashes = split_field(row[hash_field])
    require(
        len(paths) == len(hashes),
        f"{row_id}: {path_field} and {hash_field} have different item counts",
    )
    for relative, expected in zip(paths, hashes):
        require(
            len(expected) == 64 and all(character in "0123456789abcdef" for character in expected),
            f"{row_id}: invalid lowercase SHA-256 for {relative}",
        )
        path = resolve_repo_path(repo, relative, row_id)
        actual = hashlib.sha256(canonical_source_bytes(path, row_id)).hexdigest()
        require(actual == expected, f"{row_id}: stale hash for {relative}: {actual}")
    return len(paths)


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    ledger = repo / "proof" / "avx" / "CORRESPONDENCE.csv"
    with ledger.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames is not None, "correspondence ledger has no header")
        rows = list(reader)

    require(rows, "correspondence ledger is empty")
    seen: set[str] = set()
    mapped_files = 0
    for row in rows:
        row_id = row.get("id", "")
        require(row_id.startswith("AVX-P"), f"invalid correspondence id: {row_id!r}")
        require(row_id not in seen, f"duplicate correspondence id: {row_id}")
        seen.add(row_id)
        require(row["status"] in VALID_STATUSES, f"{row_id}: invalid status {row['status']!r}")

        for path_field, hash_field in HASH_FIELDS:
            mapped_files += validate_hashes(repo, row, path_field, hash_field)

        commits = split_field(row["cpp_commit"])
        sources = split_field(row["cpp_path"])
        require(
            len(commits) == len(sources),
            f"{row_id}: cpp_commit and cpp_path have different item counts",
        )
        for commit in commits:
            require(
                7 <= len(commit) <= 40 and all(character in "0123456789abcdef" for character in commit),
                f"{row_id}: malformed source commit {commit!r}",
            )

        if row["status"] == "accepted":
            for field in ACCEPTED_REQUIRED_FIELDS:
                require(row[field].strip(), f"{row_id}: accepted row has empty {field}")

    print(f"validated {len(rows)} AVX correspondence rows and {mapped_files} mapped file hashes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
