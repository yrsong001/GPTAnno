from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple, Union, Dict
from pathlib import Path
import re
import difflib

import numpy as np
import pandas as pd


# -----------------------------
# Text normalization utilities
# -----------------------------

_PLURAL_RULES = [
    (r'(.*)ies$', r'\1y'),        # "bodies" -> "body"
    (r'(.*)ves$', r'\1f'),        # "leaves" -> "leaf" (approximate)
    (r'(.*)ses$', r'\1sis'),      # "bases"->"basis" (very approximate, heuristic)
    (r'(.*)s$', r'\1'),           # generic plural -> singular
]

def _to_singular(token: str) -> str:
    for pat, rep in _PLURAL_RULES:
        if re.match(pat, token):
            return re.sub(pat, rep, token)
    return token

def _normalize_text(s: str) -> str:
    # lower, trim, remove leading list numbering like "1.  T cell"
    s = s or ""
    s = s.strip().lower()
    s = re.sub(r'^\d+\.\s*', '', s)
    # collapse whitespace and dashes/underscores
    s = re.sub(r'[_\-]+', ' ', s)
    s = re.sub(r'\s+', ' ', s)
    # remove trivial suffix "cell(s)"
    s = s.replace(" cells", " cell")
    # singularize naively
    parts = [ _to_singular(p) for p in s.split() ]
    s = " ".join(parts)
    return s


# -----------------------------
# Mapping Data
# -----------------------------

@dataclass
class MappingDB:
    names: List[str]

    @classmethod
    def from_csv(cls, csv_path: Union[str, Path], name_col: Optional[str] = None) -> "MappingDB":
        df = pd.read_csv(csv_path)
        # Heuristic: if name_col is not provided, pick the first column that looks string-like.
        if name_col is None:
            for c in df.columns:
                if df[c].dtype == object:
                    name_col = c
                    break
        if name_col is None:
            raise ValueError("Could not infer the name column from the CSV. Please specify name_col.")
        names = [ _normalize_text(x) for x in df[name_col].astype(str).tolist() ]
        # de-duplicate while preserving order
        seen = set()
        uniq = []
        for n in names:
            if n not in seen:
                seen.add(n)
                uniq.append(n)
        return cls(uniq)


# -----------------------------
# Matching logic
# -----------------------------

@dataclass
class MatchResult:
    cleaned: str
    matched: Optional[str]
    method: Optional[str]
    score: float

def _best_fuzzy_match(query: str, candidates: List[str]) -> Tuple[Optional[str], float]:
    if not candidates:
        return None, 0.0
    matches = difflib.get_close_matches(query, candidates, n=1, cutoff=0.6)
    if not matches:
        return None, 0.0
    top = matches[0]
    # difflib doesn't give a score directly; we can approximate via SequenceMatcher
    score = difflib.SequenceMatcher(None, query, top).ratio()
    return top, float(score)

def clean_and_match_annotation(annotation: str,
                               mapping_db: MappingDB,
                               try_substring: bool = True,
                               try_fuzzy: bool = True) -> MatchResult:
    """
    Normalize an annotation string and map it to a controlled vocabulary.

    Strategy:
      1) exact match on normalized strings
      2) substring containment (e.g., "naive cd4 t" -> contains "t")
      3) fuzzy match via difflib

    Returns a MatchResult with the chosen method and a rough score.
    """
    cleaned = _normalize_text(annotation)
    # trivial unknowns
    if cleaned in {"unknown", "na", "n/a", "unk", ""}:
        return MatchResult(cleaned=cleaned, matched=None, method=None, score=0.0)

    # 1) exact
    if cleaned in mapping_db.names:
        return MatchResult(cleaned=cleaned, matched=cleaned, method="exact", score=1.0)

    # 2) substring (both ways) — heuristic
    if try_substring:
        for ref in mapping_db.names:
            if cleaned in ref or ref in cleaned:
                return MatchResult(cleaned=cleaned, matched=ref, method="substring", score=0.9)

    # 3) fuzzy
    if try_fuzzy:
        ref, sc = _best_fuzzy_match(cleaned, mapping_db.names)
        if ref is not None:
            return MatchResult(cleaned=cleaned, matched=ref, method="fuzzy", score=float(sc))

    # fallback
    return MatchResult(cleaned=cleaned, matched=None, method=None, score=0.0)


def map_annotations(series: Union[pd.Series, List[str]], mapping_db: MappingDB) -> pd.DataFrame:
    """
    Vectorize clean_and_match_annotation over a list/Series.
    Returns a DataFrame with columns: raw, cleaned, matched, method, score
    """
    if not isinstance(series, pd.Series):
        series = pd.Series(series, dtype="object")

    rows = []
    for raw in series.astype(str):
        res = clean_and_match_annotation(raw, mapping_db)
        rows.append({
            "raw": raw,
            "cleaned": res.cleaned,
            "matched": res.matched,
            "method": res.method,
            "score": res.score,
        })
    return pd.DataFrame(rows)



def load_gpt_mapping(csv_path: Union[str, Path], name_col_guess: Optional[str] = None) -> MappingDB:
    """
    Convenience wrapper that guesses the cell type name column from
    `/mnt/data/GPTCelltype_mapping.csv`-like files. If there's a column named
    'x' it will be used; otherwise we pick the first object dtype column.
    """
    df = pd.read_csv(csv_path)
    name_col = name_col_guess
    if name_col is None:
        if "x" in df.columns:
            name_col = "x"
        else:
            for c in df.columns:
                if df[c].dtype == object:
                    name_col = c
                    break
    if name_col is None:
        raise ValueError("Could not infer the mapping name column. Please specify name_col_guess.")
    # Reuse MappingDB to normalize and deduplicate
    return MappingDB.from_csv(csv_path, name_col=name_col)
