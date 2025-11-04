import os
import re
import json
import argparse
import hashlib
import time
from typing import List, Dict, Tuple, Any, Iterable, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

try:
    from pypdf import PdfReader  # preferred
except Exception:  # pragma: no cover
    PdfReader = None  # type: ignore

# Optional progress bar
try:
    from tqdm import tqdm  # type: ignore
except Exception:  # pragma: no cover
    def tqdm(iterable, **kwargs):  # type: ignore
        return iterable


def stable_dumps(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def compute_cache_key(model: str, messages: List[Dict[str, str]]) -> str:
    payload = {"model": model, "messages": messages}
    s = stable_dumps(payload)
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def cache_load(cache_dir: str, key: str) -> str:
    try:
        os.makedirs(cache_dir, exist_ok=True)
        path = os.path.join(cache_dir, f"{key}.json")
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                return f.read()
    except Exception:
        return ""
    return ""


def cache_save(cache_dir: str, key: str, content: str) -> None:
    try:
        os.makedirs(cache_dir, exist_ok=True)
        path = os.path.join(cache_dir, f"{key}.json")
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
    except Exception:
        pass


def read_pdf_pages(pdf_path: str) -> List[str]:
    if PdfReader is None:
        raise RuntimeError("pypdf is not installed. Please install 'pypdf'.")
    reader = PdfReader(pdf_path)
    pages: List[str] = []
    for page in reader.pages:
        try:
            text = page.extract_text() or ""
        except Exception:
            text = ""
        pages.append(text)
    return pages


def normalize_whitespace(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def split_into_sentences(text: str) -> List[str]:
    text = normalize_whitespace(text)
    # Simple sentence split; avoids heavy NLP deps
    parts = re.split(r"(?<=[\.!?])\s+(?=[A-Z\[])", text)
    sentences = [s.strip() for s in parts if s and len(s.strip()) > 0]
    return sentences


def sentence_is_candidate(text: str) -> bool:
    lower = text.lower()
    keywords = [
        "marker", "markers", "gene", "genes", "express", "expression",
        "highly", "signature", "characterized", "specific", "enriched", "upregulated",
    ]
    if any(k in lower for k in keywords):
        return True
    # comma-separated gene-like tokens (basic heuristic)
    if re.search(r"\b[A-Za-z0-9]{2,}(?:-[A-Za-z0-9]+)?(?:\s*,\s*[A-Za-z0-9-]{2,}){1,}\b", text):
        return True
    return False


def make_windows(sentences: List[str], max_len: int = 3, limit: Optional[int] = None) -> List[str]:
    windows: List[str] = []
    n = len(sentences)
    for i in range(n):
        for L in range(1, max_len + 1):
            if i + L <= n:
                chunk = normalize_whitespace(" ".join(sentences[i:i + L]))
                if sentence_is_candidate(chunk):
                    windows.append(chunk)
    # de-duplicate while preserving order
    seen: set = set()
    uniq: List[str] = []
    for w in windows:
        if w not in seen:
            seen.add(w)
            uniq.append(w)
    print(f'number of unique windows: {len(uniq)}')
    if limit is None:
        result = uniq
    else:
        result = uniq[:limit]
    print(f'number of windows to analyze: {len(result)}')
    return result


def get_model_params() -> Tuple[str, str, str]:
    api_user = os.getenv("API_USER_ID") or os.getenv("GPT_API_USER_ID") or os.getenv("GPT_API_USER") or ""
    api_password = (
        os.getenv("OPENAI_API_KEY")
        or os.getenv("API_PASSWORD")
        or os.getenv("GPT_API_PASSWORD")
        or ""
    )
    base_url = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1")
    return api_user, api_password, base_url


def call_chat_completion(messages: List[Dict[str, str]], model: str = "gpt-4o-mini") -> str:
    api_user, api_password, base_url = get_model_params()
    # Primary path: new SDK with graceful retry if temperature unsupported
    try:
        from openai import OpenAI  # type: ignore

        client = OpenAI(api_key=api_password, base_url=base_url)
        kwargs = {"model": model, "messages": messages}
        try:
            resp = client.chat.completions.create(temperature=0, **kwargs)
        except Exception as e:
            if "temperature" in str(e).lower():
                resp = client.chat.completions.create(**kwargs)
            else:
                raise
        return resp.choices[0].message.content  # type: ignore
    except Exception:
        # Fallback: raw requests, with retry if temperature unsupported
        import requests
        from requests.auth import HTTPBasicAuth

        url = base_url.rstrip("/") + "/chat/completions"
        headers = {"Content-Type": "application/json"}
        auth = None
        if api_user and api_password:
            auth = HTTPBasicAuth(api_user, api_password)
        else:
            headers["Authorization"] = f"Bearer {api_password}"

        def do_request(payload: Dict[str, Any]):
            r = requests.post(url, json=payload, headers=headers, auth=auth, timeout=120)
            r.raise_for_status()
            return r.json()

        payload: Dict[str, Any] = {"model": model, "messages": messages, "temperature": 0}
        try:
            data = do_request(payload)
        except requests.HTTPError as http_err:
            # Try to detect 'temperature' unsupported and retry without it
            msg = ""
            try:
                msg = http_err.response.json().get("error", {}).get("message", "")  # type: ignore
            except Exception:
                try:
                    msg = http_err.response.text  # type: ignore
                except Exception:
                    msg = str(http_err)
            if "temperature" in str(msg).lower():
                payload.pop("temperature", None)
                data = do_request(payload)
            else:
                raise
        return data["choices"][0]["message"]["content"]


def extract_citation_meta(first_page_text: str, model: str = "gpt-4o-mini") -> Dict[str, str]:
    prompt = (
        "You will receive text from the first page of a scientific PDF. "
        "Extract JSON with keys: first_author_surname, journal, year (YYYY). "
    )
    messages = [
        {"role": "system", "content": "You are a meticulous scientific assistant."},
        {"role": "user", "content": prompt + "\n\nTEXT:\n" + first_page_text[:8000]},
        {"role": "user", "content": "Respond ONLY with compact JSON."},
    ]
    content = call_chat_completion(messages, model=model)
    try:
        data = json.loads(content)
        return {
            "first_author_surname": str(data.get("first_author_surname", "Unknown")),
            "journal": str(data.get("journal", "Unknown Journal")),
            "year": str(data.get("year", "Unknown")),
        }
    except Exception:
        # Fallback heuristic for year
        m = re.search(r"(20\d{2}|19\d{2})", first_page_text)
        year = m.group(1) if m else "Unknown"
        return {"first_author_surname": "Unknown", "journal": "Unknown Journal", "year": year}


def build_extraction_prompt(window_text: str, name_map: Dict[str, str]) -> List[Dict[str, str]]:
    guidelines = (
        "Task: From the provided sentence(s), extract pairs of (cell type, marker genes) ONLY if the sentence explicitly or clearly implies the marker genes are highly expressed in that cell type. "
        "Ignore generic expression context or if the high expression relation is unclear. "
        "Use full cell type names consistently. If a short name appears and a full name is known from prior context, use the known full name. If not known, infer a reasonable full name."
        "Return JSON list under key 'items', each item: {full_cell_type_name, short_name, marker_genes}. "
        "marker_genes should be an array of official gene symbols (keep case as in the paper). "
        "Return an empty list if no valid pairs."
    )
    context = json.dumps({"known_name_map": name_map}, ensure_ascii=False)
    messages = [
        {"role": "system", "content": "You are an expert in single-cell biology, extracting precise relations."},
        {"role": "user", "content": guidelines},
        {"role": "user", "content": f"Known name map (short->full): {context}"},
        {"role": "user", "content": f"TEXT:\n{window_text}"},
        {"role": "user", "content": "Respond ONLY with JSON like {\"items\": [...]}"},
    ]
    return messages


def build_batch_extraction_prompt(indexed_windows: List[Tuple[int, str]], name_map: Dict[str, str]) -> List[Dict[str, str]]:
    guidelines = (
        "Task: For each labeled window below, extract pairs of (cell type, marker genes) ONLY if the window explicitly or clearly implies the marker genes are highly expressed in that cell type. "
        "Use full cell type names consistently; if a short name appears and a full name is known from prior context, use that full name. "
        "Return a single JSON object with key 'items' as a flat list. Each item must include: {window_index, full_cell_type_name, short_name, marker_genes}. "
        "marker_genes must be an array of official gene symbols. Return an empty list if no valid pairs across all windows."
    )
    context = json.dumps({"known_name_map": name_map}, ensure_ascii=False)
    header_lines = [f"[W{idx}] {txt}" for idx, txt in indexed_windows]
    windows_block = "\n\n" + "\n\n".join(header_lines)
    messages = [
        {"role": "system", "content": "You are an expert in single-cell biology, extracting precise relations."},
        {"role": "user", "content": guidelines},
        {"role": "user", "content": f"Known name map (short->full): {context}"},
        {"role": "user", "content": f"WINDOWS:{windows_block}"},
        {"role": "user", "content": "Respond ONLY with JSON like {\"items\": [{\"window_index\": <int>, ...}]}"},
    ]
    return messages


def generate_domain_heuristic_map() -> Dict[str, str]:
    mapping = {
        # Common cardiac/vascular abbreviations
        "cms": "cardiomyocytes",
        "cm": "cardiomyocyte",
        "fbs": "fibroblasts",
        "fb": "fibroblast",
        "ecs": "endothelial cells",
        "ec": "endothelial cell",
        "epcs": "epicardial cells",
        "epc": "epicardial cell",
        "epcs": "epicardial cells",
        "smcs": "smooth muscle cells",
        "smc": "smooth muscle cell",
        "vsmcs": "vascular smooth muscle cells",
        "cmcs": "cardiomyocytes",
        # Prefix chunks
        "endo": "endocardial",
        "epi": "epicardial",
        "vasc": "vascular",
        "atr": "atrial",
        "vent": "ventricular",
        "avc": "atrioventricular canal",
        "oft": "outflow tract",
    }
    return mapping


def build_pdf_abbreviation_map(text: str) -> Dict[str, str]:
    """Find patterns like 'Full Name (ABBR)' or 'ABBR (Full Name)'. Return short->full (lowercased key)."""
    short_to_full: Dict[str, str] = {}
    # Full (Short)
    for m in re.finditer(r"([A-Za-z][A-Za-z \-/]+?)\s*\(([^)]+)\)", text):
        full = m.group(1).strip()
        short = m.group(2).strip()
        if 1 <= len(short) <= 40 and len(full) >= 3:
            key = short.strip()
            short_to_full[key.lower()] = full
    # Short (Full)
    for m in re.finditer(r"([A-Za-z0-9_\-/]+)\s*\(([^)]+)\)", text):
        short = m.group(1).strip()
        full = m.group(2).strip()
        if 1 <= len(short) <= 40 and len(full) >= 3:
            key = short.strip()
            short_to_full[key.lower()] = full
    return short_to_full


def expand_short_cell_type(short_name: str, abbrev_map: Dict[str, str]) -> str:
    if not short_name:
        return ""
    # Direct mapping
    full = abbrev_map.get(short_name.lower())
    if full:
        return full
    # Split on '_' or '-' and expand chunks if possible
    domain_map = generate_domain_heuristic_map()
    parts = re.split(r"[_\-/]+", short_name)
    expanded_parts: List[str] = []
    for p in parts:
        if not p:
            continue
        key = p.lower()
        expanded = domain_map.get(key)
        if expanded:
            expanded_parts.append(expanded)
        else:
            # Try plural normalized key (e.g., 'ECs' -> 'ecs')
            if key.endswith("s") and len(key) > 2:
                maybe = domain_map.get(key[:-1])
                if maybe:
                    expanded_parts.append(maybe)
                    continue
            # Keep original token as fallback
            expanded_parts.append(p)
    if expanded_parts:
        return " ".join(expanded_parts)
    return short_name


def singularize_cell_type(name: str) -> str:
    original = name
    name = name.strip()
    # Do not alter explicit prefixes like pre-, pro-
    # Basic plural heuristics
    tokens = name.split()
    out_tokens: List[str] = []
    for tok in tokens:
        low = tok.lower()
        if low in {"pre", "pro"}:
            out_tokens.append(tok)
            continue
        # preserve '-like' containing tokens as-is at this step; decision to drop later
        if low.endswith("cells"):
            out_tokens.append(tok[:-1])  # cells -> cell
            continue
        if low.endswith("ies") and len(tok) > 4:
            out_tokens.append(tok[:-3] + "y")
            continue
        if low.endswith("sses"):
            out_tokens.append(tok[:-2])
            continue
        if low.endswith("s") and not low.endswith("ss") and len(tok) > 3 and not low.endswith("us"):
            out_tokens.append(tok[:-1])
        else:
            out_tokens.append(tok)
    result = " ".join(out_tokens)
    return result


def should_drop_cell_type(clean_name: str) -> bool:
    low = clean_name.lower()
    # Drop vague '-like' cell types entirely
    if "-like" in low:
        return True
    return False


def polish_cell_type(name: str, citation: Dict[str, str]) -> str:
    name = singularize_cell_type(name)
    # Replace trailing " <number>" pattern
    m = re.match(r"^(.*?\S)\s+(\d+)$", name)
    if m:
        base = m.group(1)
        idx = m.group(2)
        surname = citation.get("first_author_surname", "Unknown")
        journal = citation.get("journal", "Unknown Journal")
        year = citation.get("year", "Unknown")
        return f"{base} subtype {idx} ({surname} et al., {journal} {year})"
    return name


def clean_and_filter_items(items: List[Dict[str, Any]], citation: Dict[str, str], short_to_full_map: Dict[str, str]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    kept: List[Dict[str, Any]] = []
    dropped: List[Dict[str, Any]] = []
    for it in items:
        full = str(it.get("full_cell_type_name", "")).strip()
        short = str(it.get("short_name", "")).strip()
        # If full missing but short exists, try to expand using abbreviation maps
        if (not full) and short:
            guessed = expand_short_cell_type(short, short_to_full_map)
            if guessed:
                full = guessed
                it["full_cell_type_name"] = full
        markers = it.get("marker_genes", []) or []
        if isinstance(markers, str):
            # split on comma/semicolon
            markers = [m.strip() for m in re.split(r"[,;]", markers) if m.strip()]
        # must have markers
        if not full or len(markers) == 0:
            dropped.append(it)
            continue
        full_clean = polish_cell_type(full, citation)
        if should_drop_cell_type(full_clean):
            it["_reason"] = "contains -like"
            dropped.append(it)
            continue
        # rule 1: drop names starting with a single letter followed by digits (e.g., A1, C4, C4-0 ...)
        if re.match(r"^[A-Za-z]\d+", full_clean.strip()):
            it["_reason"] = "name starts with letter+digit prefix"
            dropped.append(it)
            continue
        # Dedup markers, keep order
        seen: set = set()
        uniq_markers: List[str] = []
        for m in markers:
            if not m:
                continue
            key = m.upper()
            if key not in seen:
                seen.add(key)
                uniq_markers.append(m)
        if len(uniq_markers) == 0:
            it["_reason"] = "no valid markers after cleanup"
            dropped.append(it)
            continue
        kept.append({
            "full_cell_type_name": full_clean,
            "short_name": short,
            "marker_genes": ", ".join(uniq_markers),
        })
    # rule 2: prefer '-positive-' over '-expressing-' when marker sets are identical
    def markers_key(s: str) -> Tuple[str, ...]:
        arr = [t.strip() for t in re.split(r"[,;]", s) if t.strip()]
        arr = [t.upper() for t in arr]
        return tuple(sorted(set(arr)))

    # Group kept by marker set
    groups: Dict[Tuple[str, ...], List[int]] = {}
    for idx, it in enumerate(kept):
        key = markers_key(it["marker_genes"])
        groups.setdefault(key, []).append(idx)

    to_remove: set = set()
    for key, indices in groups.items():
        # build maps of pattern A-(expressing|positive)-B
        expr_map: Dict[Tuple[str, str], List[int]] = {}
        pos_map: Dict[Tuple[str, str], List[int]] = {}
        for i in indices:
            name_l = kept[i]["full_cell_type_name"].lower().strip()
            m_expr = re.match(r"^(.*?)-expressing\s+(.*)$", name_l)
            m_pos = re.match(r"^(.*?)-positive\s+(.*)$", name_l)
            if m_expr:
                a = m_expr.group(1).strip()
                b = m_expr.group(2).strip()
                expr_map.setdefault((a, b), []).append(i)
            if m_pos:
                a = m_pos.group(1).strip()
                b = m_pos.group(2).strip()
                pos_map.setdefault((a, b), []).append(i)
        # if same (A,B) exists in pos and expr, drop all expr
        for ab in expr_map.keys():
            if ab in pos_map:
                for i in expr_map[ab]:
                    to_remove.add(i)
    if to_remove:
        new_kept: List[Dict[str, Any]] = []
        for idx, it in enumerate(kept):
            if idx in to_remove:
                it_copy = dict(it)
                it_copy["_reason"] = "prefer -positive- over -expressing- for identical marker set"
                dropped.append(it_copy)
            else:
                new_kept.append(it)
        kept = new_kept
    return kept, dropped


def chat_with_cache(messages: List[Dict[str, str]], model: str, cache_dir: str = None) -> str:
    if cache_dir:
        key = compute_cache_key(model, messages)
        cached = cache_load(cache_dir, key)
        if cached:
            return cached
        content = call_chat_completion(messages, model=model)
        cache_save(cache_dir, key, content)
        return content
    else:
        return call_chat_completion(messages, model=model)


def extract_from_windows_batched(
    windows: List[str],
    model: str = "gpt-4o-mini",
    initial_name_map: Dict[str, str] = None,
    batch_size: int = 8,
    max_workers: int = 4,
    cache_dir: str = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    all_items: List[Dict[str, Any]] = []
    short_to_full: Dict[str, str] = {}
    if initial_name_map:
        # Normalize keys to lower for robust lookup
        for k, v in initial_name_map.items():
            if not k:
                continue
            short_to_full[k.lower()] = v

    # Prepare batches of (1-based index, text)
    indexed = list(enumerate(windows, 1))
    batches: List[List[Tuple[int, str]]]= [indexed[i:i+batch_size] for i in range(0, len(indexed), batch_size)]

    def process_batch(batch: List[Tuple[int, str]]) -> List[Dict[str, Any]]:
        messages = build_batch_extraction_prompt(batch, short_to_full)
        content = chat_with_cache(messages, model=model, cache_dir=cache_dir)
        data = json.loads(content)
        items = data.get("items", [])
        # Attach provenance text from index
        idx_to_text = {idx: txt for idx, txt in batch}
        for it in items:
            try:
                widx = int(it.get("window_index"))
            except Exception:
                widx = None
            if widx in idx_to_text:
                it["_source_text"] = idx_to_text[widx]
        return items

    # Concurrent execution over batches
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(process_batch, b) for b in batches]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Extracting(batches)", unit="batch"):
            try:
                items = fut.result()
            except Exception:
                continue
            # update name map and collect
            for it in items:
                full = str(it.get("full_cell_type_name", "")).strip()
                short = str(it.get("short_name", "")).strip()
                if full and short and short.lower() not in short_to_full:
                    short_to_full[short.lower()] = full
            all_items.extend(items)
    return all_items, short_to_full


def choose_pdf_interactively(papers_dir: str) -> str:
    pdfs = [f for f in os.listdir(papers_dir) if f.lower().endswith(".pdf")]
    if not pdfs:
        raise FileNotFoundError(f"No PDF found in {papers_dir}")
    print("Available PDFs:")
    for i, f in enumerate(pdfs, 1):
        print(f"  [{i}] {f}")
    while True:
        sel = input("Select a PDF by number: ").strip()
        if not sel.isdigit():
            print("Please enter a valid number.")
            continue
        k = int(sel)
        if 1 <= k <= len(pdfs):
            return os.path.join(papers_dir, pdfs[k - 1])
        print("Out of range. Try again.")


def run(pdf_path: str, outputs_dir: str, model: str = "gpt-4o-mini", window_limit: Optional[int] = None) -> None:
    t0 = time.perf_counter()
    pages = read_pdf_pages(pdf_path)
    if not pages:
        raise RuntimeError("Empty PDF text.")
    first_page = pages[0] if pages else ""
    citation = extract_citation_meta(first_page, model=model)
    print(f'citation: {citation}')
    # Base name used for cache and output directories
    base = os.path.splitext(os.path.basename(pdf_path))[0]
    # Build abbreviation map from full text (seed for extraction + cleanup)
    full_text = "\n".join(pages)
    abbrev_map = build_pdf_abbreviation_map(full_text)
    # Merge with domain heuristics simple expansions (only for lookup, not overwrite PDF-derived)
    domain_map = generate_domain_heuristic_map()
    for k, v in domain_map.items():
        abbrev_map.setdefault(k.lower(), v)
    t_abbrev = time.perf_counter()
    # Gather sentences across pages
    all_sentences: List[str] = []
    for pg in pages:
        all_sentences.extend(split_into_sentences(pg))
    windows = make_windows(all_sentences, max_len=3, limit=window_limit)
    t_windows = time.perf_counter()

    # print(f"Windows to analyze: {len(windows)}")
    # Use batched, concurrent extraction with on-disk cache (under outputs/.cache per paper)
    cache_dir = os.path.join(outputs_dir, ".cache", base)
    items, name_map = extract_from_windows_batched(
        windows,
        model=model,
        initial_name_map=abbrev_map,
        batch_size=8,
        max_workers=4,
        cache_dir=cache_dir,
    )
    t_extract = time.perf_counter()
    print(f"Raw extracted items: {len(items)}")

    # Save 1-reference- with provenance before cleanup
    out_dir = os.path.join(outputs_dir, base)
    os.makedirs(out_dir, exist_ok=True)
    ref_csv = os.path.join(out_dir, f"1-reference-{base}.csv")
    ref_rows: List[Dict[str, Any]] = []
    for it in items:
        markers = it.get("marker_genes", []) or []
        if isinstance(markers, list):
            markers_str = ", ".join([m for m in markers if m])
        else:
            markers_str = str(markers)
        ref_rows.append({
            "full_cell_type_name": it.get("full_cell_type_name", ""),
            "short_name": it.get("short_name", ""),
            "marker_genes": markers_str,
            "source_text": it.get("_source_text", ""),
        })
    pd.DataFrame(ref_rows, columns=["full_cell_type_name", "short_name", "marker_genes", "source_text"]).to_csv(ref_csv, index=False)
    print(f"Saved: {ref_csv}")
    t_refsave = time.perf_counter()

    # Combine extraction-learned map and prebuilt abbrev map for cleanup fallback
    cleanup_map: Dict[str, str] = {}
    cleanup_map.update(abbrev_map)
    cleanup_map.update(name_map)
    kept, dropped = clean_and_filter_items(items, citation, cleanup_map)
    t_cleanup = time.perf_counter()
    if dropped:
        print("Dropped items (reason):")
        for it in dropped[:50]:  # limit print size
            reason = it.get("_reason", "missing markers or name")
            print(f" - {it.get('full_cell_type_name', it.get('short_name', 'UNKNOWN'))}: {reason}")
        if len(dropped) > 50:
            print(f" ... and {len(dropped)-50} more dropped entries")

    # Save 2-removed- with reasons and provenance
    removed_csv = os.path.join(out_dir, f"2-removed-{base}.csv")
    removed_rows: List[Dict[str, Any]] = []
    for it in dropped:
        markers = it.get("marker_genes", []) or []
        if isinstance(markers, list):
            markers_str = ", ".join([m for m in markers if m])
        else:
            markers_str = str(markers)
        removed_rows.append({
            "full_cell_type_name": it.get("full_cell_type_name", ""),
            "short_name": it.get("short_name", ""),
            "marker_genes": markers_str,
            "reason": it.get("_reason", ""),
            "source_text": it.get("_source_text", ""),
        })
    pd.DataFrame(removed_rows, columns=["full_cell_type_name", "short_name", "marker_genes", "reason", "source_text"]).to_csv(removed_csv, index=False)
    print(f"Saved: {removed_csv}")
    t_removedsave = time.perf_counter()

    # Deduplicate by full_cell_type_name
    dedup: Dict[str, Dict[str, Any]] = {}
    for it in kept:
        key = it["full_cell_type_name"].lower()
        if key not in dedup:
            dedup[key] = it
        else:
            # merge marker genes
            prev = dedup[key]["marker_genes"]
            merged = [m.strip() for m in re.split(r"[,;]", prev) if m.strip()]
            now = [m.strip() for m in re.split(r"[,;]", it["marker_genes"]) if m.strip()]
            seen: set = set(x.upper() for x in merged)
            for g in now:
                if g.upper() not in seen:
                    merged.append(g)
                    seen.add(g.upper())
            dedup[key]["marker_genes"] = ", ".join(merged)

    out_rows = list(dedup.values())
    df = pd.DataFrame(out_rows, columns=["full_cell_type_name", "short_name", "marker_genes"])
    final_csv = os.path.join(out_dir, f"3-final-{base}.csv")
    df.to_csv(final_csv, index=False)
    print(f"Saved: {final_csv}")
    t_finalsave = time.perf_counter()

    # Save runtime statistics next to CSVs
    stats = {
        "pdf": base,
        "model": model,
        "num_pages": len(pages),
        "num_sentences": len(all_sentences),
        "num_windows": len(windows),
        "raw_items": len(items),
        "kept_items": len(kept),
        "dropped_items": len(dropped),
        "final_rows": len(out_rows),
        "batch_size": 8,
        "max_workers": 4,
        "cache_dir": cache_dir,
        "timings_seconds": {
            "build_abbrev": round(t_abbrev - t0, 3),
            "split_and_windows": round(t_windows - t_abbrev, 3),
            "extraction": round(t_extract - t_windows, 3),
            "save_reference": round(t_refsave - t_extract, 3),
            "cleanup": round(t_cleanup - t_refsave, 3),
            "save_removed": round(t_removedsave - t_cleanup, 3),
            "dedup_and_save_final": round(t_finalsave - t_removedsave, 3),
            "total": round(t_finalsave - t0, 3),
        },
    }
    runtime_json = os.path.join(out_dir, "runtime.json")
    with open(runtime_json, "w", encoding="utf-8") as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    print(f"Saved: {runtime_json}")


def main():
    parser = argparse.ArgumentParser(description="Extract cell types and marker genes from a PDF using GPT.")
    parser.add_argument("--pdf", type=str, default="", help="Path to input PDF. If omitted, choose interactively from ./papers.")
    parser.add_argument("--out", type=str, default="", help="Output directory. Default: ./outputs next to script.")
    parser.add_argument("--model", type=str, default=os.getenv("OPENAI_MODEL", "gpt-5-nano"), help="Model name.") # gpt-4o-mini
    parser.add_argument("--papers_dir", type=str, default="papers", help="Directory to look for PDFs if --pdf not provided.")
    parser.add_argument("--limit", type=int, default=None, help="Max number of windows to analyze; omit to use all.")
    args = parser.parse_args()

    if not args.pdf:
        # Resolve papers_dir relative to script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        papers_dir = args.papers_dir
        if not os.path.isabs(papers_dir):
            papers_dir = os.path.join(script_dir, papers_dir)
        pdf_path = choose_pdf_interactively(papers_dir)
    else:
        pdf_path = args.pdf
        if not os.path.isabs(pdf_path):
            pdf_path = os.path.abspath(pdf_path)
    if not os.path.exists(pdf_path):
        raise FileNotFoundError(pdf_path)

    # Resolve outputs directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.out:
        out_path = args.out
        if not os.path.isabs(out_path):
            out_path = os.path.abspath(out_path)
        if out_path.lower().endswith(".csv"):
            outputs_dir = os.path.dirname(out_path) or os.path.join(script_dir, "outputs")
        else:
            outputs_dir = out_path
    else:
        outputs_dir = os.path.join(script_dir, "outputs")

    api_user, api_password, _ = get_model_params()
    if not api_password:
        raise RuntimeError("API credentials not set. Please export API_PASSWORD or OPENAI_API_KEY.")

    run(pdf_path, outputs_dir, model=args.model, window_limit=args.limit)


if __name__ == "__main__":
    main()