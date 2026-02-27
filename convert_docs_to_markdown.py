#!/usr/bin/env python3
"""
GROMACS RST → Obsidian-ready Markdown Converter
================================================

Converts all .rst files in docs/ to clean GitHub-Flavored Markdown in docs_markdown/,
then post-processes to fix Pandoc artifacts that break Obsidian rendering.

Requirements:
  - Python 3.8+
  - Pandoc installed and on PATH (https://pandoc.org/installing.html)

Usage:
  python convert_docs_to_markdown.py            # full clean conversion
  python convert_docs_to_markdown.py --fix-only  # re-run post-processing only (skip Pandoc)
"""

import argparse
import gc
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DOCS_SRC = Path("docs")
DOCS_OUT = Path("docs_markdown")

# Max retries when a directory is locked (e.g. by Obsidian on Windows)
RMTREE_RETRIES = 3
RMTREE_DELAY = 2  # seconds


def _rmtree_safe(path: Path):
    """Remove directory tree with retries for Windows file-lock issues."""
    for attempt in range(1, RMTREE_RETRIES + 1):
        try:
            gc.collect()  # release any Python-held file handles
            shutil.rmtree(path)
            return
        except PermissionError as e:
            if attempt == RMTREE_RETRIES:
                print(f"Error: cannot remove '{path}' — a program (e.g. Obsidian) may have it open.")
                print(f"  Close the program or manually delete '{path}', then re-run.")
                print(f"  ({e})")
                sys.exit(1)
            print(f"  Folder locked, retrying in {RMTREE_DELAY}s (attempt {attempt}/{RMTREE_RETRIES}) ...")
            time.sleep(RMTREE_DELAY)


# ===========================================================================
# Stage 1 — Pandoc conversion
# ===========================================================================

def run_pandoc_conversion():
    """Find every .rst in DOCS_SRC and convert to GFM markdown in DOCS_OUT.

    Uses incremental update: converts all .rst files in-place and then removes
    any .md files in DOCS_OUT that no longer have a corresponding .rst source.
    This avoids needing to delete the entire output directory (which fails if
    Obsidian or another app has it open).
    """
    if not DOCS_SRC.exists():
        print(f"Error: source directory '{DOCS_SRC}' not found.")
        sys.exit(1)

    if shutil.which("pandoc") is None:
        print("Error: pandoc is not installed or not on PATH.")
        print("Install it from https://pandoc.org/installing.html")
        sys.exit(1)

    rst_files = sorted(DOCS_SRC.rglob("*.rst"))
    if not rst_files:
        print(f"No .rst files found in {DOCS_SRC}")
        sys.exit(1)

    print(f"Stage 1: Converting {len(rst_files)} .rst files with Pandoc ...")

    # Build set of expected output paths (relative to DOCS_OUT)
    expected_md = set()

    failed = []
    for i, rst in enumerate(rst_files, 1):
        rel = rst.relative_to(DOCS_SRC)
        out = DOCS_OUT / rel.with_suffix(".md")
        out.parent.mkdir(parents=True, exist_ok=True)
        expected_md.add(out.resolve())

        result = subprocess.run(
            ["pandoc", str(rst), "-f", "rst", "-t", "gfm", "-o", str(out)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            failed.append((rst, result.stderr))
            print(f"  [{i}/{len(rst_files)}] FAILED: {rel}")
        elif i % 50 == 0 or i == len(rst_files):
            print(f"  [{i}/{len(rst_files)}] converted ...")

    # Remove stale .md files that no longer have a .rst source
    stale = 0
    if DOCS_OUT.exists():
        for md in DOCS_OUT.rglob("*.md"):
            if md.resolve() not in expected_md:
                try:
                    md.unlink()
                    stale += 1
                except OSError:
                    pass
        # Remove empty directories (bottom-up)
        for dirpath in sorted(DOCS_OUT.rglob("*"), reverse=True):
            if dirpath.is_dir():
                try:
                    dirpath.rmdir()  # only succeeds if empty
                except OSError:
                    pass
    if stale:
        print(f"  Removed {stale} stale .md files (no matching .rst source).")

    if failed:
        print(f"\n  WARNING: {len(failed)} files failed to convert:")
        for path, err in failed:
            print(f"    {path}: {err.strip()[:120]}")
    else:
        print(f"  All {len(rst_files)} files converted successfully.")

    return len(rst_files), len(failed)


# ===========================================================================
# Stage 2 — Post-processing (Obsidian cleanup)
# ===========================================================================

# --- Individual fixers (order matters) -------------------------------------

def fix_subst_links(text):
    r"""Replace Sphinx substitution links like [\|Gromacs\|](##SUBST##...) with plain text."""
    text = re.sub(
        r'\[\\?\|Gromacs\\?\|\]\(##SUBST##\|Gromacs\|\)', 'GROMACS', text
    )
    text = re.sub(
        r'\[\|Gromacs\|\]\(##SUBST##\|Gromacs\|\)', 'GROMACS', text
    )
    text = re.sub(
        r'\[\\?\|([^|\\]+)\\?\|\]\(##SUBST##\|[^)]+\|\)',
        r'`\1`', text
    )
    text = re.sub(
        r'\[\|([^|]+)\|\]\(##SUBST##\|[^)]+\|\)',
        lambda m: 'GROMACS' if m.group(1) == 'Gromacs' else f'`{m.group(1)}`',
        text,
    )
    return text


def fix_inline_math(text):
    """<span class="math inline">$...$</span> → $...$"""
    def _replace(m):
        c = re.sub(r'\s*\n\s*', ' ', m.group(1).strip())
        if c.startswith('$') and c.endswith('$'):
            c = c[1:-1]
        return f'${c}$'

    text = re.sub(
        r'<span\s*\n?class="math inline">(.*?)</span>', _replace,
        text, flags=re.DOTALL,
    )
    text = re.sub(
        r'<span\s+class="math inline">(.*?)</span>', _replace,
        text, flags=re.DOTALL,
    )
    return text


def fix_title_ref(text):
    """<span class="title-ref">X</span> → `X`"""
    return re.sub(
        r'<span class="title-ref">(.*?)</span>', r'`\1`',
        text, flags=re.DOTALL,
    )


def fix_span_labels(text):
    """Remove <span label="..."> wrappers around math blocks."""
    text = re.sub(r'<span label="[^"]*">\s*\n?', '', text)
    text = re.sub(r'^</span>\s*$', '', text, flags=re.MULTILINE)
    return text


def fix_gfm_inline_math(text):
    r"""Convert Pandoc GFM inline math $`...\`$ to Obsidian $...$ syntax."""
    return re.sub(r'\$`([^`]*)`\$', r'$\1$', text)


def fix_math_code_blocks(text):
    """Convert all ``` math blocks to $$ blocks for Obsidian.

    Handles three variants:
    1. Plain:       ``` math ... ```
    2. Blockquoted: > > ``` math ... > ```
    3. Indented:      ``` math ...   ```
    """
    # 1. Blockquoted math: lines prefixed with > (possibly nested >, with optional leading whitespace)
    def _replace_blockquote_math(m):
        content_block = m.group(1)
        # Strip leading whitespace + > markers from each line
        lines = []
        for line in content_block.split('\n'):
            cleaned = re.sub(r'^[\s>]+', '', line)
            lines.append(cleaned)
        return '$$\n' + '\n'.join(lines) + '\n$$'

    text = re.sub(
        r'^[ \t]*(?:>\s*)+``` ?math\s*\n((?:[ \t]*(?:>\s*).*\n)*?)[ \t]*(?:>\s*)+```\s*$',
        _replace_blockquote_math,
        text, flags=re.MULTILINE,
    )

    # 2. Indented math: lines prefixed with spaces
    def _replace_indented_math(m):
        indent = m.group(1)
        content_block = m.group(2)
        # Strip the same leading whitespace from content lines
        lines = []
        for line in content_block.split('\n'):
            # Remove up to len(indent) spaces from start
            if line.startswith(indent):
                lines.append(line[len(indent):])
            else:
                lines.append(line.lstrip())
        return '$$\n' + '\n'.join(lines) + '\n$$'

    text = re.sub(
        r'^([ \t]+)``` ?math\s*\n((?:.*\n)*?)\1```\s*$',
        _replace_indented_math,
        text, flags=re.MULTILINE,
    )

    # 3. Plain (no prefix) — original handler
    text = re.sub(
        r'``` ?math\s*\n(.*?)\n```', r'$$\n\1\n$$',
        text, flags=re.DOTALL,
    )

    return text


def fix_figure_blocks(text):
    """<figure>…</figure> → ![alt](src) + *caption*"""
    def _replace(m):
        block = m.group(0)
        img = re.search(r'<img[^>]*src="([^"]*)"[^>]*/?\s*>', block)
        alt = re.search(r'alt="([^"]*)"', block)
        cap = re.search(r'<figcaption>(.*?)</figcaption>', block, re.DOTALL)
        src = img.group(1) if img else ''
        alt_t = alt.group(1) if alt else ''
        cap_t = re.sub(r'<[^>]+>', '', cap.group(1)).strip() if cap else ''
        cap_t = re.sub(r'\s+', ' ', cap_t)
        if src:
            r = f'![{alt_t or cap_t}]({src})'
            return f'{r}\n\n*{cap_t}*' if cap_t else r
        return f'*{cap_t}*' if cap_t else ''

    return re.sub(r'<figure>\s*\n.*?</figure>', _replace, text, flags=re.DOTALL)


def fix_html_tables(text):
    """Convert <table>…</table> to Markdown pipe tables."""
    def _replace(m):
        block = m.group(0)
        cap_m = re.search(r'<caption>(.*?)</caption>', block, re.DOTALL)
        caption = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', cap_m.group(1)).strip()) if cap_m else ''

        headers = []
        thead = re.search(r'<thead>(.*?)</thead>', block, re.DOTALL)
        if thead:
            headers = [re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', h).strip())
                       for h in re.findall(r'<th[^>]*>(.*?)</th>', thead.group(1), re.DOTALL)]

        rows = []
        tbody = re.search(r'<tbody>(.*?)</tbody>', block, re.DOTALL)
        body = tbody.group(1) if tbody else block
        for rhtml in re.findall(r'<tr[^>]*>(.*?)</tr>', body, re.DOTALL):
            if '<th' in rhtml and thead:
                continue
            cells = [re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', c).strip())
                     for c in re.findall(r'<td[^>]*>(.*?)</td>', rhtml, re.DOTALL)]
            if any(cells):
                rows.append(cells)

        if not headers and not rows:
            return ''
        ncols = max(len(headers), max((len(r) for r in rows), default=0))
        if ncols == 0:
            return ''

        lines = []
        if caption:
            lines.append(f'**{caption}**\n')
        hdrs = headers if headers else [''] * ncols
        while len(hdrs) < ncols:
            hdrs.append('')
        lines.append('| ' + ' | '.join(hdrs) + ' |')
        lines.append('| ' + ' | '.join(['---'] * ncols) + ' |')
        for row in rows:
            while len(row) < ncols:
                row.append('')
            lines.append('| ' + ' | '.join(row) + ' |')
        return '\n'.join(lines)

    return re.sub(r'<table>.*?</table>', _replace, text, flags=re.DOTALL)


def fix_toctree(text):
    """Remove toctree divs."""
    return re.sub(r'<div class="toctree"[^>]*>\s*\n(.*?\n)*?</div>\s*\n?', '', text)


def fix_div_id_anchors(text):
    """Remove <div id="..."> anchor lines."""
    return re.sub(r'^<div id="[^"]*">\s*\n?', '', text, flags=re.MULTILINE)


def fix_sphinx_divs(text):
    """Convert Sphinx directive divs to Markdown equivalents."""
    # todo → callout
    text = re.sub(r'<div class="todo">\s*\n(.*?)\n</div>',
                  lambda m: f'> **TODO:** {m.group(1).strip()}\n', text, flags=re.DOTALL)
    text = re.sub(r'<div class="todolist">\s*\n?</div>\s*\n?', '', text)

    # seealso
    text = re.sub(r'<div class="seealso">\s*\n(.*?)\n</div>',
                  lambda m: f'> **See also:** {m.group(1).strip()}\n', text, flags=re.DOTALL)

    # admonitions
    text = re.sub(r'<div class="admonition\s+(\w+)">\s*\n(.*?)\n</div>',
                  lambda m: f'> **{m.group(1).strip().title()}:** {m.group(2).strip()}\n',
                  text, flags=re.DOTALL)

    # versionchanged
    text = re.sub(r'<div class="versionchanged">\s*\n(.*?)\n</div>',
                  lambda m: f'> *Changed:* {m.group(1).strip()}\n', text, flags=re.DOTALL)

    # cmake / function → code blocks
    for cls in ('cmake', 'function'):
        text = re.sub(rf'<div class="{cls}">\s*\n(.*?)\n</div>',
                      lambda m: f'```cmake\n{m.group(1).strip()}\n```\n', text, flags=re.DOTALL)

    # parsed-literal → code block
    text = re.sub(r'<div class="parsed-literal">\s*\n(.*?)\n</div>',
                  lambda m: f'```\n{m.group(1).strip()}\n```\n', text, flags=re.DOTALL)

    # Simple unwrap for these classes
    for cls in ('line-block', 'mdp-value', 'mdp', 'glossary', 'tab',
                'autofunction', 'autoclass', 'automodule', 'autoexception',
                'autodecorator', 'autodata', 'autoprogram', 'argparse', 'only'):
        text = re.sub(rf'<div class="{cls}"[^>]*>\s*\n(.*?)\n</div>',
                      r'\1\n', text, flags=re.DOTALL)

    # sidebar → blockquote
    text = re.sub(r'<div class="sidebar">\s*\n(.*?)\n</div>',
                  lambda m: '\n'.join(f'> {l}' for l in m.group(1).strip().split('\n')) + '\n',
                  text, flags=re.DOTALL)

    # digraph → dot code block
    text = re.sub(r'<div class="digraph">\s*\n(.*?)\n</div>',
                  lambda m: f'```dot\n{m.group(1).strip()}\n```\n', text, flags=re.DOTALL)

    # testsetup → remove
    text = re.sub(r'<div class="testsetup">\s*\n(.*?)\n</div>\s*\n?', '', text, flags=re.DOTALL)

    return text


def fix_interpreted_text_code(text):
    """<code class="interpreted-text" ...>X</code> → `X`"""
    text = re.sub(r'<code class="interpreted-text"[^>]*>(.*?)</code>',
                  r'`\1`', text, flags=re.DOTALL)
    text = re.sub(r'<code class="interpreted-text"\s*\n[^>]*>(.*?)</code>',
                  r'`\1`', text, flags=re.DOTALL)
    return text


def fix_remaining_html(text):
    """Clean up all residual HTML tags."""
    # Remove any remaining <div ...> and </div>
    text = re.sub(r'^\s*<div[^>]*>\s*$', '', text, flags=re.MULTILINE)
    text = re.sub(r'^\s*</div>\s*$', '', text, flags=re.MULTILINE)
    text = re.sub(r'</div>', '', text)

    # Stray table fragment tags
    for tag in ('table', 'thead', 'tbody', 'tr', 'td', 'th', 'caption', 'colgroup', 'col'):
        text = re.sub(rf'</?{tag}[^>]*>', '', text)

    # Inline <div class="line-block"> (content on same line)
    text = re.sub(r'<div class="line-block">\s*', '', text)

    # Stray spans
    text = re.sub(r'<span>(.*?)</span>', r'\1', text, flags=re.DOTALL)
    text = re.sub(r'<span id="[^"]*">', '', text)
    text = re.sub(r'</span>', '', text)

    # Inline HTML → Markdown equivalents
    text = re.sub(r'</?blockquote>', '', text)
    text = re.sub(r'</?p>', '', text)
    text = re.sub(r'<sup>(.*?)</sup>', r'^{\1}', text)
    text = re.sub(r'<sub>(.*?)</sub>', r'_{\1}', text)
    text = re.sub(r'<em>(.*?)</em>', r'*\1*', text)
    text = re.sub(r'<strong>(.*?)</strong>', r'**\1**', text)
    text = re.sub(r'<code>(.*?)</code>', r'`\1`', text)
    text = re.sub(r'<a\s+href="([^"]*)"[^>]*>(.*?)</a>', r'[\2](\1)', text)
    text = re.sub(r'<a[^>]*></a>', '', text)
    text = re.sub(r'</?li>', '', text)
    text = re.sub(r'<col[^>]*/>', '', text)
    text = re.sub(r'</?colgroup>', '', text)
    text = re.sub(r'</?gmx[^>]*>', '', text)
    text = re.sub(r'</?top[^>]*>', '', text)

    # Table cross-references
    text = re.sub(r'`Table %s <([^>]+)>`', 'the table below', text)

    # Orphaned blockquote lines (just > with optional surrounding whitespace)
    text = re.sub(r'^\s*>\s*$', '', text, flags=re.MULTILINE)

    # Excessive blank lines
    text = re.sub(r'\n{4,}', '\n\n\n', text)

    return text


# --- Pipeline --------------------------------------------------------------

FIXERS = [
    fix_subst_links,
    fix_gfm_inline_math,
    fix_inline_math,
    fix_title_ref,
    fix_span_labels,
    fix_math_code_blocks,
    fix_figure_blocks,
    fix_html_tables,
    fix_toctree,
    fix_div_id_anchors,
    fix_sphinx_divs,
    fix_interpreted_text_code,
    fix_remaining_html,
    fix_subst_links,  # run again after HTML tables have been flattened
]


def postprocess_all():
    """Apply all fixers to every .md file in DOCS_OUT."""
    md_files = sorted(DOCS_OUT.rglob("*.md"))
    print(f"Stage 2: Post-processing {len(md_files)} Markdown files ...")

    modified = 0
    for filepath in md_files:
        text = filepath.read_text(encoding="utf-8")
        original = text
        for fixer in FIXERS:
            text = fixer(text)
        if text != original:
            filepath.write_text(text, encoding="utf-8")
            modified += 1

    print(f"  {modified}/{len(md_files)} files cleaned up.")
    return modified


# ===========================================================================
# Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Convert GROMACS RST docs to Obsidian-ready Markdown.",
    )
    parser.add_argument(
        "--fix-only", action="store_true",
        help="Skip Pandoc conversion; only re-run the post-processing cleanup.",
    )
    args = parser.parse_args()

    os.chdir(Path(__file__).resolve().parent)

    if args.fix_only:
        if not DOCS_OUT.exists():
            print(f"Error: '{DOCS_OUT}' not found. Run without --fix-only first.")
            sys.exit(1)
        postprocess_all()
    else:
        total, failed = run_pandoc_conversion()
        postprocess_all()
        print(f"\nDone. {total - failed} files ready in {DOCS_OUT}/")


if __name__ == "__main__":
    main()
