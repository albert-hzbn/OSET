# -*- coding: utf-8 -*-
"""
Convert OSET_manuscript.md to OSET_manuscript.tex  (pdflatex-compatible).
"""
import re

# ---------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[12pt,a4paper]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,amsthm}
\usepackage{bm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{array}
\usepackage{longtable}
\usepackage{geometry}
\usepackage{setspace}
\usepackage[hidelinks,bookmarksnumbered]{hyperref}
\usepackage{caption}
\usepackage{enumerate}
\usepackage{xcolor}
\usepackage{listings}
\usepackage{mathtools}
\usepackage{float}
\usepackage{cancel}
\usepackage{microtype}

\geometry{margin=2.5cm}
\doublespacing

\newtheorem{theorem}{Theorem}[section]
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{lemma}[theorem]{Lemma}
\theoremstyle{remark}
\newtheorem{remark}{Remark}

\lstset{basicstyle=\small\ttfamily, frame=single,
        breaklines=true, keepspaces=true, columns=flexible}

\newcommand{\PV}{\mathrm{P.V.}}
\newcommand{\beff}{b_{\mathrm{eff}}}
\newcommand{\DeltaF}{\Delta F_0}

\begin{document}

"""

POSTAMBLE = "\n\\end{document}\n"

# ---------------------------------------------------------------------------
# Unicode replacements  (use explicit code points to avoid source-encoding
# issues on any platform)
# ---------------------------------------------------------------------------
UNICODE_TEXT = [
    # Arrows
    ('↔', r'$\leftrightarrow$'),   # ↔
    ('→', r'$\rightarrow$'),       # →
    ('←', r'$\leftarrow$'),        # ←
    ('⇒', r'$\Rightarrow$'),       # ⇒
    ('⇐', r'$\Leftarrow$'),        # ⇐
    ('⇔', r'$\Leftrightarrow$'),   # ⇔
    ('≠', r'$\neq$'),              # ≠
    ('≈', r'$\approx$'),           # ≈
    ('≤', r'$\leq$'),              # ≤
    ('≥', r'$\geq$'),              # ≥
    # Dashes
    ('—', '---'),                  # em dash —
    ('–', '--'),                   # en dash –
    # Non-breaking space
    (' ', '~'),
    # Smart quotes
    ('“', '``'),                   # "
    ('”', "''"),                   # "
    ('‘', '`'),                    # '
    ('’', "'"),                    # '
    # Middle dot / bullet
    ('·', r'$\cdot$'),
    ('•', r'\textbullet{}'),
    # Misc math symbols used outside $
    ('α', r'$\alpha$'),
    ('β', r'$\beta$'),
    ('ν', r'$\nu$'),
    ('τ', r'$\tau$'),
    ('γ', r'$\gamma$'),
    ('μ', r'$\mu$'),
    ('ε', r'$\varepsilon$'),
    ('×', r'$\times$'),
    # Check / cross marks (results tables)
    ('✓', r'$\checkmark$'),
    ('✗', r'$\times$'),
    # Superscripts / subscripts not covered by utf8 inputenc
    ('⁻', r'\textsuperscript{$-$}'),
    ('⁰', r'\textsuperscript{0}'),
    ('⁴', r'\textsuperscript{4}'),
    ('⁵', r'\textsuperscript{5}'),
    ('⁶', r'\textsuperscript{6}'),
    ('⁷', r'\textsuperscript{7}'),
    ('⁸', r'\textsuperscript{8}'),
    ('⁹', r'\textsuperscript{9}'),
    ('₀', r'\textsubscript{0}'),
    ('₁', r'\textsubscript{1}'),
    ('₂', r'\textsubscript{2}'),
    ('₃', r'\textsubscript{3}'),
    ('₄', r'\textsubscript{4}'),
]

# Box-drawing characters used in the ASCII tree diagram
BOX_ASCII = [
    ('│', '|'),   # │
    ('├', '|'),   # ├
    ('└', '\\'), # └
    ('─', '-'),   # ─
    ('┬', '+'),   # ┬
    ('┤', '|'),   # ┤
    ('┌', '+'),   # ┌
    ('┐', '+'),   # ┐
    ('┘', '+'),   # ┘
    ('┼', '+'),   # ┼
    ('┴', '+'),   # ┴
]

# ASCII equivalents for Unicode inside lstlisting/verbatim (listings rejects
# raw multibyte UTF-8); applied only in code mode.
CODE_ASCII = [
    ('↔','<->'),('→','->'),('←','<-'),('⇒','=>'),('⇐','<='),('⇔','<=>'),
    ('≪','<<'),('≫','>>'),('≈','~'),('≠','!='),('≤','<='),('≥','>='),
    ('×','x'),('·','.'),('±','+/-'),('∝','~'),('√','sqrt'),('∞','inf'),
    ('α','alpha'),('β','beta'),('ν','nu'),('τ','tau'),('γ','gamma'),('μ','mu'),
    ('ε','eps'),('σ','sigma'),('Δ','Delta'),('θ','theta'),('Θ','Theta'),
    ('λ','lambda'),('ρ','rho'),('φ','phi'),('ψ','psi'),('Ω','Omega'),('π','pi'),
    ('°','deg'),('—','--'),('–','-'),('‰','o/oo'),
    ('“','"'),('”','"'),('‘',"'"),('’',"'"),('…','...'),
    ('⁻','^-'),('⁰','^0'),('¹','^1'),('²','^2'),('³','^3'),
    ('⁴','^4'),('⁵','^5'),('⁶','^6'),('⁷','^7'),('⁸','^8'),('⁹','^9'),
    ('₀','_0'),('₁','_1'),('₂','_2'),('₃','_3'),('₄','_4'),
    ('₅','_5'),('₆','_6'),('₇','_7'),('₈','_8'),('₉','_9'),
]


def apply_unicode(text, in_code=False):
    """Replace Unicode special characters with safe LaTeX or ASCII equivalents."""
    for u, r in BOX_ASCII:
        text = text.replace(u, r)
    if not in_code:
        for u, r in UNICODE_TEXT:
            text = text.replace(u, r)
    else:
        for u, r in CODE_ASCII:
            text = text.replace(u, r)
        # safety net: drop any remaining non-ASCII so listings never sees raw UTF-8
        text = ''.join(c if ord(c) < 128 else '?' for c in text)
    return text


# ---------------------------------------------------------------------------
# Math extraction / restoration  (protects $...$ from inline-markup rules)
# ---------------------------------------------------------------------------
def extract_math(text):
    store = []
    def repl(m):
        store.append(m.group(0))
        return '\x00M%d\x00' % len(store)
    text = re.sub(r'\$[^$\n]+?\$', repl, text)
    return text, store

def restore_math(text, store):
    def repl(m):
        return store[int(m.group(1)) - 1]
    return re.sub(r'\x00M(\d+)\x00', repl, text)


# ---------------------------------------------------------------------------
# Inline markup  **bold**  *italic*  `code`
# ---------------------------------------------------------------------------
def convert_inline(text):
    text = re.sub(r'\*\*\*(.+?)\*\*\*', r'\\textbf{\\textit{\1}}', text)
    text = re.sub(r'\*\*(.+?)\*\*',     r'\\textbf{\1}', text)
    text = re.sub(r'(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)', r'\\textit{\1}', text)
    text = re.sub(r'`([^`]+)`', r'\\texttt{\1}', text)
    # escape LaTeX specials in prose (math already extracted to placeholders);
    # skip any already preceded by a backslash to avoid double-escaping
    text = re.sub(r'(?<!\\)([&%#_])', r'\\\1', text)
    return text


# ---------------------------------------------------------------------------
# Process one ordinary text line
# ---------------------------------------------------------------------------
def process_line(line):
    line = apply_unicode(line)
    text, store = extract_math(line)
    text = convert_inline(text)
    return restore_math(text, store)


# ---------------------------------------------------------------------------
# Headings: protect $...$ from hyperref PDF-string complaints
# ---------------------------------------------------------------------------
def safe_heading(title):
    def repl(m):
        math  = m.group(0)
        plain = re.sub(r'[\\{}$^_]', '', m.group(1))
        return r'\texorpdfstring{' + math + '}{' + plain + '}'
    return re.sub(r'\$(.+?)\$', repl, title)


# ---------------------------------------------------------------------------
# Table → longtable
# ---------------------------------------------------------------------------
def split_row(row):
    return [c.strip() for c in row.strip().strip('|').split('|')]

def convert_table(buf):
    rows = []
    for line in buf:
        line = line.strip()
        if not line.startswith('|'):
            continue
        if re.match(r'^\|[-: |]+\|$', line):
            continue                      # separator row
        rows.append(split_row(line))
    if not rows:
        return ''
    ncols = max(len(r) for r in rows)
    for r in rows:
        while len(r) < ncols:
            r.append('')
    spec = 'l' * ncols
    out  = [r'\begin{center}', r'\begin{longtable}{' + spec + '}', r'\toprule']
    for k, row in enumerate(rows):
        cells = ' & '.join(process_line(c) for c in row)
        out.append(cells + r' \\')
        if k == 0:
            out.append(r'\midrule')
    out += [r'\bottomrule', r'\end{longtable}', r'\end{center}']
    return '\n'.join(out)


# ---------------------------------------------------------------------------
# Main converter
# ---------------------------------------------------------------------------
def convert(md_text):
    lines = md_text.splitlines()
    n     = len(lines)
    out   = []
    i     = 0

    in_item  = False     # itemize
    in_enum  = False     # enumerate
    in_code  = False
    in_disp  = False     # $$ ... $$
    in_tbl   = False
    tbl_buf  = []
    disp_buf = []

    def end_lists():
        nonlocal in_item, in_enum
        if in_item:
            out.append(r'\end{itemize}');   in_item = False
        if in_enum:
            out.append(r'\end{enumerate}'); in_enum = False

    while i < n:
        raw  = lines[i]
        line = raw.rstrip()

        # ── Code fence ──────────────────────────────────────────────────────
        if line.startswith('```'):
            if not in_code:
                end_lists()
                lang = line[3:].strip()
                out.append(r'\begin{lstlisting}' +
                           ('[language=%s]' % lang if lang else ''))
                in_code = True
            else:
                out.append(r'\end{lstlisting}')
                in_code = False
            i += 1; continue

        if in_code:
            out.append(apply_unicode(raw, in_code=True))
            i += 1; continue

        # ── Display math  (bare $$ on its own line) ──────────────────────────
        if line.strip() == '$$':
            if not in_disp:
                end_lists(); in_disp = True; disp_buf = []
            else:
                body = '\n'.join(disp_buf)
                if r'\tag' in body or r'\begin{cases' in body:
                    env = 'equation'
                else:
                    env = 'equation*'
                out.append(r'\begin{' + env + '}')
                out.append(body)
                out.append(r'\end{'   + env + '}')
                in_disp = False; disp_buf = []
            i += 1; continue

        if in_disp:
            disp_buf.append(line); i += 1; continue

        # ── Single-line display  $$body$$  ──────────────────────────────────
        m1 = re.match(r'^\$\$(.+)\$\$\s*$', line)
        if m1:
            end_lists()
            body = m1.group(1)
            env  = 'equation' if (r'\tag' in body or r'\begin{cases' in body) else 'equation*'
            out += [r'\begin{' + env + '}', body, r'\end{' + env + '}']
            i += 1; continue

        # ── Table ────────────────────────────────────────────────────────────
        if line.startswith('|'):
            if not in_tbl:
                end_lists(); in_tbl = True; tbl_buf = []
            tbl_buf.append(line); i += 1; continue
        else:
            if in_tbl:
                out.append(convert_table(tbl_buf)); in_tbl = False; tbl_buf = []

        # ── Horizontal rule ──────────────────────────────────────────────────
        if re.match(r'^-{3,}\s*$', line):
            end_lists()
            out.append(r'\medskip\noindent\rule{\linewidth}{0.4pt}\medskip')
            i += 1; continue

        # ── Heading ──────────────────────────────────────────────────────────
        hm = re.match(r'^(#{1,4})\s+(.+)', line)
        if hm:
            end_lists()
            level  = len(hm.group(1))
            htitle = hm.group(2).strip()

            # Skip the markdown Table-of-Contents block
            if htitle.strip() in ('Table of Contents', 'Contents'):
                i += 1
                while i < n and not lines[i].strip():          # skip blanks
                    i += 1
                while (i < n and lines[i].strip()              # skip entries
                       and not lines[i].strip().startswith('##')):
                    i += 1
                continue

            # Remove anchor fragments  {#section-i}
            htitle = re.sub(r'\s*\{#[^}]+\}', '', htitle)
            htitle = process_line(htitle)
            htitle = safe_heading(htitle)

            cmds = {1: r'\section', 2: r'\section',
                    3: r'\subsection', 4: r'\subsubsection'}
            out.append(cmds[level] + '{' + htitle + '}')
            i += 1; continue

        # ── Blank line ───────────────────────────────────────────────────────
        if not line.strip():
            end_lists(); out.append(''); i += 1; continue

        # ── Unordered list  - item ───────────────────────────────────────────
        um = re.match(r'^\s*[-*]\s+(.*)', line)
        if um:
            if in_enum: out.append(r'\end{enumerate}'); in_enum = False
            if not in_item: out.append(r'\begin{itemize}'); in_item = True
            out.append(r'  \item ' + process_line(um.group(1)))
            i += 1; continue

        # ── Ordered list  1. item ────────────────────────────────────────────
        om = re.match(r'^\s*\d+\.\s+(.*)', line)
        if om:
            if in_item: out.append(r'\end{itemize}'); in_item = False
            if not in_enum: out.append(r'\begin{enumerate}'); in_enum = True
            out.append(r'  \item ' + process_line(om.group(1)))
            i += 1; continue

        # ── Normal paragraph ─────────────────────────────────────────────────
        end_lists()
        out.append(process_line(line))
        i += 1

    # close any still-open environments
    end_lists()
    if in_tbl:   out.append(convert_table(tbl_buf))
    if in_disp:  out += [r'\[', '\n'.join(disp_buf), r'\]']
    if in_code:  out.append(r'\end{lstlisting}')
    return '\n'.join(out)


# ---------------------------------------------------------------------------
# Title / abstract extractor
# ---------------------------------------------------------------------------
def extract_header(md_text):
    lines = md_text.splitlines()
    title, abstract = '', ''
    abs_lines, in_abs = [], False
    remainder = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith('# ') and not title:
            title = s[2:].strip(); remainder = i + 1
        elif s == '## Abstract':
            in_abs = True
        elif in_abs:
            if s.startswith('## ') or s.startswith('**Keywords'):
                in_abs = False; remainder = i; break
            abs_lines.append(s)
    abstract = ' '.join(abs_lines).strip()
    return title, abstract, '\n'.join(lines[remainder:])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # Legacy converter for the extended markdown monograph (reports/).
    # The journal manuscript is now maintained directly in
    # manuscript/OSET_ActaMaterialia.tex (Elsevier elsarticle).
    SRC = 'd:/Work/Others/Padmanabhan/reports/OSET_manuscript.md'
    DST = 'd:/Work/Others/Padmanabhan/reports/OSET_manuscript.tex'

    with open(SRC, encoding='utf-8-sig') as f:
        md = f.read()

    # strip HTML comments (e.g. commented-out abstract) and stray BOMs
    md = re.sub(r'<!--.*?-->', '', md, flags=re.DOTALL)
    md = md.replace('﻿', '')

    title, abstract, body = extract_header(md)

    authors = (r'K.\,A.~Padmanabhan \and '
               r'H.~Gleiter \and '
               r'S.\,V.~Divinski')

    title_block  = '\\title{'  + process_line(title)   + '}\n'
    title_block += '\\author{' + authors               + '}\n'
    title_block += '\\date{}\n'

    abs_block = ''
    if abstract:
        abs_block = ('\\begin{abstract}\n'
                     + process_line(abstract) + '\n'
                     + '\\end{abstract}\n')

    body_tex = convert(body)

    # Replace the numbered reference list with thebibliography
    ref_tag = '\\section{References}'
    if ref_tag in body_tex:
        pre, post = body_tex.split(ref_tag, 1)
        post = post.replace('\\begin{enumerate}',
                            '\\begin{thebibliography}{99}', 1)
        post = post.replace('\\end{enumerate}',
                            '\\end{thebibliography}', 1)
        post = re.sub(r'  \\item ', r'  \\bibitem{ref} ', post)
        body_tex = pre + ref_tag + post

    doc = (PREAMBLE
           + title_block + '\n'
           + '\\maketitle\n\n'
           + abs_block + '\n'
           + '\\tableofcontents\n'
           + '\\newpage\n\n'
           + body_tex
           + POSTAMBLE)

    with open(DST, 'w', encoding='utf-8') as f:
        f.write(doc)

    # Sanity check
    for env in ('equation', 'equation*', 'itemize', 'enumerate',
                'longtable', 'lstlisting', 'document', 'thebibliography'):
        b = doc.count('\\begin{' + env + '}')
        e = doc.count('\\end{'   + env + '}')
        tag = 'OK' if b == e else ('MISMATCH b=%d e=%d' % (b, e))
        print('  %-22s %s' % (env, tag))
    print('\nWritten %d lines -> %s' % (len(doc.splitlines()), DST))
