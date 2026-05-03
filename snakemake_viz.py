#!/usr/bin/env python3
"""
snakemake_viz.py — Auto-generate an nf-core-style workflow diagram from any Snakefile.

Usage:
    python snakemake_viz.py Snakefile [options]
    python snakemake_viz.py Snakefile -o diagram.png
    python snakemake_viz.py Snakefile -o diagram.html   # interactive HTML only

Options:
    -o / --output   Output file (.png or .html)  [default: workflow_diagram.png]
    --title         Diagram title                [default: derived from filename]
    --dpi           PNG pixel scale factor       [default: 2]
    --no-browser    Skip auto-opening the result

Strategy
────────
Pass 1 – parse rule blocks from Snakefile + all included .smk files.
Pass 2 – if <2 real rules found, synthesise "process nodes" from the
          output path patterns listed in `rule all` (handles the common
          case where the top-level Snakefile only has rule all + includes
          pointing to files not present locally).
Edges  – inferred by matching path-stem tokens between producer outputs
          and consumer inputs, plus a shared-prefix heuristic.
Groups – keyword heuristic assigns each node to a coloured swim-lane.
Render – SVG embedded in HTML; screenshot via Playwright for PNG.

Dependencies:
    pip install playwright
    python -m playwright install chromium
"""

import re, sys, os, argparse, textwrap, time, pathlib
from collections import defaultdict


# ════════════════════════════════════════════════════════════════════════════
# 1. PARSING
# ════════════════════════════════════════════════════════════════════════════

INCLUDE_RE = re.compile(r'^\s*include\s*:\s*["\'](.+?)["\']', re.M)
RULE_RE    = re.compile(r'^rule\s+(\w+)\s*:', re.M)
INPUT_RE   = re.compile(
    r'\binput\s*:(.*?)(?=\n[ \t]*(?:output|params|threads|resources|conda|'
    r'shell|run|log|benchmark|wildcard_constraints|message|priority|group|'
    r'envmodules|cache|wrapper|notebook|\Z))', re.S)
OUTPUT_RE  = re.compile(
    r'\boutput\s*:(.*?)(?=\n[ \t]*(?:input|params|threads|resources|conda|'
    r'shell|run|log|benchmark|wildcard_constraints|message|priority|group|'
    r'envmodules|cache|wrapper|notebook|\Z))', re.S)
SHELL_RE   = re.compile(
    r'\bshell\s*:(.*?)(?=\n[ \t]*(?:input|output|params|threads|resources|'
    r'conda|log|benchmark|run|\Z))', re.S)
STR_RE     = re.compile(r'["\']([^"\']+)["\']')
PATH_RE    = re.compile(r'["\']([^"\']*\/[^"\']+)["\']')


def collect_files(snakefile: str) -> list:
    root = pathlib.Path(snakefile).parent
    seen, queue, out = set(), [snakefile], []
    while queue:
        f = queue.pop(0)
        fp = pathlib.Path(f)
        if not fp.exists():
            fp = root / fp
        key = str(fp.resolve()) if fp.exists() else f
        if key in seen:
            continue
        seen.add(key)
        if fp.exists():
            out.append(str(fp))
            try:
                src = fp.read_text(errors='replace')
            except Exception:
                continue
            for inc in INCLUDE_RE.findall(src):
                queue.append(str(root / inc))
    return out


def extract_paths(block: str) -> list:
    paths = []
    for m in STR_RE.finditer(block):
        s = m.group(1)
        if ('/' in s or '{' in s) and len(s) > 4:
            paths.append(s)
    return paths


def parse_rules(snakefile: str):
    """
    Returns (rules_dict, encounter_order).
    rules_dict: name -> {inputs, outputs, shell}
    """
    files = collect_files(snakefile)
    rules, order = {}, []
    for fpath in files:
        try:
            src = pathlib.Path(fpath).read_text(errors='replace')
        except Exception:
            continue
        positions = [(m.start(), m.group(1)) for m in RULE_RE.finditer(src)]
        for i, (pos, name) in enumerate(positions):
            if name == 'all':
                continue
            end = positions[i+1][0] if i+1 < len(positions) else len(src)
            block = src[pos:end]
            inp_m = INPUT_RE.search(block)
            out_m = OUTPUT_RE.search(block)
            sh_m  = SHELL_RE.search(block)
            rules[name] = {
                'inputs':  extract_paths(inp_m.group(1)) if inp_m else [],
                'outputs': extract_paths(out_m.group(1)) if out_m else [],
                'shell':   sh_m.group(1).strip()[:200]   if sh_m  else '',
            }
            if name not in order:
                order.append(name)
    return rules, order


# ════════════════════════════════════════════════════════════════════════════
# 2. FALLBACK: synthesise nodes from rule all output paths
# ════════════════════════════════════════════════════════════════════════════

# Maps a path segment keyword → canonical tool name
CANONICAL = {
    'fastp':'fastp', 'fastqc':'FastQC', 'multiqc':'MultiQC',
    'trimmomatic':'Trimmomatic', 'trim_galore':'Trim Galore',
    'nanoq':'NanoQ', 'filtlong':'Filtlong', 'nanostat':'NanoStat',
    'kraken':'Kraken2', 'bracken':'Bracken', 'gambit':'Gambit',
    'melon':'Melon', 'metaphlan':'MetaPhlAn', 'kaiju':'Kaiju',
    'centrifuge':'Centrifuge', 'sylph':'Sylph',
    'sample_validation':'Sample Validation', 'validation':'Validation',
    'dragonflye':'Dragonflye', 'spades':'SPAdes', 'flye':'Flye',
    'unicycler':'Unicycler', 'hifiasm':'Hifiasm', 'raven':'Raven',
    'miniasm':'Miniasm', 'canu':'Canu', 'wtdbg':'wtdbg2',
    'medaka':'Medaka', 'polypolish':'Polypolish', 'pilon':'Pilon',
    'bwa':'BWA', 'bowtie':'Bowtie2', 'minimap':'Minimap2',
    'samtools':'SAMtools', 'bamtools':'BAMtools',
    'checkm2':'CheckM2', 'checkm':'CheckM', 'compleasm':'Compleasm',
    'busco':'BUSCO', 'quast':'QUAST',
    'mlst':'MLST', 'rmlst':'rMLST',
    'wgkb':'WG Kraken+Bracken', 'wgv':'WG Validation',
    'abricate':'ABRicate', 'amrfinderplus':'AMRFinder+',
    'amrfinder':'AMRFinder+', 'resfinder':'ResFinder', 'rgi':'RGI',
    'prokka':'Prokka', 'bakta':'Bakta', 'eggnog':'eggNOG',
    'interpro':'InterPro',
    'chewbbaca':'chewBBACA', 'chewie':'chewBBACA',
    'grapetree':'GrapeTree', 'mst':'MST', 'distance_matrix':'Distance Matrix',
    'gas_mcluster':'GAS mclustering', 'contigs_dir':'Contigs Dir',
    'schema':'Schema Creation', 'allelecall':'Allelecall',
    'gatk':'GATK', 'bcftools':'BCFtools', 'freebayes':'FreeBayes',
    'deepvariant':'DeepVariant',
    'diamond':'DIAMOND', 'blast':'BLAST', 'hmmer':'HMMER',
}

SKIP_PARTS = {'output','input','assembly','reads','SR','LR','sample',
              'comparisongroup','modules','*',''}


def tool_from_path(path: str) -> str | None:
    """Derive the canonical tool name from a path pattern."""
    clean = re.sub(r'\{[^}]+\}', '*', path)
    parts = [p for p in re.split(r'[/]', clean) if p not in SKIP_PARTS]
    # walk from right: first recognisable part wins
    for p in reversed(parts):
        key = p.lower().rstrip('0123456789').rstrip('_-')
        if key in CANONICAL:
            return CANONICAL[key]
    # second pass: substring match
    flat = '/'.join(parts).lower()
    for keyword, canon in CANONICAL.items():
        if keyword in flat:
            return canon
    return None


def synthesise_from_rule_all(snakefile: str):
    """
    When real rules aren't available, build a node list + output paths
    by parsing the output paths in rule all.
    """
    try:
        src = pathlib.Path(snakefile).read_text(errors='replace')
    except Exception:
        return {}, []

    # Grab everything inside `rule all:` block
    m = re.search(r'^rule\s+all\s*:(.*?)(?=^rule\s|\Z)', src, re.M|re.S)
    block = m.group(1) if m else src

    paths = [m.group(1) for m in PATH_RE.finditer(block)]

    # Build a set of (tool_name, representative_path)
    seen_tools = {}
    order = []
    for path in paths:
        tool = tool_from_path(path)
        if tool and tool not in seen_tools:
            seen_tools[tool] = path
            order.append(tool)

    rules = {}
    for tool in order:
        rules[tool] = {
            'inputs':  [],
            'outputs': [seen_tools[tool]],
            'shell':   '',
        }
    return rules, order


# ════════════════════════════════════════════════════════════════════════════
# 3. EDGE INFERENCE
# ════════════════════════════════════════════════════════════════════════════

def stem_tokens(path: str) -> set:
    clean = re.sub(r'\{[^}]+\}', '*', path)
    parts = re.split(r'[/._\-]', clean)
    return {p.lower() for p in parts
            if len(p) > 2 and p not in ('output','input','sample','*','tsv',
                                         'txt','csv','gz','html','json','png',
                                         'tre','tab','trn','bam','vcf','fa',
                                         'fasta','fastq')}


# Hard-coded logical ordering for known tools (used when path-overlap fails)
KNOWN_ORDER = [
    'fastp','FastQC','Trimmomatic','Trim Galore',        # SR QC
    'NanoQ','Filtlong','NanoStat',                        # LR QC
    'Kraken2','Bracken','Gambit','Melon','Sample Validation',  # reads speciesID
    'Dragonflye','SPAdes','Flye','Unicycler','Hifiasm',   # assembly
    'CheckM2','CheckM','Compleasm','BUSCO','QUAST',        # asm QC
    'WG Kraken+Bracken','WG Validation','MLST','rMLST',   # asm speciesID
    'ABRicate','AMRFinder+','ResFinder','RGI',             # resistome
    'BWA','Bowtie2','Minimap2','SAMtools',                 # mapping
    'GATK','BCFtools','FreeBayes','DeepVariant',           # variants
    'Prokka','Bakta','eggNOG','InterPro',                  # annotation
    'Contigs Dir','chewBBACA','Schema Creation','Allelecall',
    'MST','Distance Matrix','GrapeTree',                   # cgMLST
    'GAS mclustering',
]


def infer_edges(rules: dict, order: list) -> list:
    names = list(rules.keys())
    out_toks = {n: [stem_tokens(p) for p in rules[n]['outputs']] for n in names}

    edges = set()

    # Pass 1: path-token overlap
    for a in names:
        for b in names:
            if a == b:
                continue
            for atoks in out_toks.get(a, []):
                for binp in rules[b]['inputs']:
                    overlap = atoks & stem_tokens(binp) - {
                        'output','reads','assembly','speciesid','resistome',
                        'cgmlst','comparisongroup','contigs'}
                    if len(overlap) >= 2:
                        edges.add((a, b))

    # Pass 2: if very few edges, use encounter order as a linear chain
    # (safe fallback — avoids empty diagram)
    if len(edges) < max(1, len(names)//4):
        # try KNOWN_ORDER positions
        idx = {n: KNOWN_ORDER.index(n) if n in KNOWN_ORDER else 9999
               for n in names}
        sorted_names = sorted(names, key=lambda n: (idx[n], order.index(n)
                                                     if n in order else 9999))
        for i in range(len(sorted_names)-1):
            edges.add((sorted_names[i], sorted_names[i+1]))

    return list(edges)


# ════════════════════════════════════════════════════════════════════════════
# 4. COLOUR / GROUP HEURISTIC
# ════════════════════════════════════════════════════════════════════════════

# (keywords, group_label, fill, stroke, text_color)
KEYWORD_GROUPS = [
    (['qc','fastp','fastqc','multiqc','trimmomatic','trim','nanoq','filtlong',
      'nanostat','checkm','compleasm','busco','quast'],
     'QC',           '#fff8e8','#c9960a','#6b4700'),
    (['kraken','bracken','metaphlan','centrifuge','gambit','melon','kaiju',
      'diamond','sylph','validation','validate','mlst','rmlst',
      'wgkb','wgv','taxon','species','classify'],
     'Species ID',   '#edfaf4','#2fa868','#0e5532'),
    (['assembly','dragonflye','spades','flye','unicycler','hifiasm','raven',
      'canu','wtdbg','miniasm','medaka','polypolish','pilon','contigs'],
     'Assembly',     '#fff4ee','#d4631a','#6b2800'),
    (['abricate','amrfinder','resfinder','card','vfdb','rgi','resistome',
      'amr','resistance','virulence'],
     'Resistome',    '#fff0f0','#d43535','#6b0000'),
    (['cgmlst','chewbbaca','chewie','grapetree','mst','phylo','tree',
      'distance','cluster','gas','schema','allelecall','trainingfiles'],
     'cgMLST',       '#f8f0ff','#8b3fd4','#3a0070'),
    (['align','bwa','bowtie','minimap','samtools','bamtools','mapping'],
     'Mapping',      '#edf6ff','#1a7fc4','#00305a'),
    (['variant','snp','gatk','bcftools','freebayes','deepvariant','vcf'],
     'Variants',     '#fff0f8','#c43580','#60003a'),
    (['annotate','prokka','bakta','annotation','eggnog','interpro','hmmer',
      'blast','diamond'],
     'Annotation',   '#f4fbec','#5a9e18','#244a00'),
]

DEFAULT_GROUP = ('Other', '#f0f2f8', '#5a6480', '#1e2535')


def classify(name: str, shell: str = '') -> tuple:
    key = (name + ' ' + shell).lower()
    for keywords, grp, fill, stroke, text in KEYWORD_GROUPS:
        if any(k in key for k in keywords):
            return grp, fill, stroke, text
    return DEFAULT_GROUP


# ════════════════════════════════════════════════════════════════════════════
# 5. TOPOLOGICAL LAYOUT
# ════════════════════════════════════════════════════════════════════════════

def topo_layout(rules: dict, edges: list, order: list) -> dict:
    """Returns pos: name -> (col, row)"""
    names = [n for n in order if n in rules]
    # longest-path column assignment
    col = {n: 0 for n in names}
    for _ in range(len(names) + 2):
        changed = False
        for a, b in edges:
            if a in col and b in col and col[b] <= col[a]:
                col[b] = col[a] + 1
                changed = True
        if not changed:
            break

    # pack into rows per column, keeping encounter order
    col_rows: dict = defaultdict(list)
    for n in names:
        col_rows[col[n]].append(n)

    pos = {}
    for c, ns in col_rows.items():
        for r, n in enumerate(ns):
            pos[n] = (c, r)
    return pos


# ════════════════════════════════════════════════════════════════════════════
# 6. HTML / SVG BUILDER
# ════════════════════════════════════════════════════════════════════════════

LANE_BG = {
    'QC':         '#fffcf2',
    'Species ID': '#f2fdf7',
    'Assembly':   '#fff8f4',
    'Resistome':  '#fff5f5',
    'cgMLST':     '#f8f2ff',
    'Mapping':    '#f2f8ff',
    'Variants':   '#fff4fc',
    'Annotation': '#f4faee',
    'Other':      '#f5f6fa',
}


def build_html(rules: dict, edges: list, pos: dict, title: str,
               groups: dict) -> str:
    BW, BH, BR         = 148, 38, 7
    COL_GAP, ROW_GAP    = 30, 14
    LANE_LBL_W          = 70
    LANE_PAD_X          = 14
    LANE_PAD_Y          = 12
    MARGIN              = 16

    # ── swim-lane assignment ─────────────────────────────────────────────
    lane_members: dict = defaultdict(list)  # group -> [(col, row, name)]
    for name, (c, r) in pos.items():
        grp = groups[name][0]
        lane_members[grp].append((c, r, name))

    lane_order = sorted(lane_members,
                        key=lambda g: min(c for c,r,n in lane_members[g]))
    for g in lane_order:
        lane_members[g].sort()

    # re-assign local rows within each lane (per column)
    local_pos: dict = {}
    lane_max_rows: dict = {}
    for g in lane_order:
        col_cnt: dict = defaultdict(int)
        for c, r, name in lane_members[g]:
            local_pos[name] = (c, col_cnt[c])
            col_cnt[c] += 1
        lane_max_rows[g] = max(col_cnt.values(), default=1)

    lane_y: dict = {}
    cur_y = MARGIN
    for g in lane_order:
        h = (lane_max_rows[g] * (BH + ROW_GAP) - ROW_GAP + 2 * LANE_PAD_Y)
        lane_y[g] = (cur_y, h)
        cur_y += h + 10

    total_h = cur_y + MARGIN
    max_col = max((c for c, r in pos.values()), default=0)
    COL_W = BW + COL_GAP
    CONTENT_X = LANE_LBL_W + LANE_PAD_X
    total_w = CONTENT_X + (max_col + 1) * COL_W + MARGIN

    def box_xy(name):
        g = groups[name][0]
        c, lr = local_pos[name]
        gy, gh = lane_y[g]
        x = CONTENT_X + c * COL_W
        y = gy + LANE_PAD_Y + lr * (BH + ROW_GAP)
        return x, y

    coords = {n: box_xy(n) for n in rules if n in local_pos}

    # ── SVG elements ─────────────────────────────────────────────────────
    parts = []

    # lane backgrounds
    for g in lane_order:
        gy, gh = lane_y[g]
        bg = LANE_BG.get(g, '#f5f6fa')
        parts.append(
            f'<rect x="4" y="{gy}" width="{total_w-8}" height="{gh}" '
            f'rx="8" fill="{bg}" stroke="#d0d6e8" stroke-width="1"/>')
        tx, ty = 4 + LANE_LBL_W/2, gy + gh/2
        ge = g.replace('&', '&amp;')
        parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" fill="#5a6480" '
            f'font-family="\'Source Sans 3\',sans-serif" '
            f'font-size="11" font-weight="700" letter-spacing="0.04em" '
            f'text-anchor="middle" dominant-baseline="middle" '
            f'transform="rotate(-90,{tx:.1f},{ty:.1f})">{ge}</text>')

    # arrowhead markers
    parts.append(textwrap.dedent("""
      <defs>
        <marker id="arr" markerWidth="7" markerHeight="7"
                refX="6" refY="3.5" orient="auto">
          <polygon points="0 0,7 3.5,0 7" fill="#8a96ad"/>
        </marker>
        <marker id="arrd" markerWidth="7" markerHeight="7"
                refX="6" refY="3.5" orient="auto">
          <polygon points="0 0,7 3.5,0 7" fill="#b0a0d4"/>
        </marker>
      </defs>"""))

    # edges
    for a, b in edges:
        if a not in coords or b not in coords:
            continue
        ax, ay = coords[a]
        bx, by = coords[b]
        same_lane = groups[a][0] == groups[b][0]
        dashed = not same_lane
        marker = 'arrd' if dashed else 'arr'
        stroke = '#b0a0d4' if dashed else '#8a96ad'
        da = 'stroke-dasharray="5,3"' if dashed else ''

        ar_x, ar_y = ax + BW, ay + BH // 2
        al_x, al_y = bx, by + BH // 2
        ab_x, ab_y = ax + BW // 2, ay + BH
        bt_x, bt_y = bx + BW // 2, by

        if abs(ay - by) < 4 and bx > ax:
            d = f'M{ar_x},{ar_y} L{al_x-6},{al_y}'
        elif bx > ax:
            mx = ar_x + (al_x - ar_x) * 0.5
            d = (f'M{ar_x},{ar_y} H{mx:.0f} '
                 f'V{al_y:.0f} H{al_x-6:.0f}')
        else:
            d = (f'M{ab_x},{ab_y} V{ab_y+12:.0f} '
                 f'H{bt_x:.0f} V{bt_y+6:.0f}')

        parts.append(
            f'<path d="{d}" fill="none" stroke="{stroke}" '
            f'stroke-width="1.4" {da} marker-end="url(#{marker})" '
            f'stroke-linecap="round" stroke-linejoin="round" opacity="0.85"/>')

    # boxes
    for name, (bx, by) in coords.items():
        grp, fill, stroke, txt_col = groups[name]
        label = name.replace('_', ' ')
        if len(label) > 20:
            label = label[:19] + '…'
        le = label.replace('&', '&amp;').replace('<', '&lt;')
        data = rules[name]
        outs = '; '.join(data['outputs'][:3]) or '—'
        tip = (name + ': ' + outs[:140]).replace('"', '&quot;')
        parts.append(
            f'<g class="rbox" data-tip="{tip}">'
            f'<rect x="{bx}" y="{by}" width="{BW}" height="{BH}" '
            f'rx="{BR}" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>'
            f'<text x="{bx+BW//2}" y="{by+BH//2}" fill="{txt_col}" '
            f'font-family="\'Source Sans 3\',sans-serif" '
            f'font-size="11.5" font-weight="600" '
            f'text-anchor="middle" dominant-baseline="middle">'
            f'{le}</text></g>')

    svg_body = '\n'.join(parts)
    title_e = title.replace('&', '&amp;').replace('<', '&lt;')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title_e}</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@400;600;700&display=swap');
  *{{box-sizing:border-box;margin:0;padding:0}}
  body{{background:#f5f6f8;font-family:'Source Sans 3',sans-serif;
        display:flex;flex-direction:column;align-items:center;
        padding:28px 16px 48px;min-height:100vh}}
  h1{{font-size:1.25rem;font-weight:700;color:#1e2535;
      margin-bottom:4px;letter-spacing:-.01em}}
  p.sub{{font-size:.75rem;color:#7a8499;margin-bottom:22px;letter-spacing:.03em}}
  .wrap{{background:#fff;border-radius:12px;
         box-shadow:0 2px 16px rgba(0,0,0,.08);
         padding:24px 20px;overflow-x:auto}}
  svg{{display:block}}
  .rbox{{cursor:default}}
  .rbox:hover rect{{filter:brightness(.93)}}
  #tip{{position:fixed;pointer-events:none;background:#1e2535;color:#e8ecf4;
        font-family:'Source Sans 3',sans-serif;font-size:12px;line-height:1.5;
        padding:8px 12px;border-radius:6px;max-width:280px;opacity:0;
        transition:opacity .1s;z-index:99;
        box-shadow:0 4px 12px rgba(0,0,0,.25)}}
  #tip.on{{opacity:1}}
  #tip strong{{display:block;color:#fff;font-size:12.5px;margin-bottom:3px}}
</style>
</head>
<body>
<h1>{title_e}</h1>
<p class="sub">Auto-generated Snakemake workflow diagram</p>
<div class="wrap">
<svg id="diagram" xmlns="http://www.w3.org/2000/svg"
     width="{total_w}" height="{total_h}">
{svg_body}
</svg>
</div>
<div id="tip"></div>
<script>
const tip = document.getElementById('tip');
document.querySelectorAll('.rbox').forEach(g => {{
  const t = g.dataset.tip || '';
  const [name, ...rest] = t.split(': ');
  g.addEventListener('mouseenter', () => {{
    tip.innerHTML = '<strong>' + name + '</strong>' + (rest.join(': ')||'');
    tip.classList.add('on');
  }});
  g.addEventListener('mousemove', e => {{
    tip.style.left = (e.clientX + 14) + 'px';
    tip.style.top  = (e.clientY - 40) + 'px';
  }});
  g.addEventListener('mouseleave', () => tip.classList.remove('on'));
}});
</script>
</body>
</html>"""


# ════════════════════════════════════════════════════════════════════════════
# 7. PNG RENDERING VIA PLAYWRIGHT
# ════════════════════════════════════════════════════════════════════════════

def render_png(html_path: str, png_path: str, dpi: int = 2):
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        sys.exit("ERROR: playwright not installed.\n"
                 "Run: pip install playwright && python -m playwright install chromium")
    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={'width': 2600, 'height': 1400})
        page.goto(f'file://{pathlib.Path(html_path).resolve()}',
                  wait_until='networkidle')
        time.sleep(1.5)  # let fonts load
        clip = page.evaluate("""() => {
            const el = document.querySelector('.wrap');
            const r = el.getBoundingClientRect();
            return {x:r.left-4, y:r.top-4,
                    width:r.width+8, height:r.height+8};
        }""")
        page.screenshot(path=png_path, clip=clip, full_page=False,
                        scale='device' if dpi >= 2 else 'css')
        browser.close()


# ════════════════════════════════════════════════════════════════════════════
# 8. MAIN
# ════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description='Auto-generate an nf-core-style diagram from a Snakefile.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python snakemake_viz.py Snakefile
              python snakemake_viz.py Snakefile -o diagram.png --title "My Pipeline"
              python snakemake_viz.py workflow.smk -o diagram.html
        """))
    ap.add_argument('snakefile')
    ap.add_argument('-o', '--output', default='workflow_diagram.png')
    ap.add_argument('--title', default='')
    ap.add_argument('--dpi', type=int, default=2)
    ap.add_argument('--no-browser', action='store_true')
    args = ap.parse_args()

    sf = args.snakefile
    if not pathlib.Path(sf).exists():
        sys.exit(f"ERROR: not found: {sf}")

    title = args.title or pathlib.Path(sf).stem.replace('_', ' ').title()

    print(f"[1/4] Parsing  {sf} …")
    rules, order = parse_rules(sf)

    if len(rules) < 2:
        print(f"      Only {len(rules)} explicit rule(s) found in accessible files.")
        print("      Falling back to synthesising nodes from rule all output paths …")
        rules, order = synthesise_from_rule_all(sf)

    if not rules:
        sys.exit("ERROR: No rules or output paths found. Is this a valid Snakefile?")

    print(f"      {len(rules)} nodes: {', '.join(list(rules)[:10])}"
          f"{'…' if len(rules) > 10 else ''}")

    print("[2/4] Inferring edges …")
    edges = infer_edges(rules, order)
    print(f"      {len(edges)} connections")

    groups = {n: classify(n, rules[n].get('shell', '')) for n in rules}

    print("[3/4] Computing layout …")
    pos = topo_layout(rules, edges, order)

    print("[4/4] Rendering …")
    html = build_html(rules, edges, pos, title, groups)

    out = pathlib.Path(args.output)
    if out.suffix.lower() == '.html':
        out.write_text(html, encoding='utf-8')
        print(f"✓ Saved HTML → {out}")
        if not args.no_browser:
            import webbrowser
            webbrowser.open(out.resolve().as_uri())
    else:
        tmp = out.with_suffix('._tmp.html')
        tmp.write_text(html, encoding='utf-8')
        render_png(str(tmp), str(out), dpi=args.dpi)
        tmp.unlink(missing_ok=True)
        print(f"✓ Saved PNG  → {out}")
        if not args.no_browser:
            import webbrowser
            webbrowser.open(out.resolve().as_uri())


if __name__ == '__main__':
    main()
