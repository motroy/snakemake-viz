"""
snakemake_viz_wasm.py — Snakemake workflow visualizer for browser / Pyodide use.

Adapted from snakemake_viz.py: all file-system I/O removed; exposes
visualize_string(content, title) -> full HTML string as the public API.
"""

import re, textwrap
from collections import defaultdict


# ════════════════════════════════════════════════════════════════════════════
# 1. REGEXES & PARSING
# ════════════════════════════════════════════════════════════════════════════

RULE_RE   = re.compile(r'^rule\s+(\w+)\s*:', re.M)
INPUT_RE  = re.compile(
    r'\binput\s*:(.*?)(?=\n[ \t]*(?:output|params|threads|resources|conda|'
    r'shell|run|log|benchmark|wildcard_constraints|message|priority|group|'
    r'envmodules|cache|wrapper|notebook|\Z))', re.S)
OUTPUT_RE = re.compile(
    r'\boutput\s*:(.*?)(?=\n[ \t]*(?:input|params|threads|resources|conda|'
    r'shell|run|log|benchmark|wildcard_constraints|message|priority|group|'
    r'envmodules|cache|wrapper|notebook|\Z))', re.S)
SHELL_RE  = re.compile(
    r'\bshell\s*:(.*?)(?=\n[ \t]*(?:input|output|params|threads|resources|'
    r'conda|log|benchmark|run|\Z))', re.S)
STR_RE    = re.compile(r'["\']([^"\']+)["\']')
PATH_RE   = re.compile(r'["\']([^"\']*\/[^"\']+)["\']')


def extract_paths(block: str) -> list:
    paths = []
    for m in STR_RE.finditer(block):
        s = m.group(1)
        if ('/' in s or '{' in s) and len(s) > 4:
            paths.append(s)
    return paths


def parse_rules_from_string(content: str):
    """Parse rule blocks directly from Snakefile content string."""
    rules, order = {}, []
    positions = [(m.start(), m.group(1)) for m in RULE_RE.finditer(content)]
    for i, (pos, name) in enumerate(positions):
        if name == 'all':
            continue
        end = positions[i + 1][0] if i + 1 < len(positions) else len(content)
        block = content[pos:end]
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


def tool_from_path(path: str):
    clean = re.sub(r'\{[^}]+\}', '*', path)
    parts = [p for p in re.split(r'[/]', clean) if p not in SKIP_PARTS]
    for p in reversed(parts):
        key = p.lower().rstrip('0123456789').rstrip('_-')
        if key in CANONICAL:
            return CANONICAL[key]
    flat = '/'.join(parts).lower()
    for keyword, canon in CANONICAL.items():
        if keyword in flat:
            return canon
    return None


def synthesise_from_string(content: str):
    """Synthesise nodes from rule all output paths when real rules are sparse."""
    m = re.search(r'^rule\s+all\s*:(.*?)(?=^rule\s|\Z)', content, re.M | re.S)
    block = m.group(1) if m else content
    paths = [pm.group(1) for pm in PATH_RE.finditer(block)]
    seen_tools, order = {}, []
    for path in paths:
        tool = tool_from_path(path)
        if tool and tool not in seen_tools:
            seen_tools[tool] = path
            order.append(tool)
    rules = {t: {'inputs': [], 'outputs': [seen_tools[t]], 'shell': ''}
             for t in order}
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


KNOWN_ORDER = [
    'fastp','FastQC','Trimmomatic','Trim Galore',
    'NanoQ','Filtlong','NanoStat',
    'Kraken2','Bracken','Gambit','Melon','Sample Validation',
    'Dragonflye','SPAdes','Flye','Unicycler','Hifiasm',
    'CheckM2','CheckM','Compleasm','BUSCO','QUAST',
    'WG Kraken+Bracken','WG Validation','MLST','rMLST',
    'ABRicate','AMRFinder+','ResFinder','RGI',
    'BWA','Bowtie2','Minimap2','SAMtools',
    'GATK','BCFtools','FreeBayes','DeepVariant',
    'Prokka','Bakta','eggNOG','InterPro',
    'Contigs Dir','chewBBACA','Schema Creation','Allelecall',
    'MST','Distance Matrix','GrapeTree',
    'GAS mclustering',
]


def infer_edges(rules: dict, order: list) -> list:
    names = list(rules.keys())
    out_toks = {n: [stem_tokens(p) for p in rules[n]['outputs']] for n in names}
    edges = set()

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

    if len(edges) < max(1, len(names) // 4):
        idx = {n: KNOWN_ORDER.index(n) if n in KNOWN_ORDER else 9999
               for n in names}
        sorted_names = sorted(names, key=lambda n: (idx[n], order.index(n)
                                                     if n in order else 9999))
        for i in range(len(sorted_names) - 1):
            edges.add((sorted_names[i], sorted_names[i + 1]))

    return list(edges)


# ════════════════════════════════════════════════════════════════════════════
# 4. COLOUR / GROUP HEURISTIC
# ════════════════════════════════════════════════════════════════════════════

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
    names = [n for n in order if n in rules]
    col = {n: 0 for n in names}
    for _ in range(len(names) + 2):
        changed = False
        for a, b in edges:
            if a in col and b in col and col[b] <= col[a]:
                col[b] = col[a] + 1
                changed = True
        if not changed:
            break

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
    BW, BH, BR        = 148, 38, 7
    COL_GAP, ROW_GAP  = 30, 14
    LANE_LBL_W        = 70
    LANE_PAD_X        = 14
    LANE_PAD_Y        = 12
    MARGIN            = 16

    lane_members: dict = defaultdict(list)
    for name, (c, r) in pos.items():
        grp = groups[name][0]
        lane_members[grp].append((c, r, name))

    lane_order = sorted(lane_members,
                        key=lambda g: min(c for c, r, n in lane_members[g]))
    for g in lane_order:
        lane_members[g].sort()

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

    parts = []

    for g in lane_order:
        gy, gh = lane_y[g]
        bg = LANE_BG.get(g, '#f5f6fa')
        parts.append(
            f'<rect x="4" y="{gy}" width="{total_w-8}" height="{gh}" '
            f'rx="8" fill="{bg}" stroke="#d0d6e8" stroke-width="1"/>')
        tx, ty = 4 + LANE_LBL_W / 2, gy + gh / 2
        ge = g.replace('&', '&amp;')
        parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" fill="#5a6480" '
            f'font-family="\'Source Sans 3\',sans-serif" '
            f'font-size="11" font-weight="700" letter-spacing="0.04em" '
            f'text-anchor="middle" dominant-baseline="middle" '
            f'transform="rotate(-90,{tx:.1f},{ty:.1f})">{ge}</text>')

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
        padding:28px 16px 48px;min-height:100vh;text-align:center}}
  h1{{font-size:1.25rem;font-weight:700;color:#1e2535;
      margin-bottom:4px;letter-spacing:-.01em}}
  p.sub{{font-size:.75rem;color:#7a8499;margin-bottom:22px;letter-spacing:.03em}}
  .wrap{{background:#fff;border-radius:12px;
         box-shadow:0 2px 16px rgba(0,0,0,.08);
         padding:24px 20px;overflow-x:auto;
         display:inline-block;max-width:calc(100vw - 32px);
         text-align:left}}
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
# 7. PUBLIC API (called from JavaScript via Pyodide)
# ════════════════════════════════════════════════════════════════════════════

_last_node_count = 0
_last_edge_count = 0


def visualize_string(content: str, title: str = 'Workflow') -> str:
    """
    Main entry point for browser / Pyodide use.
    Takes raw Snakefile text, returns a complete standalone HTML string.
    """
    global _last_node_count, _last_edge_count

    if not title:
        title = 'Workflow'

    rules, order = parse_rules_from_string(content)
    if len(rules) < 2:
        rules, order = synthesise_from_string(content)
    if not rules:
        raise ValueError(
            "No rules or output paths found. "
            "Please check that this is a valid Snakefile."
        )

    edges  = infer_edges(rules, order)
    groups = {n: classify(n, rules[n].get('shell', '')) for n in rules}
    pos    = topo_layout(rules, edges, order)

    _last_node_count = len(rules)
    _last_edge_count = len(edges)

    return build_html(rules, edges, pos, title, groups)
