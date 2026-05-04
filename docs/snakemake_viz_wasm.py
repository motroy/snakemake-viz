"""
snakemake_viz_wasm.py — Snakemake workflow visualizer for browser / Pyodide use.

Adapted from snakemake_viz.py: all file-system I/O removed; exposes
visualize_string(content, title) -> full HTML string as the public API.
"""

import re, textwrap
from collections import defaultdict


# ════════════════════════════════════════════════════════════════════════════
# 0. ARCHIVE EXTRACTION  (zip / tar — called from JS via Pyodide)
# ════════════════════════════════════════════════════════════════════════════

def _is_snakefile_path(name: str) -> bool:
    base = name.replace('\\', '/').split('/')[-1]
    return (base in ('Snakefile', 'snakefile') or
            base.endswith('.smk') or base.endswith('.snakefile'))


def extract_snakefiles_from_archive(data, filename: str) -> str:
    """
    Extract and concatenate all Snakefile/.smk content from a zip or tar archive.
    `data` is a bytes-like object (Uint8Array from JS).
    """
    import io
    buf   = io.BytesIO(bytes(data))
    parts = []
    fname = filename.lower()

    if fname.endswith('.zip'):
        import zipfile
        try:
            with zipfile.ZipFile(buf) as zf:
                for name in sorted(zf.namelist()):
                    if not name.endswith('/') and _is_snakefile_path(name):
                        try:
                            parts.append(zf.read(name).decode('utf-8', errors='replace'))
                        except Exception:
                            pass
        except zipfile.BadZipFile as exc:
            raise ValueError(f'Invalid zip archive: {exc}') from exc

    elif any(fname.endswith(ext)
             for ext in ('.tar.gz', '.tgz', '.tar.bz2', '.tar.xz', '.tar')):
        import tarfile
        try:
            with tarfile.open(fileobj=buf) as tf:
                for member in sorted(tf.getmembers(), key=lambda m: m.name):
                    if member.isfile() and _is_snakefile_path(member.name):
                        try:
                            fobj = tf.extractfile(member)
                            if fobj:
                                parts.append(fobj.read().decode('utf-8', errors='replace'))
                        except Exception:
                            pass
        except tarfile.TarError as exc:
            raise ValueError(f'Invalid tar archive: {exc}') from exc
    else:
        parts.append(bytes(data).decode('utf-8', errors='replace'))

    if not parts:
        raise ValueError('No Snakefile or .smk files found in the archive.')
    return '\n\n'.join(parts)


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
               groups: dict, direction: str = 'LR') -> str:
    BW, BH, BR        = 148, 38, 7
    COL_GAP, ROW_GAP  = 30, 14
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

    if direction == 'TB':
        # Vertical strips — nodes flow top-to-bottom by DAG depth
        LANE_LBL_H = 26

        lane_x: dict = {}
        cur_x = MARGIN
        for g in lane_order:
            w = lane_max_rows[g] * (BW + COL_GAP) - COL_GAP + 2 * LANE_PAD_X
            lane_x[g] = (cur_x, w)
            cur_x += w + 10

        total_w   = cur_x + MARGIN
        max_col   = max((c for c, r in pos.values()), default=0)
        ROW_H     = BH + ROW_GAP
        CONTENT_Y = LANE_LBL_H + LANE_PAD_Y
        total_h   = CONTENT_Y + (max_col + 1) * ROW_H + MARGIN

        def box_xy(name):
            g = groups[name][0]
            c, lr = local_pos[name]
            gx, gw = lane_x[g]
            x = gx + LANE_PAD_X + lr * (BW + COL_GAP)
            y = CONTENT_Y + c * ROW_H
            return x, y

    else:
        # Horizontal strips — nodes flow left-to-right by DAG depth (default)
        LANE_LBL_W = 70

        lane_y: dict = {}
        cur_y = MARGIN
        for g in lane_order:
            h = lane_max_rows[g] * (BH + ROW_GAP) - ROW_GAP + 2 * LANE_PAD_Y
            lane_y[g] = (cur_y, h)
            cur_y += h + 10

        total_h   = cur_y + MARGIN
        max_col   = max((c for c, r in pos.values()), default=0)
        COL_W     = BW + COL_GAP
        CONTENT_X = LANE_LBL_W + LANE_PAD_X
        total_w   = CONTENT_X + (max_col + 1) * COL_W + MARGIN

        def box_xy(name):
            g = groups[name][0]
            c, lr = local_pos[name]
            gy, gh = lane_y[g]
            x = CONTENT_X + c * COL_W
            y = gy + LANE_PAD_Y + lr * (BH + ROW_GAP)
            return x, y

    coords = {n: box_xy(n) for n in rules if n in local_pos}

    parts = []

    # Lane backgrounds + labels
    for g in lane_order:
        bg = LANE_BG.get(g, '#f5f6fa')
        ge = g.replace('&', '&amp;')
        if direction == 'TB':
            gx, gw = lane_x[g]
            parts.append(
                f'<rect x="{gx}" y="4" width="{gw}" height="{total_h-8}" '
                f'rx="8" fill="{bg}" stroke="#d0d6e8" stroke-width="1"/>')
            parts.append(
                f'<text x="{gx + gw/2:.1f}" y="{LANE_LBL_H/2:.1f}" fill="#5a6480" '
                f'font-family="\'Source Sans 3\',sans-serif" '
                f'font-size="11" font-weight="700" letter-spacing="0.04em" '
                f'text-anchor="middle" dominant-baseline="middle">{ge}</text>')
        else:
            gy, gh = lane_y[g]
            parts.append(
                f'<rect x="4" y="{gy}" width="{total_w-8}" height="{gh}" '
                f'rx="8" fill="{bg}" stroke="#d0d6e8" stroke-width="1"/>')
            tx, ty = 4 + LANE_LBL_W / 2, gy + gh / 2
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

        if direction == 'TB':
            sx, sy = ax + BW // 2, ay + BH      # source bottom-centre
            tx, ty = bx + BW // 2, by           # target top-centre
            if abs(sx - tx) < 4 and ty > sy:
                d = f'M{sx},{sy} L{tx},{ty-6}'
            elif ty > sy:
                my = sy + (ty - sy) // 2
                d = f'M{sx},{sy} V{my} H{tx} V{ty-6}'
            else:
                d = f'M{sx},{sy} V{sy+12} H{tx} V{ty-6}'
        else:
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
  body{{background:#f5f6f8;font-family:'Source Sans 3',sans-serif;min-height:100vh}}
  h1{{font-size:1.25rem;font-weight:700;color:#1e2535;margin-bottom:4px;letter-spacing:-.01em}}
  p.sub{{font-size:.75rem;color:#7a8499;margin-bottom:22px;letter-spacing:.03em}}
  #main{{text-align:center;padding:20px 16px 48px}}
  .wrap{{background:#fff;border-radius:12px;box-shadow:0 2px 16px rgba(0,0,0,.08);
         padding:24px 20px;overflow-x:auto;display:inline-block;
         max-width:calc(100vw - 32px);text-align:left}}
  svg{{display:block}}
  .rbox{{cursor:default}}
  .rbox:hover rect{{filter:brightness(.93)}}
  #tip{{position:fixed;pointer-events:none;background:#1e2535;color:#e8ecf4;
        font-family:'Source Sans 3',sans-serif;font-size:12px;line-height:1.5;
        padding:8px 12px;border-radius:6px;max-width:280px;opacity:0;
        transition:opacity .1s;z-index:99;box-shadow:0 4px 12px rgba(0,0,0,.25)}}
  #tip.on{{opacity:1}}
  #tip strong{{display:block;color:#fff;font-size:12.5px;margin-bottom:3px}}
  #ctrl{{position:sticky;top:0;z-index:10;background:rgba(245,246,248,.96);
         backdrop-filter:blur(4px);border-bottom:1px solid #e2e8f0;
         padding:7px 16px;display:flex;align-items:center;gap:6px;flex-wrap:wrap}}
  .cbtn{{display:inline-flex;align-items:center;padding:4px 10px;
         border:1.5px solid #d0d6e4;border-radius:5px;background:#fff;
         font-size:.78rem;font-weight:600;color:#3a4460;cursor:pointer;
         font-family:inherit;transition:border-color .12s,color .12s}}
  .cbtn:hover{{border-color:#2563eb;color:#2563eb}}
  #zp{{font-size:.78rem;min-width:3.2em;text-align:center;color:#5a6480;font-family:inherit}}
  .csep{{width:1px;height:18px;background:#d0d6e4;margin:0 3px;flex-shrink:0}}
  .clbl{{font-size:.75rem;font-weight:600;color:#94a3b8;letter-spacing:.03em}}
</style>
</head>
<body>
<div id="ctrl">
  <span class="clbl">Zoom</span>
  <button class="cbtn" onclick="zBy(-0.2)" title="Zoom out">−</button>
  <span id="zp">100%</span>
  <button class="cbtn" onclick="zBy(0.2)" title="Zoom in">+</button>
  <button class="cbtn" onclick="zReset()" title="Reset zoom">↺</button>
  <div class="csep"></div>
  <button class="cbtn" onclick="dlSvg()">↓ SVG</button>
  <button class="cbtn" onclick="dlPng()">↓ PNG</button>
</div>
<div id="main">
<h1>{title_e}</h1>
<p class="sub">Auto-generated Snakemake workflow diagram</p>
<div class="wrap">
<svg id="diagram" xmlns="http://www.w3.org/2000/svg"
     width="{total_w}" height="{total_h}">
{svg_body}
</svg>
</div>
</div>
<div id="tip"></div>
<script>
const svg = document.getElementById('diagram');
const ow = +svg.getAttribute('width'), oh = +svg.getAttribute('height');
svg.setAttribute('viewBox', '0 0 ' + ow + ' ' + oh);
let z = 1;
function zBy(d){{ setZ(z + d); }}
function zReset(){{ setZ(1); }}
function setZ(v){{
  z = Math.max(0.1, Math.min(5, v));
  svg.setAttribute('width',  Math.round(ow * z));
  svg.setAttribute('height', Math.round(oh * z));
  document.getElementById('zp').textContent = Math.round(z * 100) + '%';
}}
function dlSvg(){{
  const cl = svg.cloneNode(true);
  cl.setAttribute('width', ow); cl.setAttribute('height', oh);
  const b = new Blob([new XMLSerializer().serializeToString(cl)], {{type:'image/svg+xml'}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(b); a.download = 'workflow_diagram.svg';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}
function dlPng(){{
  const cl = svg.cloneNode(true);
  cl.setAttribute('width', ow); cl.setAttribute('height', oh);
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(cl)], {{type:'image/svg+xml'}}));
  const img = new Image();
  img.onload = () => {{
    const c = document.createElement('canvas');
    c.width = ow * 2; c.height = oh * 2;
    const ctx = c.getContext('2d');
    ctx.scale(2, 2); ctx.fillStyle = '#f5f6f8'; ctx.fillRect(0, 0, ow, oh);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    c.toBlob(blob => {{
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob); a.download = 'workflow_diagram.png';
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
    }});
  }};
  img.src = url;
}}
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
# 7. MERMAID / EXCALIDRAW BUILDER
# ════════════════════════════════════════════════════════════════════════════

def _mermaid_id(name: str) -> str:
    return re.sub(r'[^a-zA-Z0-9_]', '_', name)


def _build_mermaid_def(rules: dict, edges: list, groups: dict,
                       order: list, direction: str) -> str:
    dir_str = 'LR' if direction == 'LR' else 'TD'

    group_members: dict = defaultdict(list)
    for name in order:
        if name in rules:
            group_members[groups[name][0]].append(name)

    lines = [f'graph {dir_str}']
    for grp, members in group_members.items():
        grp_id    = re.sub(r'[^a-zA-Z0-9]', '_', grp)
        grp_label = grp.replace('"', "'")
        lines.append(f'    subgraph {grp_id}["{grp_label}"]')
        for name in members:
            nid   = _mermaid_id(name)
            label = name.replace('"', "'").replace('[', '(').replace(']', ')')
            lines.append(f'        {nid}["{label}"]')
        lines.append('    end')

    for a, b in edges:
        if a in rules and b in rules:
            lines.append(f'    {_mermaid_id(a)} --> {_mermaid_id(b)}')

    return '\n'.join(lines)


def _build_mermaid_html(mermaid_def: str, title: str) -> str:
    import json as _json
    title_e = title.replace('&', '&amp;').replace('<', '&lt;')
    mermaid_def_js = _json.dumps(mermaid_def)  # safe JS string literal
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title_e}</title>
<script src="https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.min.js"></script>
<link href="https://fonts.googleapis.com/css2?family=Caveat:wght@400;600;700&display=swap" rel="stylesheet">
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    background: #fafaf8;
    font-family: 'Caveat', cursive;
    min-height: 100vh;
    display: flex;
    flex-direction: column;
    align-items: center;
    padding: 0 16px 48px;
  }}
  h1 {{
    font-family: 'Caveat', cursive;
    font-size: 2rem;
    font-weight: 700;
    color: #1e293b;
    margin: 20px 0 4px;
    text-align: center;
  }}
  p.sub {{
    font-size: 1rem;
    color: #7a8499;
    margin-bottom: 24px;
    text-align: center;
  }}
  .mermaid {{
    background: #fff;
    border: 2px solid #d0d4dc;
    border-radius: 14px;
    padding: 28px 24px;
    box-shadow: 4px 4px 0 #d0d4dc;
    overflow-x: auto;
    max-width: calc(100vw - 32px);
  }}
  #ctrl {{
    position: sticky; top: 0; z-index: 10;
    background: rgba(250,250,248,.96);
    backdrop-filter: blur(4px);
    border-bottom: 2px solid #e2e8f0;
    padding: 7px 16px;
    display: flex; align-items: center; gap: 6px; flex-wrap: wrap;
    width: 100%;
    margin-bottom: 20px;
  }}
  .cbtn {{
    display: inline-flex; align-items: center; padding: 3px 12px;
    border: 2px solid #333; border-radius: 5px; background: #fff;
    font-family: 'Caveat', cursive; font-size: 1.05rem; font-weight: 600;
    color: #333; cursor: pointer; box-shadow: 2px 2px 0 #333;
    transition: transform .1s, box-shadow .1s;
  }}
  .cbtn:hover {{ transform: translate(-1px,-1px); box-shadow: 3px 3px 0 #333; }}
  .cbtn:active {{ transform: translate(1px,1px); box-shadow: 1px 1px 0 #333; }}
  #zp {{ font-family:'Caveat',cursive; font-size:1rem; min-width:3.2em;
         text-align:center; color:#5a6480; font-weight:600; }}
  .csep {{ width:1px; height:18px; background:#d0d6e4; margin:0 3px; flex-shrink:0; }}
  .clbl {{ font-size:.85rem; font-weight:700; color:#94a3b8; letter-spacing:.03em; }}
</style>
</head>
<body>
<div id="ctrl">
  <span class="clbl">Zoom</span>
  <button class="cbtn" onclick="zBy(-0.2)" title="Zoom out">−</button>
  <span id="zp">100%</span>
  <button class="cbtn" onclick="zBy(0.2)" title="Zoom in">+</button>
  <button class="cbtn" onclick="zReset()" title="Reset zoom">↺</button>
  <div class="csep"></div>
  <button class="cbtn" onclick="dlSvg()">↓ SVG</button>
  <button class="cbtn" onclick="dlPng()">↓ PNG</button>
  <button class="cbtn" onclick="dlCode()">↓ Code</button>
</div>
<h1>{title_e}</h1>
<p class="sub">Auto-generated Snakemake workflow · Excalidraw style</p>
<div class="mermaid" id="mermaid-wrap">
%%{{init: {{'look': 'handDrawn', 'theme': 'default', 'themeVariables': {{'fontFamily': "'Caveat', cursive", 'fontSize': '16px'}}}}}}%%
{mermaid_def}
</div>
<script>
const _MERMAID_DEF = {mermaid_def_js};

mermaid.initialize({{
  startOnLoad: true,
  look: 'handDrawn',
  theme: 'default',
  themeVariables: {{ fontFamily: "'Caveat', cursive", fontSize: '16px' }}
}});

let z = 1;
const wrap = document.getElementById('mermaid-wrap');

function zBy(d) {{ setZ(z + d); }}
function zReset() {{ setZ(1); }}
function setZ(v) {{
  z = Math.max(0.1, Math.min(5, v));
  wrap.style.transform = 'scale(' + z + ')';
  wrap.style.transformOrigin = 'top center';
  wrap.style.marginBottom = Math.round((z - 1) * wrap.scrollHeight) + 'px';
  document.getElementById('zp').textContent = Math.round(z * 100) + '%';
}}

function _getSvgCloneWithFills() {{
  const svg = document.querySelector('.mermaid svg');
  if (!svg) return null;

  // Snapshot computed fills for native SVG text nodes BEFORE cloning.
  const fills = [];
  svg.querySelectorAll('text, tspan').forEach(el => {{
    const f = window.getComputedStyle(el).getPropertyValue('fill');
    fills.push((!f || f === 'none' || f === 'rgba(0, 0, 0, 0)') ? '#1a1a1a' : f);
  }});

  const clone = svg.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');

  // Apply snapshotted fills to native SVG text elements.
  [...clone.querySelectorAll('text, tspan')].forEach((el, i) => {{
    if (fills[i] !== undefined) el.setAttribute('fill', fills[i]);
  }});

  // Mermaid 11 renders node labels inside <foreignObject> (HTML-in-SVG).
  // foreignObject does NOT render in standalone SVG files or when an SVG
  // is loaded via img.src for canvas drawing, which breaks both SVG and PNG
  // downloads.  Replace each one with a plain SVG <text> element so the
  // exported file is fully self-contained.
  clone.querySelectorAll('foreignObject').forEach(fo => {{
    const label = fo.textContent.trim();
    const fw = parseFloat(fo.getAttribute('width'))  || 100;
    const fh = parseFloat(fo.getAttribute('height')) || 30;
    const fx = parseFloat(fo.getAttribute('x'))      || 0;
    const fy = parseFloat(fo.getAttribute('y'))      || 0;
    const textEl = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    textEl.setAttribute('x', String(fx + fw / 2));
    textEl.setAttribute('y', String(fy + fh / 2));
    textEl.setAttribute('text-anchor', 'middle');
    textEl.setAttribute('dominant-baseline', 'middle');
    textEl.setAttribute('fill', '#1a1a1a');
    textEl.setAttribute('font-size', '14');
    textEl.setAttribute('font-family', 'sans-serif');
    textEl.textContent = label;
    fo.parentNode.replaceChild(textEl, fo);
  }});

  return clone;
}}

function dlSvg() {{
  const clone = _getSvgCloneWithFills();
  if (!clone) return;
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], {{type: 'image/svg+xml'}}));
  a.download = 'workflow_diagram.svg';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

function dlPng() {{
  const svg = document.querySelector('.mermaid svg');
  if (!svg) return;
  const clone = _getSvgCloneWithFills();
  if (!clone) return;
  const w = svg.viewBox.baseVal.width  || +svg.getAttribute('width')  || svg.getBoundingClientRect().width;
  const h = svg.viewBox.baseVal.height || +svg.getAttribute('height') || svg.getBoundingClientRect().height;
  clone.setAttribute('width',  w);
  clone.setAttribute('height', h);
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], {{type: 'image/svg+xml'}}));
  const img = new Image();
  img.onload = () => {{
    const canvas = document.createElement('canvas');
    canvas.width = w * 2; canvas.height = h * 2;
    const ctx = canvas.getContext('2d');
    ctx.scale(2, 2);
    ctx.fillStyle = '#fafaf8';
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0);
    URL.revokeObjectURL(url);
    canvas.toBlob(blob => {{
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'workflow_diagram.png';
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
    }});
  }};
  img.onerror = () => {{ URL.revokeObjectURL(url); }};
  img.src = url;
}}

function dlCode() {{
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([_MERMAID_DEF], {{type: 'text/plain'}}));
  a.download = 'workflow_diagram.mmd';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

// After Mermaid renders, expose SVG as id="diagram" for outer frame zoom
// and notify parent to resize the iframe
setTimeout(() => {{
  const svg = document.querySelector('.mermaid svg');
  if (svg) svg.id = 'diagram';
  try {{
    window.parent.postMessage({{
      type: 'snakeviz-rendered',
      height: document.documentElement.scrollHeight
    }}, '*');
  }} catch(e) {{}}
}}, 800);
</script>
</body>
</html>"""


# ════════════════════════════════════════════════════════════════════════════
# 8. PUBLIC API (called from JavaScript via Pyodide)
# ════════════════════════════════════════════════════════════════════════════

_last_node_count = 0
_last_edge_count = 0


def visualize_string(content: str, title: str = 'Workflow',
                     direction: str = 'LR') -> str:
    """
    Main entry point for browser / Pyodide use.
    Takes raw Snakefile text, returns a complete standalone HTML string.
    direction: 'LR' (left-to-right, default) or 'TB' (top-to-bottom).
    """
    global _last_node_count, _last_edge_count

    if not title:
        title = 'Workflow'
    if direction not in ('LR', 'TB'):
        direction = 'LR'

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

    return build_html(rules, edges, pos, title, groups, direction)


def visualize_mermaid_string(content: str, title: str = 'Workflow',
                             direction: str = 'LR') -> str:
    """
    Like visualize_string but returns a Mermaid.js / Excalidraw-styled HTML.
    """
    global _last_node_count, _last_edge_count

    if not title:
        title = 'Workflow'
    if direction not in ('LR', 'TB'):
        direction = 'LR'

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

    _last_node_count = len(rules)
    _last_edge_count = len(edges)

    mermaid_def = _build_mermaid_def(rules, edges, groups, order, direction)
    return _build_mermaid_html(mermaid_def, title)
