from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from xml.sax.saxutils import escape


OUT = Path(__file__).resolve().parents[1] / "EpiLink_method_schematic_editable.pptx"
EMU = 914400
SW, SH = 13.333, 7.5

NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
NS_P = "http://schemas.openxmlformats.org/presentationml/2006/main"
NS_R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

COL = {
    "navy": "173A69", "ink": "16243A", "gray": "9AA3AE", "lightgray": "E8ECF1",
    "pink": "C81B73", "blue": "2C6FD4", "teal": "5EB9D1", "purple": "6D3CA6",
    "green": "27A56A", "orange": "F39A46", "red": "D43B62", "lav": "EDE9F7",
    "pale_pink": "FBEAF4", "pale_blue": "E9F3FD", "pale_green": "EAF7F0",
}

sid = 1


def emu(v: float) -> int:
    return round(v * EMU)


def xel(tag: str, attrs: dict[str, str] | None = None, body: str = "") -> str:
    attrs = attrs or {}
    att = "".join(f' {k}="{escape(str(v))}"' for k, v in attrs.items())
    return f"<{tag}{att}>{body}</{tag}>"


def selfxel(tag: str, attrs: dict[str, str] | None = None) -> str:
    attrs = attrs or {}
    att = "".join(f' {k}="{escape(str(v))}"' for k, v in attrs.items())
    return f"<{tag}{att}/>"


def text_body(text: str, size=12, color="16243A", bold=False, align="l", font="Arial", valign="mid") -> str:
    paras = str(text).split("\n")
    ps = []
    for i, line in enumerate(paras):
        ppr = f'<a:pPr algn="{align}"><a:defRPr sz="{size*100}" b="{str(bold).lower()}"><a:solidFill><a:srgbClr val="{color}"/></a:solidFill><a:latin typeface="{font}"/></a:defRPr></a:pPr>'
        run = f'<a:r><a:rPr lang="en-US" sz="{size*100}" b="{str(bold).lower()}"><a:solidFill><a:srgbClr val="{color}"/></a:solidFill><a:latin typeface="{font}"/></a:rPr><a:t>{escape(line)}</a:t></a:r>'
        ps.append(f"<a:p>{ppr}{run}<a:endParaRPr lang=\"en-US\" sz=\"{size*100}\"/></a:p>")
    return f'<p:txBody><a:bodyPr anchor="{valign}"/><a:lstStyle/>{"".join(ps)}</p:txBody>'


def solid_fill(color: str, alpha: int | None = None) -> str:
    if color == "none":
        return "<a:noFill/>"
    if alpha is None:
        return f'<a:solidFill><a:srgbClr val="{color}"/></a:solidFill>'
    return f'<a:solidFill><a:srgbClr val="{color}"><a:alpha val="{alpha*1000}"/></a:srgbClr></a:solidFill>'


def line_props(color="16243A", width=1.0, dash=None, begin=None, end=None) -> str:
    body = f'<a:solidFill><a:srgbClr val="{color}"/></a:solidFill>'
    extra = ""
    if dash:
        extra += f'<a:prstDash val="{dash}"/>'
    if begin:
        extra += f'<a:headEnd type="{begin}" w="med" len="med"/>'
    if end:
        extra += f'<a:tailEnd type="{end}" w="med" len="med"/>'
    return f'<a:ln w="{emu(width/72)}">{body}{extra}</a:ln>'


def next_id() -> int:
    global sid
    sid += 1
    return sid


def nvsp(name: str, shape_id: int) -> str:
    return f'<p:nvSpPr><p:cNvPr id="{shape_id}" name="{escape(name)}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>'


def rect(x, y, w, h, fill="none", line="none", radius=False, text=None, size=12, color="16243A", bold=False,
         align="l", name="Rectangle", transparency=None):
    shape_id = next_id()
    geom = "roundRect" if radius else "rect"
    sppr = f'<p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm><a:prstGeom prst="{geom}"><a:avLst/></a:prstGeom>'
    if fill == "none":
        sppr += "<a:noFill/>"
    else:
        sppr += solid_fill(fill, transparency)
    if line == "none":
        sppr += "<a:ln><a:noFill/></a:ln>"
    else:
        sppr += line_props(line, 0.8)
    sppr += "</p:spPr>"
    tx = text_body(text, size, color, bold, align) if text is not None else ""
    return f'<p:sp>{nvsp(name, shape_id)}{sppr}{tx}</p:sp>'


def ellipse(x, y, w, h, fill="none", line="none", text=None, size=10, color="16243A", bold=False, name="Ellipse"):
    shape_id = next_id()
    sppr = f'<p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm><a:prstGeom prst="ellipse"><a:avLst/></a:prstGeom>{solid_fill(fill) if fill != "none" else "<a:noFill/>"}'
    sppr += f'{line_props(line, 0.7) if line != "none" else "<a:ln><a:noFill/></a:ln>"}</p:spPr>'
    tx = text_body(text, size, color, bold, "ctr") if text is not None else ""
    return f'<p:sp>{nvsp(name, shape_id)}{sppr}{tx}</p:sp>'


def line(x1, y1, x2, y2, color="16243A", width=1.0, end=None, begin=None, dash=None, name="Connector"):
    shape_id = next_id()
    sppr = f'<p:spPr><a:xfrm><a:off x="{emu(x1)}" y="{emu(y1)}"/><a:ext cx="{emu(x2-x1)}" cy="{emu(y2-y1)}"/></a:xfrm><a:prstGeom prst="line"><a:avLst/></a:prstGeom>{line_props(color, width, dash, begin, end)}</p:spPr>'
    return f'<p:cxnSp><p:nvCxnSpPr><p:cNvPr id="{shape_id}" name="{name}"/><p:cNvCxnSpPr/><p:nvPr/></p:nvCxnSpPr>{sppr}</p:cxnSp>'


def curve(x, y, w, h, color="16243A", width=1.5, end="triangle", name="Curved arrow"):
    shape_id = next_id()
    # Cubic Bezier in a 1000x1000 local coordinate system.
    geom = ('<a:custGeom><a:avLst/><a:gdLst/><a:ahLst/><a:cxnLst/>'
            '<a:pathLst><a:path w="1000" h="1000">'
            '<a:moveTo><a:pt x="0" y="620"/></a:moveTo>'
            '<a:cubicBezTo><a:pt x="300" y="80"/><a:pt x="730" y="80"/><a:pt x="1000" y="450"/></a:cubicBezTo>'
            '</a:path></a:pathLst></a:custGeom>')
    sppr = f'<p:spPr><a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm>{geom}{line_props(color, width, None, None, end)}<a:noFill/></p:spPr>'
    return f'<p:sp>{nvsp(name, shape_id)}{sppr}</p:sp>'


def text_box(x, y, w, h, text, size=12, color="16243A", bold=False, align="l", name="Text"):
    return rect(x, y, w, h, fill="none", line="none", text=text, size=size, color=color, bold=bold, align=align, name=name)


def matrix(x, y, cols, rows, cw, ch, colors):
    out = []
    for r in range(rows):
        for c in range(cols):
            col = colors[(r*3+c*2) % len(colors)]
            out.append(rect(x+c*cw, y+r*ch, cw-0.015, ch-0.015, fill=col, line="FFFFFF", name="Matrix cell"))
    return "".join(out)


def dots(x, y, colors, counts, dx=0.08, dy=0.07):
    out = []
    for i, (color, count) in enumerate(zip(colors, counts)):
        for j in range(count):
            xx = x + (j % 7) * dx + ((i*3) % 3) * 0.02
            yy = y + i * 0.43 + (j // 7) * dy + ((j*5+i) % 4) * 0.01
            out.append(ellipse(xx, yy, 0.045, 0.045, fill=color, line="none", name="Cell dot"))
    return "".join(out)


def bar_track(x, y, heights, color, scale=0.015):
    out = [line(x, y+0.24, x+2.20, y+0.24, color=color, width=0.7)]
    for i, v in enumerate(heights):
        out.append(rect(x+0.05+i*0.13, y+0.24-v*scale, 0.07, v*scale, fill=color, line="none", name="Accessibility peak"))
    return "".join(out)


def plot_curve(x, y, w, h, color, dashed=False, label=None):
    out = [line(x, y+h, x+w, y+h, color="77818E", width=0.7), line(x, y, x, y+h, color="77818E", width=0.7)]
    pts = [(0.00, .90), (.08, .72), (.18, .48), (.30, .30), (.44, .22), (.58, .17), (.74, .13), (1.0, .10)]
    for (a, b), (c, d) in zip(pts, pts[1:]):
        out.append(line(x+a*w, y+b*h, x+c*w, y+d*h, color=color, width=1.5, dash="dash" if dashed else None))
    if label:
        out.append(text_box(x+w*0.52, y+h*0.44, w*0.45, 0.22, label, size=8, color=color, name="Plot label"))
    return "".join(out)


def slide_xml() -> str:
    parts = []
    # Background
    parts.append(rect(0, 0, SW, SH, fill="FFFFFF", line="none", name="White background"))
    # Panel a header
    parts.append(text_box(0.16, 0.08, 0.40, 0.36, "a", size=24, color=COL["ink"], bold=True, name="Panel label a"))
    parts.append(text_box(0.70, 0.08, 1.65, 0.36, "EpiLink", size=25, color=COL["navy"], bold=True, name="Title"))
    parts.append(line(2.44, 0.13, 2.44, 0.40, color=COL["gray"], width=0.8))
    parts.append(text_box(2.60, 0.10, 6.0, 0.30, "Cell-type-specific enhancer–gene mapping from scMultiome", size=15, color=COL["ink"], bold=True, name="Panel a title"))

    # scMultiome cluster map
    parts.append(text_box(0.78, 0.50, 1.35, 0.24, "scMultiome", size=12, color=COL["ink"], bold=True, name="scMultiome label"))
    parts.append(line(0.45, 1.83, 0.45, 0.88, color=COL["ink"], width=0.9, end="triangle"))
    parts.append(line(0.45, 1.83, 1.87, 1.83, color=COL["ink"], width=0.9, end="triangle"))
    parts.append(text_box(0.22, 0.83, 0.33, 0.20, "UMAP2", size=7, color=COL["gray"], name="UMAP2"))
    parts.append(text_box(1.43, 1.88, 0.42, 0.20, "UMAP1", size=7, color=COL["gray"], name="UMAP1"))
    parts.append(dots(0.70, 0.90, [COL["green"], COL["pink"], COL["blue"], COL["orange"], COL["purple"], COL["teal"]], [22, 20, 20, 19, 20, 18], dx=0.08, dy=0.06))
    parts.append(text_box(1.68, 0.98, 0.85, 0.20, "Cell type X", size=8, color=COL["pink"], bold=True, name="Cell type X label"))
    parts.append(text_box(1.88, 1.36, 0.85, 0.20, "Cell type Y", size=8, color=COL["blue"], bold=True, name="Cell type Y label"))
    parts.append(line(2.15, 1.22, 2.55, 1.22, color=COL["gray"], width=1.2, end="triangle"))

    # Cell-type profile cards
    parts.append(rect(2.68, 0.56, 2.44, 0.86, fill=COL["pale_pink"], line="none", radius=True, name="Cell type X profile card"))
    parts.append(text_box(2.82, 0.65, 0.95, 0.20, "Cell type X", size=10, color=COL["pink"], bold=True, name="Cell type X"))
    parts.append(text_box(2.82, 0.91, 0.40, 0.18, "ATAC", size=8, color=COL["ink"], name="ATAC X"))
    parts.append(bar_track(3.32, 0.83, [3, 10, 4, 18, 35, 7, 12, 4, 28, 6, 8, 25, 4], COL["pink"], scale=0.012))
    parts.append(text_box(2.82, 1.16, 0.40, 0.18, "RNA", size=8, color=COL["ink"], name="RNA X"))
    parts.append(bar_track(3.32, 1.08, [2, 4, 8, 3, 19, 6, 4, 24, 4], COL["pink"], scale=0.009))
    parts.append(rect(2.68, 1.51, 2.44, 0.86, fill=COL["pale_blue"], line="none", radius=True, name="Cell type Y profile card"))
    parts.append(text_box(2.82, 1.60, 0.95, 0.20, "Cell type Y", size=10, color=COL["blue"], bold=True, name="Cell type Y"))
    parts.append(text_box(2.82, 1.86, 0.40, 0.18, "ATAC", size=8, color=COL["ink"], name="ATAC Y"))
    parts.append(bar_track(3.32, 1.78, [2, 24, 6, 4, 14, 4, 5, 19, 2, 15], COL["blue"], scale=0.012))
    parts.append(text_box(2.82, 2.11, 0.40, 0.18, "RNA", size=8, color=COL["ink"], name="RNA Y"))
    parts.append(bar_track(3.32, 2.03, [2, 3, 4, 5, 3, 20, 7, 3, 27], COL["blue"], scale=0.009))

    # Genomic locus and arcs
    parts.append(line(5.38, 1.65, 12.95, 1.65, color=COL["ink"], width=1.2, name="Genomic locus"))
    for xx, ww, cc, lab in [(5.78, .50, COL["pink"], "Enhancer"), (6.95, .50, COL["blue"], "Enhancer"), (8.42, .50, COL["gray"], "candidate"), (10.72, .52, COL["green"], "Promoter")]:
        parts.append(rect(xx, 1.47, ww, 0.36, fill=cc, line=COL["ink"], name=lab))
        if lab != "candidate":
            parts.append(text_box(xx-0.05, 1.88, ww+0.10, 0.18, lab, size=7, color=COL["ink"], align="ctr", name=lab+" label"))
    parts.append(rect(11.32, 1.46, 1.35, 0.38, fill=COL["purple"], line=COL["ink"], name="Gene A"))
    parts.append(text_box(11.40, 1.90, 1.20, 0.18, "Gene A", size=9, color=COL["purple"], bold=True, align="ctr", name="Gene A label"))
    parts.append(line(10.98, 1.42, 10.98, 1.18, color=COL["green"], width=1.0, end="triangle", name="Promoter direction"))
    parts.append(curve(5.80, 0.77, 5.05, 0.75, color=COL["pink"], width=1.8, end="triangle", name="Cell type X link"))
    parts.append(curve(6.72, 0.96, 4.42, 0.56, color=COL["blue"], width=1.8, end="triangle", name="Cell type Y link"))
    parts.append(text_box(10.90, 0.60, 0.86, 0.18, "X link", size=8, color=COL["pink"], bold=True, name="X link legend"))
    parts.append(line(10.63, 0.70, 10.88, 0.70, color=COL["pink"], width=1.7))
    parts.append(text_box(10.90, 0.85, 0.86, 0.18, "Y link", size=8, color=COL["blue"], bold=True, name="Y link legend"))
    parts.append(line(10.63, 0.95, 10.88, 0.95, color=COL["blue"], width=1.7))
    parts.append(text_box(7.60, 2.23, 2.75, 0.28, "cell-type-specific enhancer–gene mapping", size=9, color=COL["navy"], bold=True, align="ctr", name="Mapping annotation"))

    # Divider
    parts.append(line(0.12, 2.60, 13.18, 2.60, color=COL["lightgray"], width=1.3, name="Panel divider"))

    # Panel b header
    parts.append(text_box(0.16, 2.70, 0.40, 0.36, "b", size=24, color=COL["ink"], bold=True, name="Panel label b"))
    parts.append(text_box(0.70, 2.70, 1.65, 0.36, "EpiLink", size=25, color=COL["navy"], bold=True, name="Title b"))
    parts.append(line(2.44, 2.75, 2.44, 3.02, color=COL["gray"], width=0.8))
    parts.append(text_box(2.60, 2.72, 6.2, 0.30, "Nullstrap-based inference with empirical FDR calibration", size=15, color=COL["ink"], bold=True, name="Panel b title"))

    # Workflow cards
    card_y, card_h = 3.22, 3.05
    cards = [(0.10, 2.35, "Observed scMultiome data"), (2.60, 1.56, "Poisson/NB-GLM"), (4.30, 1.62, "Nullstrap null data"), (6.05, 1.47, "Refit same GLM"), (7.65, 1.70, "Compare test statistics"), (9.50, 1.58, "empirical FDR\ncalibration"), (11.22, 1.98, "EpiLink links")]
    for x, w, title in cards:
        parts.append(rect(x, card_y, w, card_h, fill="FFFFFF", line=COL["lightgray"], radius=True, name=title.replace("\n", " ")+" card"))
        parts.append(text_box(x+0.10, card_y+0.12, w-0.20, 0.30 if "\n" not in title else 0.45, title, size=9, color=COL["ink"], bold=True, align="ctr", name=title.replace("\n", " ")+" title"))

    # Card 1
    parts.append(matrix(0.58, 4.05, 4, 4, 0.27, 0.27, ["B8D7F5", "7EB3E6", "DCCEF0", "B89ED8", "C9EAD9", "F2C6DD"]))
    parts.append(text_box(0.66, 5.24, 1.25, 0.20, "Cells × genomic features", size=7, color=COL["gray"], align="ctr", name="Matrix caption"))
    parts.append(text_box(0.64, 3.73, 0.50, 0.18, "ATAC", size=8, color=COL["ink"], align="ctr", name="ATAC matrix label"))
    parts.append(text_box(1.25, 3.73, 0.50, 0.18, "RNA", size=8, color=COL["ink"], align="ctr", name="RNA matrix label"))
    for i, cc in enumerate([COL["pink"], COL["blue"], COL["orange"], COL["green"]]):
        parts.append(ellipse(0.27, 4.05+i*0.29, 0.11, 0.11, fill=cc, line="none", name="Cell type indicator"))

    # Card 2: GLM plot
    parts.append(line(2.91, 5.20, 3.92, 5.20, color=COL["ink"], width=0.8, end="triangle"))
    parts.append(line(2.91, 5.20, 2.91, 4.08, color=COL["ink"], width=0.8, end="triangle"))
    for i, (xx, yy) in enumerate([(3.03,5.02),(3.16,4.83),(3.28,4.91),(3.42,4.65),(3.54,4.57),(3.68,4.37),(3.80,4.23)]):
        parts.append(ellipse(xx, yy, .08, .08, fill=COL["blue"], line="none", name="GLM observation"))
    for a,b,c,d in [(2.99,5.02,3.28,4.84),(3.28,4.84,3.58,4.51),(3.58,4.51,3.90,4.18)]:
        parts.append(line(a,b,c,d,color=COL["ink"],width=1.3))
    parts.append(text_box(3.03, 5.31, 0.84, 0.18, "ATAC", size=7, color=COL["ink"], align="ctr", name="ATAC axis"))
    parts.append(text_box(2.64, 4.26, 0.30, 0.18, "RNA", size=7, color=COL["ink"], name="RNA axis"))

    # Card 3 null matrix and crossed link
    parts.append(matrix(4.62, 4.04, 4, 4, 0.24, 0.24, ["C8DDF6", "B6D1EE", "E3D8F2", "D6C7E8", "D2EFE0"]))
    parts.append(matrix(5.24, 4.04, 4, 4, 0.24, 0.24, ["E7D7F0", "DBC9E9", "F0DCE7", "E3CBE1", "EFE7D0"]))
    parts.append(line(4.70, 5.20, 5.89, 5.20, color=COL["pink"], width=1.3, dash="dash"))
    parts.append(line(5.10, 4.86, 5.54, 5.34, color=COL["pink"], width=1.7, begin="triangle", end="triangle"))
    parts.append(text_box(4.54, 5.39, 1.08, 0.30, "association removed", size=7, color=COL["pink"], align="ctr", name="Association removed"))

    # Card 4 refit plot
    parts.append(line(6.34, 5.20, 7.30, 5.20, color=COL["ink"], width=0.8, end="triangle"))
    parts.append(line(6.34, 5.20, 6.34, 4.08, color=COL["ink"], width=0.8, end="triangle"))
    for i, (xx, yy) in enumerate([(6.48,5.03),(6.58,4.91),(6.70,4.76),(6.82,4.68),(6.95,4.53),(7.08,4.38),(7.19,4.30)]):
        parts.append(ellipse(xx, yy, .07, .07, fill=COL["gray"], line="none", name="Null GLM observation"))
    for a,b,c,d in [(6.46,5.05,6.75,4.80),(6.75,4.80,7.04,4.55),(7.04,4.55,7.28,4.32)]:
        parts.append(line(a,b,c,d,color=COL["gray"],width=1.2,dash="dash"))
    parts.append(text_box(6.48, 5.31, 0.74, 0.18, "ATAC", size=7, color=COL["ink"], align="ctr", name="Refit ATAC"))
    parts.append(text_box(6.08, 4.27, 0.30, 0.18, "RNA", size=7, color=COL["ink"], name="Refit RNA"))

    # Card 5 distributions
    parts.append(line(7.90, 5.20, 9.05, 5.20, color=COL["ink"], width=0.8, end="triangle"))
    parts.append(line(7.90, 5.20, 7.90, 4.08, color=COL["ink"], width=0.8, end="triangle"))
    parts.append(curve(7.98, 4.17, 0.75, 0.87, color=COL["blue"], width=1.5, end=None, name="Observed distribution"))
    parts.append(curve(8.37, 4.34, 0.76, 0.70, color=COL["gray"], width=1.5, end=None, name="Null distribution"))
    parts.append(line(8.55, 4.08, 8.55, 5.20, color=COL["ink"], width=1.0, dash="dash"))
    parts.append(text_box(8.03, 5.31, 0.95, 0.18, "Test statistic", size=7, color=COL["ink"], align="ctr", name="Test statistic axis"))
    parts.append(text_box(8.65, 4.16, 0.55, 0.18, "Observed", size=7, color=COL["blue"], name="Observed legend"))
    parts.append(text_box(8.65, 4.38, 0.42, 0.18, "Null", size=7, color=COL["gray"], name="Null legend"))

    # Card 6 FDR curve
    parts.append(line(9.82, 5.20, 10.85, 5.20, color=COL["ink"], width=0.8, end="triangle"))
    parts.append(line(9.82, 5.20, 9.82, 4.08, color=COL["ink"], width=0.8, end="triangle"))
    parts.append(curve(9.90, 4.14, 0.88, 0.90, color=COL["green"], width=1.7, end=None, name="Empirical FDR curve"))
    parts.append(line(9.88, 4.64, 10.78, 4.64, color=COL["green"], width=0.9, dash="dash"))
    parts.append(text_box(9.96, 4.68, 0.20, 0.18, "q", size=8, color=COL["green"], bold=True, name="FDR q"))
    parts.append(text_box(10.02, 5.31, 0.64, 0.18, "Threshold", size=7, color=COL["ink"], align="ctr", name="Threshold axis"))
    parts.append(text_box(9.56, 4.20, 0.20, 0.40, "FDR", size=7, color=COL["ink"], name="FDR axis"))

    # Card 7 output link
    parts.append(line(11.54, 5.20, 12.92, 5.20, color=COL["ink"], width=1.0))
    parts.append(rect(11.68, 4.84, 0.38, 0.28, fill=COL["pink"], line=COL["ink"], name="Selected enhancer"))
    parts.append(rect(12.17, 4.84, 0.38, 0.28, fill=COL["green"], line=COL["ink"], name="Selected promoter"))
    parts.append(rect(12.54, 4.84, 0.85, 0.29, fill=COL["purple"], line=COL["ink"], name="Selected gene"))
    parts.append(curve(11.86, 4.25, 1.05, 0.48, color=COL["pink"], width=1.5, end="triangle", name="Selected EpiLink"))
    parts.append(text_box(11.62, 5.28, 0.54, 0.18, "Enhancer", size=7, color=COL["ink"], align="ctr", name="Output enhancer label"))
    parts.append(text_box(12.10, 5.28, 0.64, 0.18, "Promoter", size=7, color=COL["ink"], align="ctr", name="Output promoter label"))
    parts.append(text_box(12.66, 5.28, 0.58, 0.18, "Gene A", size=7, color=COL["purple"], bold=True, align="ctr", name="Output gene label"))
    parts.append(rect(11.45, 5.76, 1.55, 0.35, fill=COL["pale_green"], line="none", radius=True, name="Robustness badge"))
    parts.append(ellipse(11.58, 5.84, 0.18, 0.18, fill=COL["green"], line="none", text="✓", size=10, color="FFFFFF", bold=True, name="Check mark"))
    parts.append(text_box(11.83, 5.79, 1.05, 0.22, "robust to model\nmisspecification", size=7, color=COL["ink"], bold=True, name="Robustness text"))
    # Flow arrows between cards
    for xx in [2.48, 4.18, 5.93, 7.53, 9.38, 11.10]:
        parts.append(line(xx, 4.74, xx+0.12, 4.74, color=COL["ink"], width=1.1, end="triangle", name="Workflow arrow"))

    # Footer note
    parts.append(text_box(0.28, 6.92, 12.70, 0.25, "Concept reconstruction: native PowerPoint shapes and text; refine exact model notation and link-generation details as needed.", size=7, color=COL["gray"], align="ctr", name="Footer note"))

    sp_tree = '<p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name="Group 1"/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>' + "".join(parts) + "</p:spTree>"
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld xmlns:a="{NS_A}" xmlns:r="{NS_R}" xmlns:p="{NS_P}"><p:cSld name="EpiLink method schematic"><p:spTree>{sp_tree.split('<p:spTree>',1)[1].rsplit('</p:spTree>',1)[0]}</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>'''


def package_pptx(path: Path):
    slide = slide_xml()
    content_types = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>
<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>
<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>
<Override PartName="/ppt/slides/slide1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>
<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>'''
    root_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="ppt/presentation.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>'''
    pres = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentation xmlns:a="{NS_A}" xmlns:r="{NS_R}" xmlns:p="{NS_P}" saveSubsetFonts="1">
<p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst>
<p:sldIdLst><p:sldId id="256" r:id="rId2"/></p:sldIdLst>
<p:sldSz cx="{emu(SW)}" cy="{emu(SH)}" type="screen16x9"/><p:notesSz cx="6858000" cy="9144000"/>
</p:presentation>'''
    pres_rels = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="slideMasters/slideMaster1.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide1.xml"/>
</Relationships>'''
    slide_rels = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/></Relationships>'''
    master_rels = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme" Target="../theme/theme1.xml"/></Relationships>'''
    layout_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="../slideMasters/slideMaster1.xml"/></Relationships>'''
    common_inner = f'''<p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>'''
    master = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldMaster xmlns:a="{NS_A}" xmlns:r="{NS_R}" xmlns:p="{NS_P}"><p:cSld name="Master"><p:spTree>{common_inner}</p:spTree></p:cSld><p:sldLayoutIdLst><p:sldLayoutId id="1" r:id="rId1"/></p:sldLayoutIdLst><p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/></p:sldMaster>'''
    layout = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldLayout xmlns:a="{NS_A}" xmlns:r="{NS_R}" xmlns:p="{NS_P}" type="blank" preserve="1"><p:cSld name="Blank"><p:spTree>{common_inner}</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>'''
    theme = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="{NS_A}" name="EpiLink Theme"><a:themeElements><a:clrScheme name="EpiLink"><a:dk1><a:srgbClr val="000000"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1><a:dk2><a:srgbClr val="16243A"/></a:dk2><a:lt2><a:srgbClr val="FFFFFF"/></a:lt2><a:accent1><a:srgbClr val="173A69"/></a:accent1><a:accent2><a:srgbClr val="C81B73"/></a:accent2><a:accent3><a:srgbClr val="2C6FD4"/></a:accent3><a:accent4><a:srgbClr val="27A56A"/></a:accent4><a:accent5><a:srgbClr val="F39A46"/></a:accent5><a:accent6><a:srgbClr val="6D3CA6"/></a:accent6><a:hlink><a:srgbClr val="0563C1"/></a:hlink><a:folHlink><a:srgbClr val="954F72"/></a:folHlink></a:clrScheme><a:fontScheme name="EpiLink"><a:majorFont><a:latin typeface="Arial"/></a:majorFont><a:minorFont><a:latin typeface="Arial"/></a:minorFont></a:fontScheme><a:fmtScheme name="EpiLink"><a:fillStyleLst/><a:lnStyleLst/><a:effectStyleLst/><a:bgFillStyleLst/></a:fmtScheme></a:themeElements></a:theme>'''
    core = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:title>EpiLink method schematic</dc:title><dc:subject>Editable reconstruction of EpiLink methods overview</dc:subject><dc:creator>OpenAI</dc:creator></cp:coreProperties>'''
    app = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"><Application>PowerPoint</Application><PresentationFormat>On-screen Show (16:9)</PresentationFormat><Slides>1</Slides></Properties>'''
    with ZipFile(path, "w", ZIP_DEFLATED) as z:
        files = {
            "[Content_Types].xml": content_types,
            "_rels/.rels": root_rels,
            "ppt/presentation.xml": pres,
            "ppt/_rels/presentation.xml.rels": pres_rels,
            "ppt/slides/slide1.xml": slide,
            "ppt/slides/_rels/slide1.xml.rels": slide_rels,
            "ppt/slideMasters/slideMaster1.xml": master,
            "ppt/slideMasters/_rels/slideMaster1.xml.rels": master_rels,
            "ppt/slideLayouts/slideLayout1.xml": layout,
            "ppt/slideLayouts/_rels/slideLayout1.xml.rels": layout_rels,
            "ppt/theme/theme1.xml": theme,
            "docProps/core.xml": core,
            "docProps/app.xml": app,
        }
        for name, data in files.items():
            z.writestr(name, data)


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    package_pptx(OUT)
    print(OUT)
