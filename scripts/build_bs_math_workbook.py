"""Build eight printable Chinese lessons and a bookmarked combined workbook.

Source of truth: docs/bs_math/*.md. Equations are rendered by Matplotlib at
400 dpi; Chinese text is embedded using a local TrueType font. Production
pricing code is not imported or changed. PDFs remain generated local outputs.
"""
from pathlib import Path
import argparse
import hashlib
import html
import os
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "docs" / "bs_math"
OUTPUT = ROOT / "output" / "pdf" / "black_scholes_math"
TEMP = ROOT / "out" / "bs-math-qa"

# Use the desktop's bundled ReportLab when absent from the active interpreter.
try:
    import reportlab
except ImportError:
    runtime = Path(os.environ.get("USERPROFILE", str(Path.home()))) / ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/Lib/site-packages"
    if not runtime.is_dir():
        raise SystemExit("ReportLab is missing. Install reportlab in the active Python environment.")
    sys.path.append(str(runtime))

from PIL import Image as PILImage
import matplotlib
matplotlib.use("Agg")
from matplotlib.mathtext import math_to_image
from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                              Table, TableStyle, HRFlowable, KeepTogether)

PAGE_W, PAGE_H = A4
WIDTH = PAGE_W - 100
INK = colors.HexColor("#183249")
MUTED = colors.HexColor("#547080")
ACCENT = colors.HexColor("#276d83")


def register_fonts():
    selected = os.environ.get("BS_WORKBOOK_FONT")
    candidates = [Path(selected)] if selected else []
    candidates += [Path("C:/Windows/Fonts/simsun.ttc"),
                   Path("/usr/share/fonts/truetype/arphic/uming.ttc")]
    for path in candidates:
        if path.is_file():
            pdfmetrics.registerFont(TTFont("Chinese",str(path), subfontIndex=0))
            pdfmetrics.registerFontFamily("Chinese",normal="Chinese",bold="Chinese",italic="Chinese",boldItalic="Chinese")
            return
    raise SystemExit("Missing Chinese TrueType font. Set BS_WORKBOOK_FONT to a compatible font path.")


def styles():
    # ReportLab's CJK splitter rejects empty text fragments adjacent to inline
    # equation images. Ordinary wrapping supports long Chinese runs via character
    # splitting and handles those image fragments correctly.
    base = dict(fontName="Chinese", textColor=INK, wordWrap="LTR", splitLongWords=1, autoLeading="max")
    return {
        "body": ParagraphStyle("body",fontSize=10.4,leading=17,spaceAfter=7,**base),
        "h1": ParagraphStyle("h1",fontSize=21,leading=30,spaceBefore=12,spaceAfter=18,**base),
        "h2": ParagraphStyle("h2",fontSize=12.5,leading=21,spaceBefore=15,spaceAfter=8,keepWithNext=True,**base),
        "h3": ParagraphStyle("h3",fontSize=11.3,leading=18,spaceBefore=10,spaceAfter=6,keepWithNext=True,**base),
        "caption": ParagraphStyle("caption",fontSize=8.8,leading=13,spaceAfter=9,**base),
        "cell": ParagraphStyle("cell",fontSize=9,leading=14,spaceAfter=0,**base),
        "code": ParagraphStyle("code",fontSize=8.6,leading=14,spaceAfter=6,**base),
        "meta": ParagraphStyle("meta",fontSize=9,leading=14,spaceAfter=12,textColor=MUTED,fontName="Chinese",wordWrap="CJK"),
    }


def equation(latex, display=False, max_width=WIDTH):
    # mathtext provides dependable math rendering without a terminal TeX install.
    latex = latex.replace(r"\!", "").replace(r"\,", " ")
    latex = re.sub(r"\\frac\s*([0-9a-zA-Z])\s*([0-9a-zA-Z])", r"\\frac{\1}{\2}", latex)
    latex = re.sub(r"\\frac(\{[^{}]*\})\s*([0-9a-zA-Z])", r"\\frac\1{\2}", latex)
    latex = re.sub(r"\\frac\s*([0-9a-zA-Z])(?=\{)", r"\\frac{\1}", latex)
    latex = re.sub(r"\\sqrt\s*([0-9a-zA-Z])", r"\\sqrt{\1}", latex)
    latex = re.sub(r"\\le(?![a-zA-Z])", r"\\leq", latex)
    latex = re.sub(r"\\ge(?![a-zA-Z])", r"\\geq", latex)
    digest = hashlib.sha256(latex.encode("utf-8")).hexdigest()[:20]
    path = TEMP / "equations" / f"{digest}.png"
    if not path.exists():
        math_to_image(f"${latex}$", str(path), prop=matplotlib.font_manager.FontProperties(size=12),
                      dpi=400, format="png",color="#183249")
    with PILImage.open(path) as img:
        pxw,pxh=img.size
    # In points, proportional to the original 12 pt equation.
    w,h=pxw*72/400,pxh*72/400
    scale=(1.0 if display else .82)
    if w*scale>max_width:
        scale=max_width/w
    return path,w*scale,h*scale


def rich(text, width=WIDTH):
    # Convert source Markdown into ReportLab's limited paragraph markup.
    pieces=[]
    offset=0
    for match in re.finditer(r"\$([^$\n]+)\$",text):
        pieces.append(html.escape(text[offset:match.start()]))
        path,w,h=equation(match.group(1), max_width=width-8)
        pieces.append(f'<img src="{html.escape(str(path),quote=True)}" width="{w:.2f}" height="{h:.2f}" valign="-2"/>')
        offset=match.end()
    pieces.append(html.escape(text[offset:]))
    result="".join(pieces)
    def link(match):
        label,url=match.groups()
        if url.startswith(("https://","http://")):
            return f'<link href="{html.escape(html.unescape(url),quote=True)}" color="#276d83">{label}</link>'
        return label
    result=re.sub(r"\[([^\]]+)\]\(([^)]+)\)",link,result)
    result=re.sub(r"\*\*([^*]+)\*\*",r"<b>\1</b>",result)
    result=re.sub(r"`([^`]+)`",r"\1",result)
    return result


def paper_pause():
    return [Spacer(1,3),Paragraph("纸笔区：先在这里写出中间步骤，再阅读答案。",STYLES["caption"])] + sum(
        ([HRFlowable(width=WIDTH,color=colors.HexColor("#d8e1e7"),thickness=.4),Spacer(1,21)] for _ in range(3)), [])


def table(lines):
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in lines]
    rows = [row for row in rows if not all(re.fullmatch(r"[:\- ]+",cell) for cell in row)]
    n=len(rows[0])
    # Prefer wider explanatory columns; keep the original order intact.
    if n==3:
        colwidths=[WIDTH*.24,WIDTH*.36,WIDTH*.40]
    elif n==4:
        colwidths=[WIDTH*.15,WIDTH*.24,WIDTH*.24,WIDTH*.37]
    elif n==5:
        colwidths=[WIDTH*.12]+[WIDTH*.22]*4
    else:
        colwidths=[WIDTH/n]*n
    cells=[[Paragraph(rich(cell,colwidths[i]-12),STYLES["cell"]) for i,cell in enumerate(row)] for row in rows]
    t=Table(cells,colWidths=colwidths,repeatRows=1,hAlign="LEFT",splitByRow=int(len(rows)>10))
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#e8f0f4")),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f7fafb")]),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),6),("RIGHTPADDING",(0,0),(-1,-1),6),
        ("TOPPADDING",(0,0),(-1,-1),7),("BOTTOMPADDING",(0,0),(-1,-1),7),
        ("LINEBELOW",(0,0),(-1,0),.7,ACCENT),
        ("LINEBELOW",(0,1),(-1,-1),.3,colors.HexColor("#d8e1e7")),
    ]))
    return [t,Spacer(1,9)]


def source_to_story(path):
    lines=path.read_text(encoding="utf-8").splitlines()
    story=[]; i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:
            i+=1; continue
        if line=="$$":
            i+=1; formula=[]
            while i<len(lines) and lines[i].strip()!="$$":
                formula.append(lines[i].strip()); i+=1
            if i==len(lines):
                raise ValueError(f"Unclosed display formula in {path}")
            eq,w,h=equation(" ".join(formula),display=True)
            story.append(KeepTogether([Spacer(1,5),Image(str(eq),width=w,height=h,hAlign="CENTER"),Spacer(1,11)]))
            i+=1; continue
        if line.startswith("|"):
            chunk=[]
            while i<len(lines) and lines[i].strip().startswith("|"):
                chunk.append(lines[i]); i+=1
            story.extend(table(chunk)); continue
        image=re.fullmatch(r"!\[([^\]]*)\]\(([^)]+)\)",line)
        if image:
            imgpath=(path.parent/image.group(2)).resolve()
            with PILImage.open(imgpath) as img:
                w,h=img.size
            wpt=WIDTH; hpt=h/w*wpt
            if hpt>330:
                wpt*=330/hpt; hpt=330
            story.append(KeepTogether([Image(str(imgpath),width=wpt,height=hpt),Spacer(1,4),Paragraph(image.group(1),STYLES["caption"])]))
            i+=1; continue
        if line.startswith("```"):
            i+=1; chunks=[]
            while i<len(lines) and not lines[i].strip().startswith("```"):
                chunks.append(lines[i]); i+=1
            # Wrap a machine-readable CSV header visually after comma tokens.
            for raw in chunks:
                parts=[]
                if len(raw)>85 and "," in raw:
                    current=""
                    for token in raw.split(","):
                        addition=token+"," if token!=raw.split(",")[-1] else token
                        if len(current)+len(addition)>75:
                            parts.append(current); current=""
                        current+=addition
                    parts.append(current)
                else:
                    parts=[raw]
                story.extend(Paragraph(html.escape(part),STYLES["code"]) for part in parts)
            i+=1; continue
        heading=re.match(r"^(#{1,3}) (.*)$",line)
        if heading:
            kind=f"h{len(heading.group(1))}"
            story.append(Paragraph(rich(heading.group(2)),STYLES[kind]))
            if kind=="h1":
                story.append(Paragraph("Black-Scholes 数学学习资料 | 先推导，再对照数值",STYLES["meta"]))
            i+=1; continue
        chunk=[line]; i+=1
        while i<len(lines) and lines[i].strip() and not re.match(r"^(#|\$\$|\||```|!\[)",lines[i].strip()):
            if re.match(r"^(\d+\.|-) ",lines[i].strip()):
                break
            chunk.append(lines[i].strip()); i+=1
        paragraph=" ".join(chunk)
        paragraph_flow=Paragraph(rich(paragraph),STYLES["body"])
        following = next((line.strip() for line in lines[i:] if line.strip()), "")
        if following == "$$":
            paragraph_flow.keepWithNext = True
        if paragraph.startswith("**纸笔停顿"):
            story.append(KeepTogether([paragraph_flow]+paper_pause()))
        else:
            story.append(paragraph_flow)
    # Keep clickable source labels in each lesson and print the complete source
    # addresses in the course index, avoiding a nearly empty references-only page.
    links=list(dict.fromkeys(re.findall(r"\[[^\]]*\]\((https?://[^)]+)\)",path.read_text(encoding="utf-8"))))
    if links and path.name=="README.md":
        addresses="<br/>".join(f'<link href="{html.escape(link,quote=True)}" color="#276d83">{html.escape(link)}</link>' for link in links)
        story.append(Spacer(1,7))
        story.append(Paragraph("本篇在线出处<br/>"+addresses,STYLES["caption"]))
    return story


def footer(title):
    def draw(canvas,doc):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#d8e1e7")); canvas.line(50, PAGE_H-36,PAGE_W-50,PAGE_H-36)
        canvas.setFont("Chinese",8); canvas.setFillColor(MUTED)
        canvas.drawString(50,PAGE_H-29,title[:42])
        canvas.drawString(50,29,"欧式期权 | 连续股息 | Theta = 日历时间导数")
        canvas.drawRightString(PAGE_W-50,29,str(doc.page))
        canvas.restoreState()
    return draw


def build_pdf(source,target):
    title=source.read_text(encoding="utf-8").splitlines()[0].removeprefix("# ")
    doc=SimpleDocTemplate(str(target),pagesize=A4,leftMargin=50,rightMargin=50,
                          topMargin=51,bottomMargin=48,title=title,
                          author="Quant_interview 数学学习资料",pageCompression=1)
    doc.build(source_to_story(source),onFirstPage=footer(title),onLaterPages=footer(title))
    return title,len(PdfReader(target).pages)


def main():
    global STYLES
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lessons", nargs="+", choices=[f"{i:02}" for i in range(1,9)],
                        help="Rebuild selected lessons and the combined workbook; other lesson PDFs must already exist.")
    parser.add_argument("--standalone", action="store_true",
                        help="Build lesson PDFs without regenerating the combined workbook.")
    args=parser.parse_args()
    sources=sorted(SOURCES.glob("[0-9][0-9]_*.md"))
    selected=set(args.lessons) if args.lessons else {source.stem[:2] for source in sources}
    missing=[OUTPUT/f"{source.stem}.pdf" for source in sources
             if source.stem[:2] not in selected and not (OUTPUT/f"{source.stem}.pdf").is_file()]
    if missing and not args.standalone:
        parser.error("Unbuilt lesson PDFs: "+", ".join(path.name for path in missing)+". Rebuild all lessons first.")
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (TEMP/"equations").mkdir(parents=True,exist_ok=True)
    register_fonts(); STYLES=styles()
    if args.standalone:
        for source in sources:
            if source.stem[:2] in selected:
                target=OUTPUT/f"{source.stem}.pdf"
                _, pages=build_pdf(source,target)
                print(f"{target.name}: {pages} pages")
        return
    items=[]
    index=TEMP/"course_index.pdf"
    title,pages=build_pdf(SOURCES/"README.md",index)
    items.append((index,title,pages))
    for source in sources:
        target=OUTPUT/f"{source.stem}.pdf"
        if source.stem[:2] in selected:
            title,pages=build_pdf(source,target)
            print(f"{target.name}: {pages} pages (rebuilt)")
        else:
            title=source.read_text(encoding="utf-8").splitlines()[0].removeprefix("# ")
            pages=len(PdfReader(target).pages)
        items.append((target,title,pages))
    writer=PdfWriter()
    for path,title,_ in items:
        writer.append(str(path),outline_item=title)
    writer.add_metadata({"/Title":"Black-Scholes 与 Greeks：纸笔推导合订本", "/Author":"Quant_interview"})
    combined=OUTPUT/"black_scholes_paper_workbook.pdf"
    with combined.open("wb") as stream:
        writer.write(stream)
    print(f"Combined: {len(writer.pages)} pages -> {combined}")


if __name__=="__main__":
    main()
