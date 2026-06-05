#!/usr/bin/env python
"""Render the canonical Markdown docs in this folder to PDF artifacts.

Simple reportlab (Platypus) Markdown renderer — handles the subset used in our
docs: ATX headings (#/##/###), paragraphs, fenced code blocks (```), pipe tables,
bullet/numbered lists, blockquotes, and inline **bold** / `code` / [links](url).
The `.md` files remain canonical; the PDFs are generated artifacts.

Output: scripts/greigite/artifacts/<name>.pdf  (one per *.md here).
Run:    uv run python scripts/greigite/md_to_pdf.py
        (or: pip install reportlab; python scripts/greigite/md_to_pdf.py)

DejaVu TTFs are registered for Greek/subscript/arrow glyphs (Δ, σ, Θ, ₃, →, …);
falls back to Helvetica/Courier if they are not installed.
"""

from __future__ import annotations

import re
from pathlib import Path

from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    ListFlowable,
    ListItem,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "artifacts"

# --- fonts ----------------------------------------------------------------- #
_DEJAVU = Path("/usr/share/fonts/truetype/dejavu")
FONT, FONT_B, FONT_M = "Helvetica", "Helvetica-Bold", "Courier"
try:
    pdfmetrics.registerFont(TTFont("DejaVu", str(_DEJAVU / "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVu-Bold", str(_DEJAVU / "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVuMono", str(_DEJAVU / "DejaVuSansMono.ttf")))
    FONT, FONT_B, FONT_M = "DejaVu", "DejaVu-Bold", "DejaVuMono"
except Exception:  # noqa: BLE001 - fall back to built-in fonts
    pass

_ss = getSampleStyleSheet()
BODY = ParagraphStyle(
    "body",
    parent=_ss["BodyText"],
    fontName=FONT,
    fontSize=9,
    leading=12.5,
    alignment=TA_LEFT,
    spaceAfter=4,
)
H1 = ParagraphStyle(
    "h1",
    parent=BODY,
    fontName=FONT_B,
    fontSize=17,
    leading=21,
    spaceBefore=10,
    spaceAfter=8,
)
H2 = ParagraphStyle(
    "h2",
    parent=BODY,
    fontName=FONT_B,
    fontSize=13,
    leading=16,
    spaceBefore=10,
    spaceAfter=5,
)
H3 = ParagraphStyle(
    "h3",
    parent=BODY,
    fontName=FONT_B,
    fontSize=11,
    leading=14,
    spaceBefore=8,
    spaceAfter=4,
)
CODE = ParagraphStyle(
    "code",
    parent=BODY,
    fontName=FONT_M,
    fontSize=7.6,
    leading=9.4,
    backColor=colors.HexColor("#f4f4f4"),
    borderPadding=4,
    textColor=colors.HexColor("#202020"),
)
CELL = ParagraphStyle("cell", parent=BODY, fontSize=7.4, leading=9)
CELL_H = ParagraphStyle("cellh", parent=CELL, fontName=FONT_B)


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def inline(s: str) -> str:
    """Markdown inline -> reportlab mini-markup (escape first, then re-add tags)."""
    s = esc(s)
    s = re.sub(r"`([^`]+)`", rf'<font face="{FONT_M}" size="8">\1</font>', s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<link href="\2" color="blue">\1</link>', s)
    return s


def split_row(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def make_table(rows: list[str]):
    header = split_row(rows[0])
    body = [split_row(r) for r in rows[2:]]  # rows[1] is the |---| separator
    data = [[Paragraph(inline(c), CELL_H) for c in header]]
    data += [[Paragraph(inline(c), CELL) for c in r] for r in body]
    t = Table(data, repeatRows=1, hAlign="LEFT")
    t.setStyle(
        TableStyle(
            [
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#999999")),
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8e8e8")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ]
        )
    )
    return t


def md_to_flowables(text: str) -> list:
    lines = text.splitlines()
    flow: list = []
    i, n = 0, len(lines)
    bullets: list = []

    def flush_bullets():
        nonlocal bullets
        if bullets:
            flow.append(
                ListFlowable(
                    [
                        ListItem(Paragraph(inline(b), BODY), leftIndent=10)
                        for b in bullets
                    ],
                    bulletType="bullet",
                    start="•",
                    leftIndent=12,
                )
            )
            flow.append(Spacer(1, 3))
            bullets = []

    while i < n:
        ln = lines[i]
        # fenced code
        if ln.lstrip().startswith("```"):
            flush_bullets()
            i += 1
            buf = []
            while i < n and not lines[i].lstrip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1  # closing fence
            flow.append(Preformatted("\n".join(buf) or " ", CODE))
            flow.append(Spacer(1, 4))
            continue
        # table (line is a pipe row and next line is a |---| separator)
        if (
            ln.strip().startswith("|")
            and i + 1 < n
            and re.match(r"^\s*\|?[ :|-]+\|", lines[i + 1])
        ):
            flush_bullets()
            tbl = [ln]
            i += 1
            while i < n and lines[i].strip().startswith("|"):
                tbl.append(lines[i])
                i += 1
            flow.append(make_table(tbl))
            flow.append(Spacer(1, 5))
            continue
        # headings
        m = re.match(r"^(#{1,3})\s+(.*)$", ln)
        if m:
            flush_bullets()
            lvl = len(m.group(1))
            flow.append(Paragraph(inline(m.group(2)), {1: H1, 2: H2, 3: H3}[lvl]))
            i += 1
            continue
        # bullets / numbered
        mb = re.match(r"^\s*[-*]\s+(.*)$", ln)
        mn = re.match(r"^\s*\d+\.\s+(.*)$", ln)
        if mb or mn:
            bullets.append((mb or mn).group(1))
            i += 1
            continue
        # blockquote
        if ln.strip().startswith(">"):
            flush_bullets()
            q = ln.strip()[1:].strip()
            flow.append(
                Paragraph(
                    inline(q),
                    ParagraphStyle(
                        "q",
                        parent=BODY,
                        leftIndent=10,
                        textColor=colors.HexColor("#555555"),
                        borderColor=colors.HexColor("#cccccc"),
                    ),
                )
            )
            i += 1
            continue
        # blank
        if not ln.strip():
            flush_bullets()
            flow.append(Spacer(1, 4))
            i += 1
            continue
        # paragraph
        flush_bullets()
        flow.append(Paragraph(inline(ln), BODY))
        i += 1

    flush_bullets()
    return flow


def convert(md: Path) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf = OUT_DIR / (md.stem + ".pdf")
    doc = SimpleDocTemplate(
        str(pdf),
        pagesize=letter,
        title=md.name,
        leftMargin=0.7 * inch,
        rightMargin=0.7 * inch,
        topMargin=0.7 * inch,
        bottomMargin=0.7 * inch,
    )
    doc.build(md_to_flowables(md.read_text()))
    return pdf


def main() -> None:
    for md in sorted(HERE.glob("*.md")):
        pdf = convert(md)
        print(f"{md.name}  ->  {pdf.relative_to(HERE)}  ({pdf.stat().st_size} B)")


if __name__ == "__main__":
    main()
