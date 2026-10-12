#!/usr/bin/env python3
"""Figs. S5 and S6 as the manuscript prints them: panel A (central) above panel B (least-stable bound).

Stacks the outputs of make_Figure_S5_predominance_300K.py and make_Figure_S6_predominance_600K.py
and puts each panel letter in a white strip above its panel, so the letter never touches the plot.
Writes artifacts/figures/Figure_S5.png and Figure_S6.png.
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import matplotlib

FIG = Path(__file__).resolve().parents[1] / "artifacts" / "figures"
FONT = ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts/ttf/DejaVuSans-Bold.ttf"), 72)
STRIP = 96

for fig in ("S5", "S6"):
    panels = [Image.open(FIG / f"Figure_{fig}{p}.png").convert("RGB") for p in "AB"]
    width = max(p.width for p in panels)
    out = Image.new("RGB", (width, sum(p.height + STRIP for p in panels)), "white")
    draw, y = ImageDraw.Draw(out), 0
    for letter, panel in zip("AB", panels):
        draw.text((18, y + 10), letter, font=FONT, fill="black")
        out.paste(panel, (0, y + STRIP))
        y += STRIP + panel.height
    out.save(FIG / f"Figure_{fig}.png", dpi=(190, 190))
    print("wrote", FIG / f"Figure_{fig}.png")
