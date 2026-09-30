#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Export the dashboard's pipeline flow chart as a poster figure (SVG + PDF + 300-dpi PNG).

The chart is built in docs/js/about.js from the release's own numbers, so the poster and the website show
the same diagram. This serves docs/ on a loopback port, lets the About tab render in the light theme, takes
the <svg id="pipeline"> as-is (its <style> falls back to the light palette when the page's theme tokens are
absent), and prints it to vector PDF and PNG with headless Chromium.

Needs a Python with Playwright (the system python3 here; not the atac_hdac env):
  python3 phase2/figures/export_pipeline.py [--out poster/figures]
"""
from __future__ import annotations

import argparse
import asyncio
import functools
import http.server
import os
import threading

from playwright.async_api import async_playwright

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.abspath(os.path.join(HERE, "..", ".."))
DOCS = os.environ.get("SECACTS_DOCS", os.path.join(SECACTS, "docs"))   # a staged copy, e.g. the next release


def serve():
    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *a):
            pass
    handler = functools.partial(Quiet, directory=DOCS)
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)       # loopback, OS-assigned port
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv


async def main(out):
    os.makedirs(out, exist_ok=True)
    srv = serve()
    url = f"http://127.0.0.1:{srv.server_address[1]}/#about"
    async with async_playwright() as p:
        b = await p.chromium.launch()
        ctx = await b.new_context(viewport={"width": 1400, "height": 1000})
        await ctx.add_init_script("try{localStorage.setItem('secacts-theme','light')}catch(e){}")
        pg = await ctx.new_page()
        await pg.goto(url)
        await pg.wait_for_selector("#pipeline")
        svg = await pg.eval_on_selector("#pipeline", "e => e.outerHTML")
        vb = await pg.eval_on_selector("#pipeline", "e => ({width: e.viewBox.baseVal.width, height: e.viewBox.baseVal.height})")
        await ctx.close()
        # standalone: fixed size in user units (px), no page tokens -> the <style> fallbacks apply
        svg = svg.replace('width="100%"', f'width="{vb["width"]}" height="{vb["height"]}"', 1)
        open(os.path.join(out, "fig0_pipeline.svg"), "w").write('<?xml version="1.0" encoding="UTF-8"?>\n' + svg)
        page = await b.new_page(viewport={"width": int(vb["width"]), "height": int(vb["height"])},
                                device_scale_factor=300 / 96 * 2.4)   # >= 300 dpi up to ~50 cm wide
        await page.set_content(f"<html><body style='margin:0;background:#fff'>{svg}</body></html>")
        await page.locator("#pipeline").screenshot(path=os.path.join(out, "fig0_pipeline.png"))
        await page.pdf(path=os.path.join(out, "fig0_pipeline.pdf"), width=f"{vb['width']}px",
                       height=f"{vb['height'] + 2}px", print_background=True, page_ranges="1")
        await b.close()
    srv.shutdown()
    print(f"[fig] fig0_pipeline -> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(SECACTS, "poster", "figures"))
    asyncio.run(main(ap.parse_args().out))
