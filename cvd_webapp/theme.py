# Presentation layer for the CVD Lifestyle Risk Explorer.
#
# Holds the "Vivid Lime Tech" design tokens, the stylesheet, and small HTML
# component helpers. Static blocks (hero, stat cards, advice cards) are emitted
# as plain HTML rather than styled Streamlit containers, so their appearance
# does not depend on Streamlit's internal DOM. Tokens mirror
# .streamlit/config.toml - keep the two in sync.

import html as _html

import streamlit as st

# --- Design tokens ---------------------------------------------------------
INK = "#0B3B2E"        # headings, primary text
INK_SOFT = "#4A6B5D"   # secondary text
INK_MUTED = "#7C9488"  # captions, footnotes
PRIMARY = "#00A86B"    # actions
PRIMARY_DARK = "#00875A"
LIME = "#B8F04A"       # accent bars, hero kicker
CANVAS = "#F4FBF0"
SURFACE = "#FFFFFF"
BORDER = "#DCEFD2"
HERO_FROM = "#0B3B2E"
HERO_TO = "#14634A"

# Chart ink. The forest plot encodes POLARITY, so it uses a diverging pair
# (blue = risk down, red = risk up) with gray as the neutral "no clear effect"
# midpoint. Blue/red passes CVD separation on white (dE 21.6 protan) where a
# green/red pair would fail, and it stays distinct from the green UI chrome.
BLUE, RED, GRAY, GRID = "#2a78d6", "#e34948", "#898781", "#e1e0d9"

_CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

/* ---- page shell ---- */
[data-testid="stMainBlockContainer"] {{
  padding-top: 1.25rem;
  padding-bottom: 4rem;
  max-width: 900px;
}}
[data-testid="stHeader"] {{ background: transparent; }}
h1, h2, h3, h4, h5 {{ letter-spacing: -0.02em; color: {INK}; }}

/* ---- hero ---- */
.hero {{
  background: linear-gradient(135deg, {HERO_FROM} 0%, {HERO_TO} 100%);
  border-radius: 24px;
  padding: 2.6rem 2.2rem 2.4rem;
  margin: 0 0 1.6rem;
  position: relative;
  overflow: hidden;
}}
.hero::after {{  /* soft lime glow, purely decorative */
  content: ""; position: absolute; right: -70px; top: -70px;
  width: 240px; height: 240px; border-radius: 50%;
  background: radial-gradient(circle, rgba(184,240,74,.22), transparent 68%);
}}
.hero-kicker {{
  color: {LIME}; font-size: .74rem; font-weight: 700;
  letter-spacing: .13em; text-transform: uppercase; margin-bottom: .7rem;
}}
.hero h1 {{
  color: #fff !important; font-size: 2.6rem; font-weight: 800;
  line-height: 1.08; letter-spacing: -0.035em; margin: 0 0 .6rem;
}}
.hero p {{
  color: rgba(255,255,255,.80); font-size: 1.06rem;
  line-height: 1.5; margin: 0; max-width: 30rem;
}}
.hero-tags {{ display: flex; flex-wrap: wrap; gap: .5rem; margin-top: 1.4rem; }}
.hero-tag {{
  background: rgba(255,255,255,.11); border: 1px solid rgba(184,240,74,.30);
  color: #EAFBDD; border-radius: 999px;
  padding: .34rem .82rem; font-size: .82rem; font-weight: 500;
}}

/* ---- section headings ---- */
.sec {{ margin: 0 0 .9rem; }}
.sec-head {{ display: flex; align-items: center; gap: .6rem; }}
.sec-num {{
  background: {LIME}; color: {INK}; font-size: .8rem; font-weight: 800;
  width: 1.55rem; height: 1.55rem; border-radius: 8px;
  display: flex; align-items: center; justify-content: center; flex: none;
}}
.sec-bar {{ background: {LIME}; width: .3rem; height: 1.5rem; border-radius: 99px; flex: none; }}
.sec-title {{ font-size: 1.4rem; font-weight: 800; color: {INK}; letter-spacing: -0.025em; }}
.sec-lead {{
  color: {INK_SOFT}; font-size: .96rem; line-height: 1.5;
  margin: .5rem 0 0; max-width: 40rem;
}}

/* ---- questionnaire group headings (own markup, so the rule above the first
       group can be suppressed without depending on Streamlit's DOM) ---- */
.qsec {{
  font-size: 1.06rem; font-weight: 800; color: {INK};
  letter-spacing: -0.02em; margin: 1.7rem 0 .2rem;
  padding-top: 1.1rem; border-top: 1px solid #EEF6E8;
}}
.qsec.first {{ margin-top: 0; padding-top: 0; border-top: none; }}
.qsec-hint {{ font-size: .84rem; font-weight: 400; color: {INK_MUTED}; margin-top: .15rem; }}

/* ---- generic card ---- */
.card {{
  background: {SURFACE}; border-radius: 20px; padding: 1.15rem 1.3rem;
  border: 1px solid {BORDER}; border-left: 4px solid {LIME};
  box-shadow: 0 1px 3px rgba(11,59,46,.05); height: 100%;
}}
.card-title {{ font-size: 1rem; font-weight: 700; color: {INK}; margin-bottom: .3rem; }}
.card-body {{ font-size: .92rem; line-height: 1.5; color: {INK_SOFT}; }}

/* ---- stat card (hero number) ---- */
.stat {{
  background: {SURFACE}; border-radius: 20px; padding: 1.25rem 1.4rem;
  border: 1px solid {BORDER}; border-left: 4px solid {LIME};
  box-shadow: 0 1px 3px rgba(11,59,46,.05); height: 100%;
}}
.stat.plain {{ border-left-color: {BORDER}; background: #FBFEF8; }}
.stat-label {{
  font-size: .74rem; font-weight: 700; letter-spacing: .09em;
  text-transform: uppercase; color: {INK_MUTED}; margin-bottom: .45rem;
}}
.stat-value {{
  font-size: 2.9rem; font-weight: 800; line-height: 1;
  letter-spacing: -0.04em; color: {INK}; font-variant-numeric: tabular-nums;
}}
.stat.plain .stat-value {{ color: {INK_SOFT}; font-size: 2.3rem; }}
.stat-sub {{ font-size: .87rem; color: {INK_SOFT}; margin-top: .5rem; line-height: 1.45; }}
.stat-chip {{
  display: inline-block; margin-top: .6rem; border-radius: 999px;
  padding: .2rem .6rem; font-size: .78rem; font-weight: 600;
}}
.chip-up {{ background: #FDECEC; color: #A32B2B; }}
.chip-down {{ background: #E8F8DC; color: {PRIMARY_DARK}; }}

/* ---- advice cards ---- */
.rcard {{
  background: {SURFACE}; border-radius: 18px; padding: 1.05rem 1.25rem;
  border: 1px solid {BORDER}; border-left: 4px solid {LIME};
  box-shadow: 0 1px 3px rgba(11,59,46,.05); margin-bottom: .7rem;
}}
.rcard.done {{ border-left-color: {PRIMARY}; }}
.rcard-top {{
  display: flex; align-items: baseline; gap: .55rem;
  flex-wrap: wrap; margin-bottom: .45rem;
}}
.rcard-title {{ font-size: 1.06rem; font-weight: 800; color: {INK}; letter-spacing: -0.02em; }}
.rcard-now {{
  background: #F1F7EC; color: {INK_SOFT}; border-radius: 999px;
  padding: .16rem .58rem; font-size: .76rem; font-weight: 500;
}}
.rcard-do {{ font-size: .96rem; font-weight: 600; color: {INK}; line-height: 1.45; }}
.rcard-tip {{ font-size: .89rem; color: {INK_SOFT}; line-height: 1.5; margin-top: .4rem; }}
.rcard-num {{
  font-size: .82rem; color: {INK_MUTED}; margin-top: .35rem;
  font-variant-numeric: tabular-nums;
}}

/* ---- pill buttons ---- */
.stButton button, [data-testid="stFormSubmitButton"] button {{
  border-radius: 999px !important; font-weight: 700 !important;
  padding: .58rem 1.6rem !important; border: none !important;
  transition: transform .12s ease, box-shadow .12s ease;
}}
[data-testid="stBaseButton-primary"],
[data-testid="stBaseButton-primaryFormSubmit"] {{
  background: {PRIMARY} !important; color: #fff !important;
  box-shadow: 0 2px 10px rgba(0,168,107,.30) !important;
}}
[data-testid="stBaseButton-primary"]:hover,
[data-testid="stBaseButton-primaryFormSubmit"]:hover {{
  background: {PRIMARY_DARK} !important;
  transform: translateY(-1px);
  box-shadow: 0 5px 16px rgba(0,168,107,.36) !important;
}}
[data-testid="stBaseButton-secondary"] {{
  background: {SURFACE} !important; color: {PRIMARY_DARK} !important;
  box-shadow: inset 0 0 0 1.5px {BORDER} !important;
}}

/* ---- questionnaire card / inputs ----
   The questionnaire is a bordered st.container rather than an st.form (a form
   would stop follow-up questions appearing in response to the answer above).
   The child chain anchors the rule to that one container, so other vertical
   blocks - which share the same test id - keep their default styling. */
[data-testid="stVerticalBlockBorderWrapper"]:has(
  > div > [data-testid="stVerticalBlock"] > [data-testid="stElementContainer"] .qsec.first
) {{
  background: {SURFACE}; border: 1px solid {BORDER} !important;
  border-radius: 22px !important; padding: 1.5rem 1.6rem !important;
  box-shadow: 0 1px 3px rgba(11,59,46,.05);
}}
[data-testid="stWidgetLabel"] p {{ font-weight: 500 !important; color: {INK} !important; }}

/* ---- sidebar (model-file status; dev-facing) ---- */
[data-testid="stSidebarContent"] {{ padding-top: 2.5rem; }}
[data-testid="stSidebar"] code {{
  white-space: normal; overflow-wrap: anywhere; font-size: .74rem;
}}

/* ---- expanders ("Learn more" drawers) ---- */
[data-testid="stExpander"] details {{
  background: transparent; border: 1px solid {BORDER} !important;
  border-radius: 14px !important;
}}
[data-testid="stExpander"] summary {{
  font-size: .89rem !important; font-weight: 600 !important;
  color: {PRIMARY_DARK} !important;
}}
[data-testid="stExpander"] summary:hover {{ color: {INK} !important; }}
[data-testid="stExpander"] p, [data-testid="stExpander"] li {{
  font-size: .88rem; color: {INK_SOFT}; line-height: 1.55;
}}

/* ---- misc ---- */
[data-testid="stCaptionContainer"] p {{ color: {INK_MUTED} !important; font-size: .84rem; }}
hr, [data-testid="stDivider"] hr {{ border-color: #E3F1D9 !important; }}
.footnote {{ font-size: .78rem; color: {INK_MUTED}; line-height: 1.5; margin-top: 1.5rem; }}

/* ---- chart legend (colour + shape, so identity is never colour-alone) ---- */
.legend {{
  display: flex; flex-wrap: wrap; gap: 1.1rem;
  font-size: .82rem; color: {INK_SOFT}; margin: .1rem 0 .2rem .2rem;
}}
.legend span {{ display: inline-flex; align-items: center; gap: .38rem; }}
.legend i {{
  width: .62rem; height: .62rem; border-radius: 50%; flex: none;
  display: inline-block;
}}
.lg-down {{ background: {BLUE}; }}
.lg-up {{ background: {RED}; }}
.lg-null {{ background: transparent; border: 2px solid {GRAY}; }}
</style>
"""


def inject():
    """Install the stylesheet. Call once, right after st.set_page_config."""
    st.markdown(_CSS, unsafe_allow_html=True)


def _esc(text):
    return _html.escape(str(text))


# --- component helpers -----------------------------------------------------
# Each returns nothing and writes directly. HTML is emitted as one line: blank
# lines or 4-space indents inside st.markdown would be parsed as markdown.

def hero(kicker, title, subtitle, tags=()):
    tag_html = "".join(
        f'<span class="hero-tag">{_esc(t)}</span>' for t in tags)
    st.markdown(
        f'<div class="hero"><div class="hero-kicker">{_esc(kicker)}</div>'
        f'<h1>{_esc(title)}</h1><p>{_esc(subtitle)}</p>'
        f'<div class="hero-tags">{tag_html}</div></div>',
        unsafe_allow_html=True,
    )


def section(title, lead=None, number=None):
    """Results-page heading. Numbered badge when `number` is given, else a bar."""
    badge = (f'<span class="sec-num">{_esc(number)}</span>' if number
             else '<span class="sec-bar"></span>')
    lead_html = f'<p class="sec-lead">{lead}</p>' if lead else ""
    st.markdown(
        f'<div class="sec"><div class="sec-head">{badge}'
        f'<span class="sec-title">{_esc(title)}</span></div>'
        f'{lead_html}</div>',
        unsafe_allow_html=True,
    )


def form_section(label, hint=None, first=False):
    """Group heading inside the questionnaire, with a rule above (not the first)."""
    hint_html = f'<div class="qsec-hint">{_esc(hint)}</div>' if hint else ""
    st.markdown(
        f'<div class="qsec{" first" if first else ""}">{_esc(label)}'
        f'{hint_html}</div>',
        unsafe_allow_html=True,
    )


def card(title, body):
    """Small feature card. `body` may contain inline HTML."""
    st.markdown(
        f'<div class="card"><div class="card-title">{_esc(title)}</div>'
        f'<div class="card-body">{body}</div></div>',
        unsafe_allow_html=True,
    )


def stat(label, value, sub=None, chip=None, chip_tone="down", plain=False):
    """Hero number. `chip_tone` is "up" (worse) or "down" (better)."""
    parts = [
        f'<div class="stat{" plain" if plain else ""}">',
        f'<div class="stat-label">{_esc(label)}</div>',
        f'<div class="stat-value">{_esc(value)}</div>',
    ]
    if chip:
        parts.append(
            f'<span class="stat-chip chip-{chip_tone}">{_esc(chip)}</span>')
    if sub:
        parts.append(f'<div class="stat-sub">{sub}</div>')
    parts.append("</div>")
    st.markdown("".join(parts), unsafe_allow_html=True)


def advice_card(title, now_label, action_html, tip, detail=None, done=False):
    """One lifestyle domain: where you are, what to change, why."""
    parts = [
        f'<div class="rcard{" done" if done else ""}">',
        f'<div class="rcard-top"><span class="rcard-title">{_esc(title)}</span>'
        f'<span class="rcard-now">Now: {_esc(now_label)}</span></div>',
        f'<div class="rcard-do">{action_html}</div>',
    ]
    if detail:
        parts.append(f'<div class="rcard-num">{detail}</div>')
    parts.append(f'<div class="rcard-tip">{_esc(tip)}</div></div>')
    st.markdown("".join(parts), unsafe_allow_html=True)


def footnote(text):
    st.markdown(f'<div class="footnote">{_esc(text)}</div>',
                unsafe_allow_html=True)


def forest_legend():
    """Legend for the forest plot: colour AND marker shape."""
    st.markdown(
        '<div class="legend">'
        '<span><i class="lg-down"></i>Lowers risk</span>'
        '<span><i class="lg-up"></i>Raises risk</span>'
        '<span><i class="lg-null"></i>No clear effect</span>'
        '</div>',
        unsafe_allow_html=True,
    )
