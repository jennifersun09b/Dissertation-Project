# CVD lifestyle risk tool - Streamlit front end for the dissertation's
# predictive (Stream 1) and causal (Stream 2) analyses.
#
# Three stages, shown one at a time: overview -> questionnaire -> results.
#
# Copy is written for patients: each section leads with one short line, and the
# methodological detail (cohort, C-statistic, null-effect caveats, references)
# lives in "Learn more" expanders so it stays available without being in the way.
# Visual tokens and HTML components are in theme.py.
#
# Run with:  streamlit run app.py
# Model files are supplied by the analyst (see README.md); the app renders a
# clear notice for any file that is missing or unreadable.

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import model_io
import questions as q
import recommendations as rec
import theme

# Cohort context for the risk readout (Stream 1 primary run: def_CVD_AFTER,
# N=458,840, 27,454 events). Update if the supplied model used another cohort.
COHORT_EVENT_RATE = 0.0598
MODEL_AUC_NOTE = "internal validation C-statistic 0.718 (95% CI 0.712-0.724)"

st.set_page_config(
    page_title="CVD Lifestyle Risk Explorer",
    page_icon="🫀",
    layout="centered",
)
theme.inject()


@st.cache_resource
def load_models():
    predictive, p_err = model_io.load_predictive()
    single, s_err = model_io.load_single_ate()
    combined, c_err = model_io.load_combined_ate()
    return {
        "predictive": predictive, "predictive_error": p_err,
        "single": single, "single_error": s_err,
        "combined": combined, "combined_error": c_err,
    }


def model_notice(error):
    st.warning(f"**Model not loaded.** {error}", icon="📦")


# ---------------------------------------------------------------------------
# Stage 1 - overview
# ---------------------------------------------------------------------------
def render_hero():
    theme.hero(
        kicker="UK Biobank · 500,000 adults",
        title="Your heart, your habits",
        subtitle="See your cardiovascular risk, and which lifestyle changes "
                 "the evidence says would actually move it.",
        tags=("3 minutes", "9 questions", "Nothing is saved"),
    )


def render_intro():
    c1, c2 = st.columns(2, gap="medium")
    with c1:
        theme.card("🔮 Your risk", "Your chance of developing cardiovascular "
                                  "disease over the coming years.")
    with c2:
        theme.card("⚖️ What helps", "Which single change is estimated to lower "
                                   "that risk — and by how much.")

    st.markdown("")
    if st.button("Start  →", type="primary"):
        st.session_state.stage = "questions"
        st.rerun()

    with st.expander("How this works"):
        st.markdown(
            "**Two analyses, two different questions.**\n\n"
            "- **Prediction** — a statistical model estimates your probability "
            "of developing CVD from nine characteristics: age, sex, BMI, and "
            "six lifestyle factors. It tells you *how likely*, not *why*.\n"
            "- **Causation** — a causal machine learning analysis (causal "
            "forests) estimates the *effect* of changing each lifestyle factor, "
            "based on what actually happened to study participants who changed "
            "their habits. It tells you *what is worth changing*.\n\n"
            "Both are built on roughly half a million UK Biobank participants "
            f"followed for up to ~10 years. Predictive model: {MODEL_AUC_NOTE}. "
            f"Cohort event rate: {COHORT_EVENT_RATE * 100:.2f}%."
        )
    with st.expander("Limits and disclaimer"):
        st.markdown(rec.DISCLAIMER)


# ---------------------------------------------------------------------------
# Stage 2 - questionnaire
# ---------------------------------------------------------------------------
def render_questionnaire():
    theme.section("About you", "Six lifestyle areas, plus your basic details.")

    # Deliberately not an st.form: widgets inside a form do not rerun the script
    # when they change, so a follow-up question cannot appear in response to the
    # answer above it. The bordered container reproduces the form's card look.
    with st.container(border=True):
        theme.form_section("Basic details", first=True)
        c1, c2 = st.columns(2)
        age = c1.number_input("Age (years)", 18, 100, 55)
        sex_label = c2.radio("Sex", ["Female", "Male"], horizontal=True)
        h1, h2 = st.columns(2)
        height = h1.number_input("Height (cm)", 120.0, 220.0, 170.0, step=0.5)
        weight = h2.number_input("Weight (kg)", 35.0, 200.0, 72.0, step=0.5)

        theme.form_section("😴 Sleep")
        sleep_hours = st.number_input("Hours of sleep in 24 hours", 2, 16, 7)
        chronotype = st.radio(
            "Are you a morning person?",
            ["Yes, definitely or more so than evening", "No, more of an evening person"],
        )
        insomnia = st.radio(
            "Trouble falling asleep, or waking in the night?",
            ["Never or rarely", "Sometimes", "Usually"],
        )
        snore = st.radio("Do you snore?", ["No", "Yes"], horizontal=True)
        doze = st.radio("Doze off unintentionally during the day?",
                        ["Never or rarely", "Sometimes or often"])

        theme.form_section("🚬 Smoking")
        smoking = st.radio(
            "Do you smoke tobacco?",
            ["I have never smoked", "I used to smoke, but stopped", "I currently smoke"],
        )

        theme.form_section("🍷 Alcohol")
        alc_status = st.radio(
            "Do you drink alcohol?",
            ["I have never drunk alcohol", "I used to drink, but stopped", "I currently drink"],
        )
        # Frequency is only asked of current drinkers, as in UK Biobank (field
        # 1558 follows on from field 20117). score_alcohol ignores it for the
        # other two answers, so nothing is lost by not asking.
        alc_freq_label = None
        if alc_status == "I currently drink":
            alc_freq_label = st.selectbox(
                "How often do you drink?",
                list(q.ALCOHOL_FREQ_OPTIONS.values()),
                index=2,
            )

        theme.form_section("🥗 Diet")
        fruit_veg = st.number_input(
            "Fruit and vegetables per day (servings)", 0.0, 20.0, 3.0, step=0.5,
            help="1 serving ≈ 1 piece of fruit, or 3 heaped tablespoons of vegetables.")
        fish = st.number_input("Fish per week (servings)", 0.0, 21.0, 1.0, step=0.5)
        proc_meat = st.radio(
            "Processed meat", ["Twice a week or less", "More than twice a week"],
            help="Bacon, ham, sausages.")
        red_meat = st.radio(
            "Red meat", ["Five times a week or less", "More than five times a week"],
            help="Beef, lamb, pork.")

        theme.form_section("🏃 Physical activity", hint="In a typical week.")
        pa1, pa2 = st.columns(2)
        walk_days = pa1.number_input("Days walking ≥10 minutes", 0, 7, 3)
        walk_mins = pa2.number_input("Minutes walking, on those days", 0, 600, 30)
        mod_days = pa1.number_input(
            "Days of moderate activity", 0, 7, 2,
            help="For example cycling, or heavy housework.")
        mod_mins = pa2.number_input("Minutes of moderate activity", 0, 600, 30)
        vig_days = pa1.number_input(
            "Days of vigorous activity", 0, 7, 0,
            help="For example running, or fast sports.")
        vig_mins = pa2.number_input("Minutes of vigorous activity", 0, 600, 0)

        theme.form_section("🧠 Mental wellbeing", hint="Over the last two weeks.")
        phq1 = st.select_slider("Feeling down, depressed, or hopeless", q.FREQUENCY_SCALE)
        phq2 = st.select_slider("Little interest or pleasure in doing things", q.FREQUENCY_SCALE)
        gad1 = st.select_slider("Feeling nervous, anxious, or on edge", q.FREQUENCY_SCALE)
        gad2 = st.select_slider("Unable to stop or control worrying", q.FREQUENCY_SCALE)
        lonely = st.radio("Do you often feel lonely?", ["No", "Yes"], horizontal=True)

        submitted = st.button("See my results  →", type="primary")

    if submitted:
        bmi = weight / (height / 100) ** 2
        freq_key = (
            {v: k for k, v in q.ALCOHOL_FREQ_OPTIONS.items()}[alc_freq_label]
            if alc_freq_label is not None else None
        )
        raw = {
            "sleep": q.score_sleep(
                chronotype_morning=chronotype.startswith("Yes"),
                hours=sleep_hours,
                insomnia_rare=insomnia == "Never or rarely",
                snores=snore == "Yes",
                dozes=doze != "Never or rarely",
            ),
            "smoking": q.score_smoking(
                {"I have never smoked": "never",
                 "I used to smoke, but stopped": "previous",
                 "I currently smoke": "current"}[smoking]),
            "alcohol": q.score_alcohol(
                {"I have never drunk alcohol": "never",
                 "I used to drink, but stopped": "previous",
                 "I currently drink": "current"}[alc_status],
                freq_key),
            "diet": q.score_diet(
                fruit_veg_servings=fruit_veg,
                fish_weekly=fish,
                processed_meat_low=proc_meat.startswith("Twice"),
                red_meat_low=red_meat.startswith("Five"),
            ),
            "pa": q.score_physical(walk_days, walk_mins, mod_days, mod_mins,
                                   vig_days, vig_mins),
            "mental": q.score_mental(
                q.FREQUENCY_SCALE.index(phq1), q.FREQUENCY_SCALE.index(phq2),
                q.FREQUENCY_SCALE.index(gad1), q.FREQUENCY_SCALE.index(gad2),
                lonely == "Yes"),
        }
        st.session_state.answers = {
            "age": age, "sex": 1 if sex_label == "Male" else 0, "bmi": bmi,
            "raw": raw, "causal": q.raw_to_causal_levels(raw),
        }
        st.session_state.stage = "results"
        st.rerun()


# ---------------------------------------------------------------------------
# Stage 3 - results
# ---------------------------------------------------------------------------
def forest_figure(rows, x_title):
    """Horizontal point + 95% CI chart for effect estimates (in pp).

    `rows` are (short_label, detail, est, low, high, significant): the short
    label sits on the axis to keep the plot area wide, and the detail appears
    on hover.

    Polarity is encoded with a diverging pair - blue below zero (risk down),
    red above (risk up) - and gray for estimates whose interval crosses zero.
    Marker shape is a second, non-colour channel for that same distinction:
    filled = statistically clear, hollow = no clear effect.
    """
    fig = go.Figure()
    fig.add_vline(x=0, line_color=theme.GRAY, line_width=1)
    for label, detail, est, low, high, significant in rows:
        color = (theme.RED if est > 0 else theme.BLUE) if significant else theme.GRAY
        fig.add_trace(go.Scatter(
            x=[est], y=[label], mode="markers",
            marker=dict(
                size=11, color=color,
                symbol="circle" if significant else "circle-open",
                # filled marks get a surface ring so overlaps stay readable
                line=dict(width=2 if not significant else 1.5,
                          color=color if not significant else "#ffffff"),
            ),
            error_x=dict(type="data", array=[high - est],
                         arrayminus=[est - low], color=color,
                         thickness=2 if significant else 1.5),
            hovertemplate=(f"<b>{label}</b><br>{detail}<br>"
                           f"%{{x:+.2f}} pp "
                           f"(95% CI {low:+.2f} to {high:+.2f})<extra></extra>"),
            showlegend=False,
        ))
    fig.update_layout(
        height=max(200, 52 * len(rows) + 86),
        margin=dict(l=10, r=18, t=10, b=40),
        font=dict(family="Inter, -apple-system, sans-serif",
                  color=theme.INK_SOFT, size=12),
        xaxis=dict(title=x_title, gridcolor=theme.GRID, zeroline=False),
        yaxis=dict(autorange="reversed",
                   tickfont=dict(size=13, color=theme.INK)),
        plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
        hoverlabel=dict(font_family="Inter, -apple-system, sans-serif"),
    )
    return fig


def render_risk(models, answers):
    theme.section("Your risk", number="1")
    if models["predictive"] is None:
        model_notice(models["predictive_error"])
        return
    row = pd.DataFrame([q.raw_scores_to_model_row(
        answers["age"], answers["sex"], answers["bmi"], answers["raw"])])
    try:
        risk = float(models["predictive"].predict_proba(row)[0, 1])
    except Exception as exc:
        st.error(f"The predictive model failed on your inputs: {exc}")
        return

    ratio = risk / COHORT_EVENT_RATE
    higher = risk > COHORT_EVENT_RATE
    c1, c2 = st.columns([3, 2], gap="medium")
    with c1:
        theme.stat(
            "Your estimated risk", f"{risk * 100:.1f}%",
            sub="Chance of developing cardiovascular disease over about 10 years.",
            chip=f"{ratio:.1f}× the study average",
            chip_tone="up" if higher else "down",
        )
    with c2:
        theme.stat("Study average", f"{COHORT_EVENT_RATE * 100:.1f}%",
                   sub="UK Biobank participants.", plain=True)

    with st.expander("How to read this number"):
        st.markdown(
            f"This is an **associational** estimate of absolute risk over the "
            f"study's follow-up (up to ~10 years) — not a diagnosis, and not a "
            f"statement about what caused your risk. The model's "
            f"{MODEL_AUC_NOTE}, meaning it ranks higher- and lower-risk people "
            f"correctly about 72% of the time. Individual estimates carry more "
            f"uncertainty than that figure suggests."
        )


def render_single_causal(models, answers):
    theme.section("What could help",
                  "Estimated effect of making <b>one</b> change. "
                  "Left of the line means lower risk.", number="2")
    if models["single"] is None:
        model_notice(models["single_error"])
        return None

    table = models["single"]
    rows, cards = [], []
    for domain, level in answers["causal"].items():
        move = rec.improvement_transition(domain, level)
        if move is None:
            cards.append((domain, level, None, None))
            continue
        match = table[(table["domain"] == domain)
                      & (table["baseline"] == move[0])
                      & (table["destination"] == move[1])]
        if match.empty:
            cards.append((domain, level, move, None))
            continue
        r = match.iloc[0]
        desc = rec.TRANSITION_DESCRIPTIONS.get(
            (domain, *move), f"moving from level {move[0]} to {move[1]}")
        rows.append((rec.domain_title(domain), desc.capitalize(),
                     r["estimate"] * 100, r["low"] * 100,
                     r["high"] * 100, bool(r["sig"])))
        cards.append((domain, level, move, r))

    if rows:
        st.plotly_chart(
            forest_figure(rows, "Change in CVD risk (pp)"),
            use_container_width=True)
        theme.forest_legend()
        st.caption("Hover a row for the exact change and confidence interval.")
        with st.expander("What the chart shows"):
            st.markdown(
                "Each row is one lifestyle change, compared with people like "
                "you who kept that factor unchanged. Estimates are absolute "
                "risk differences in **percentage points (pp)**: −1 pp means "
                "one fewer case of CVD per 100 people. The horizontal bar is "
                "the 95% confidence interval — the range the true effect "
                "plausibly falls in. When that bar crosses zero, the study "
                "could not distinguish the effect from no effect at all."
            )
    else:
        st.success(
            "You are already at the healthiest measured level in every "
            "lifestyle area — there is no single change left to model.",
            icon="✅",
        )
    return cards


def _arm_label(raw):
    """Short axis label for a combined-analysis arm.

    The stored arm_label uses raw domain keys ("pa improves only"), which read
    badly on an axis. Map them to display titles and compress to "Diet only" /
    "Both improve"; the full wording stays in the hover detail.
    """
    text = str(raw)
    doms = [t for t in text.replace("+", " ").split() if t in q.DOMAIN_TITLES]
    titles = [rec.domain_title(d) for d in doms]
    if len(titles) >= 2:
        return "Both improve"
    if titles:
        return f"{titles[0]} only"
    return text[:1].upper() + text[1:]


def render_combined_causal(models, answers):
    theme.section("Two changes together",
                  "Does improving two factors at once beat improving one?",
                  number="3")
    if models["combined"] is None:
        model_notice(models["combined_error"])
        return

    improvable = [d for d, lvl in answers["causal"].items()
                  if rec.improvement_transition(d, lvl) is not None]
    table = models["combined"]
    scenario = table[table["model_type"] == "improvement"]
    pairs = [
        p for p in scenario["pair"].unique()
        if all(part in improvable for part in str(p).split("_"))
    ]
    if not pairs:
        st.info("Fewer than two of your lifestyle areas have room to improve, "
                "so no pair applies to you.")
        return

    nice = {p: " + ".join(rec.domain_title(part)
                          for part in str(p).split("_")) for p in pairs}
    choice = st.selectbox("Pick a pair", pairs, format_func=lambda p: nice[p])
    subset = scenario[scenario["pair"] == choice]
    rows = []
    for _, r in subset.iterrows():
        significant = (r["low"] > 0) or (r["high"] < 0)
        arm = str(r["arm_label"])
        rows.append((_arm_label(arm), f"{arm}, vs keeping both unchanged",
                     r["estimate"] * 100, r["low"] * 100,
                     r["high"] * 100, significant))
    st.plotly_chart(
        forest_figure(rows, "Change in CVD risk (pp)"),
        use_container_width=True)
    theme.forest_legend()
    with st.expander("What the chart shows"):
        st.markdown(
            "Each arm is compared with people who kept **both** factors "
            "unchanged: one factor improving, the other improving, or both at "
            "once. If 'both improve' is larger than the two single changes "
            "added together, the changes reinforce each other. In this study "
            "most combined effects were statistically indistinguishable "
            "from zero."
        )


def render_recommendations(cards):
    theme.section("Your plan", number="4")
    # Actionable domains first: what to change leads, what to maintain follows.
    ordered = sorted(cards or [], key=lambda c: c[2] is None)
    for domain, level, move, ate_row in ordered:
        title = rec.domain_title(domain)
        now = rec.level_label(domain, level)
        if move is None:
            theme.advice_card(title, now, "✅ Keep this up.",
                              rec.MAINTAIN_ADVICE[domain], done=True)
            continue
        desc = rec.TRANSITION_DESCRIPTIONS.get(
            (domain, *move), f"moving to level {move[1]}")
        if ate_row is not None:
            est, lo, hi = (ate_row["estimate"] * 100, ate_row["low"] * 100,
                           ate_row["high"] * 100)
            verdict = ("statistically clear" if bool(ate_row["sig"])
                       else "no clear effect in this study")
            detail = (f"Estimated effect {est:+.2f} pp "
                      f"(95% CI {lo:+.2f} to {hi:+.2f}) — {verdict}")
        else:
            detail = "No causal estimate available for this change"
        theme.advice_card(title, now, f"🎯 Try {desc}.",
                          rec.DOMAIN_ADVICE[domain], detail=detail)

    with st.expander("Reading these numbers honestly"):
        st.markdown(
            "In this study most single lifestyle changes showed small effects "
            "with confidence intervals crossing zero over the follow-up "
            "window, and a few apparent effects ran opposite to expectation — "
            "likely reverse causation, since people often change habits "
            "*because* their health changed. The clearest signal was harm from "
            "**taking up smoking**.\n\n"
            "Null short-term effects do not mean lifestyle is irrelevant: "
            "healthier baseline levels were consistently associated with lower "
            "risk, and the advice above reflects the wider evidence base and "
            "published guidelines rather than this study alone."
        )
    with st.expander("References and disclaimer"):
        st.markdown("\n".join(f"{i}. {r}"
                              for i, r in enumerate(rec.REFERENCES, 1)))
        st.markdown(f"**Disclaimer.** {rec.DISCLAIMER}")


def render_results(models):
    answers = st.session_state.answers

    with st.expander("Your profile as the models see it"):
        st.markdown("\n".join(
            f"- **{rec.domain_title(d)}**: {rec.level_label(d, lvl)}"
            for d, lvl in answers["causal"].items()))
        st.markdown(f"- **BMI**: {answers['bmi']:.1f} kg/m²")

    render_risk(models, answers)
    st.divider()
    cards = render_single_causal(models, answers)
    st.divider()
    render_combined_causal(models, answers)
    st.divider()
    render_recommendations(cards)

    st.markdown("")
    if st.button("←  Change my answers"):
        st.session_state.stage = "questions"
        st.rerun()
    theme.footnote(
        "Research demonstration, not a medical device. It does not provide "
        "medical advice, diagnosis, or treatment. Please discuss any health "
        "concerns with a qualified clinician."
    )


# ---------------------------------------------------------------------------
def main():
    if "stage" not in st.session_state:
        st.session_state.stage = "intro"
    models = load_models()

    missing = [k for k in ("predictive", "single", "combined") if models[k] is None]
    if missing:
        st.sidebar.warning(
            "Waiting for model files:\n" + "\n".join(f"- `{m}`" for m in (
                "predictive_model.joblib" if "predictive" in missing else None,
                "single_variable_ate.joblib" if "single" in missing else None,
                "combined_variable_ate.joblib" if "combined" in missing else None,
            ) if m),
            icon="📦",
        )
        if st.sidebar.button("Reload model files"):
            load_models.clear()
            st.rerun()

    render_hero()
    stage = st.session_state.stage
    if stage == "intro":
        render_intro()
    elif stage == "questions":
        render_questionnaire()
    elif stage == "results" and "answers" in st.session_state:
        render_results(models)


main()
