# Recommendation content and transition logic for the results page.

from questions import CAUSAL_LEVEL_LABELS, DOMAIN_TITLES, HIGHER_IS_HEALTHIER


def improvement_transition(domain, level):
    """The (baseline, destination) improvement move available from a level.

    Returns None when the user is already at the healthiest level, or when the
    only move is not causally interpretable (smoking previous -> never is a
    reporting artefact, not a behaviour change).
    """
    if HIGHER_IS_HEALTHIER[domain]:
        return (level, level + 1) if level < 2 else None
    if domain == "smoking":
        return (2, 1) if level == 2 else None  # cessation is the only real move
    return (level, level - 1) if level > 0 else None


# What the improvement move means in plain words, keyed by (domain, baseline, destination)
TRANSITION_DESCRIPTIONS = {
    ("sleep", 0, 1): "moving from a poor to an intermediate sleep pattern",
    ("sleep", 1, 2): "moving from an intermediate to a healthy sleep pattern",
    ("diet", 0, 1): "meeting more of the healthy diet targets",
    ("diet", 1, 2): "meeting most of the healthy diet targets",
    ("pa", 0, 1): "moving from low to moderate physical activity",
    ("pa", 1, 2): "moving from moderate to high physical activity",
    ("mental", 1, 0): "resolving current mental health symptoms",
    ("mental", 2, 1): "reducing mental health symptom burden",
    ("alcohol", 1, 0): "stopping regular alcohol use",
    ("alcohol", 2, 1): "cutting down from regular to occasional drinking",
    ("smoking", 2, 1): "quitting smoking",
}

# Guideline-grounded advice shown per domain. Kept separate from the causal
# estimates on purpose: the dissertation found most transition effects null or
# paradoxical, so day-to-day advice leans on established guidelines while the
# causal panel reports the study's own estimates honestly.
DOMAIN_ADVICE = {
    "sleep": (
        "Aim for 7-8 hours per night with a regular schedule. Address frequent "
        "insomnia, loud snoring, or daytime sleepiness with your GP, as these "
        "can indicate treatable sleep disorders."
    ),
    "diet": (
        "Work toward at least 4.5 servings of fruit and vegetables a day, two "
        "servings of fish a week, and limited processed and red meat, in line "
        "with American Heart Association targets."
    ),
    "pa": (
        "Guidelines recommend at least 150 minutes of moderate or 75 minutes "
        "of vigorous activity per week. Regular brisk walking counts, and any "
        "increase from a low baseline helps."
    ),
    "mental": (
        "Persistent low mood, anxiety, or loneliness are worth discussing with "
        "your GP. Mental wellbeing is linked with cardiovascular health, and "
        "effective support is available."
    ),
    "alcohol": (
        "UK guidance advises no more than 14 units per week, spread over "
        "several days. If you drink regularly, cutting frequency is a "
        "reasonable first step."
    ),
    "smoking": (
        "Quitting smoking is the single most protective lifestyle change for "
        "cardiovascular health. NHS stop-smoking services roughly triple the "
        "chance of quitting successfully."
    ),
}

MAINTAIN_ADVICE = {
    "sleep": "Your sleep pattern is already in the healthy range. Keep your schedule consistent.",
    "diet": "You already meet most healthy diet targets. Keep it up.",
    "pa": "You are already in the high activity group. Keep it up.",
    "mental": "You report no current symptom areas. Stay connected and seek help early if that changes.",
    "alcohol": "You are in the lowest alcohol category. Maintaining this is protective.",
    "smoking": (
        "You do not currently smoke. Staying smoke-free matters: in this "
        "study, taking up smoking showed the single largest harmful effect "
        "of any lifestyle change examined."
    ),
}

REFERENCES = [
    "UK Biobank cohort: Sudlow C, et al. UK Biobank: an open access resource "
    "for identifying the causes of a wide range of complex diseases of middle "
    "and old age. PLoS Med. 2015;12(3):e1001779.",
    "Causal forests: Wager S, Athey S. Estimation and inference of "
    "heterogeneous treatment effects using random forests. JASA. "
    "2018;113(523):1228-1242.",
    "Generalized random forests: Athey S, Tibshirani J, Wager S. Generalized "
    "random forests. Ann Statist. 2019;47(2):1148-1178.",
    "Prediction model reporting: Collins GS, et al. TRIPOD+AI statement. BMJ. "
    "2024;385:e078378.",
    "Diet targets: Lloyd-Jones DM, et al. Defining and setting national goals "
    "for cardiovascular health promotion (Life's Simple 7). Circulation. "
    "2010;121(4):586-613.",
    "Physical activity: WHO guidelines on physical activity and sedentary "
    "behaviour. Geneva: World Health Organization; 2020.",
    "Alcohol: UK Chief Medical Officers' low risk drinking guidelines. "
    "Department of Health; 2016.",
    "Smoking cessation: NICE guideline NG209, Tobacco: preventing uptake, "
    "promoting quitting and treating dependence; 2021.",
]

DISCLAIMER = (
    "This tool is a research demonstration built on a dissertation analysis "
    "of UK Biobank data. It is not a medical device and does not provide "
    "medical advice, diagnosis, or treatment. Estimates are for a middle-aged "
    "UK cohort and may not transfer to other populations. Please discuss any "
    "health concerns with a qualified clinician."
)


def domain_title(domain):
    return DOMAIN_TITLES[domain]


def level_label(domain, level):
    return CAUSAL_LEVEL_LABELS[domain][level]
