# Questionnaire definitions and scoring for the CVD risk tool.
#
# Every question here reproduces the UK Biobank touchscreen items that the
# dissertation's data_merge_new.ipynb used to build the nine modelling
# variables, so the answers map onto exactly the scales the models were
# trained on. Question wording is deliberately kept in this one module so it
# can be revised without touching the app logic.
#
# Two codings coexist (see cohort/build_table1_causal_cohorts.py):
#   raw scores    -> fed to the predictive model
#                    (sleep/smoking/physical 0-2, alcohol 0-4, diet 0-3, mental 0-3)
#   causal levels -> all six domains collapsed to 0/1/2, used to look up
#                    lifestyle-transition effects in the causal tables

# Raw score -> causal 0/1/2 level (identity for sleep/smoking/physical)
CAUSAL_RECODES = {
    "alcohol": {0: 0, 1: 0, 2: 1, 3: 2, 4: 2},
    "diet": {0: 0, 1: 0, 2: 1, 3: 2},
    "mental": {0: 0, 1: 1, 2: 2, 3: 2},
    "sleep": {0: 0, 1: 1, 2: 2},
    "smoking": {0: 0, 1: 1, 2: 2},
    "pa": {0: 0, 1: 1, 2: 2},
}

# Direction convention per domain: is a HIGHER causal level healthier?
HIGHER_IS_HEALTHIER = {
    "sleep": True, "diet": True, "pa": True,
    "mental": False, "alcohol": False, "smoking": False,
}

DOMAIN_TITLES = {
    "sleep": "Sleep", "diet": "Diet", "pa": "Physical activity",
    "mental": "Mental wellbeing", "alcohol": "Alcohol", "smoking": "Smoking",
}

# Human-readable label for each causal level, worded from the participant's side
CAUSAL_LEVEL_LABELS = {
    "sleep": {0: "poor sleep pattern", 1: "intermediate sleep pattern", 2: "healthy sleep pattern"},
    "diet": {0: "few healthy diet targets met", 1: "some healthy diet targets met", 2: "most healthy diet targets met"},
    "pa": {0: "low activity", 1: "moderate activity", 2: "high activity"},
    "mental": {0: "no current symptoms", 1: "one symptom area", 2: "two or more symptom areas"},
    "alcohol": {0: "non-drinker / former drinker", 1: "occasional drinker", 2: "regular drinker"},
    "smoking": {0: "never smoked", 1: "previous smoker", 2: "current smoker"},
}


def score_sleep(chronotype_morning, hours, insomnia_rare, snores, dozes):
    """UKB fields 1180/1160/1200/1210/1220 -> healthy-sleep category 0/1/2.

    Five components score 1 point each when healthy; the 0-5 score maps to
    category 2 when >= 4, category 0 when <= 2, else 1.
    """
    points = 0
    points += 1 if chronotype_morning else 0
    points += 1 if 7 <= hours <= 8 else 0
    points += 1 if insomnia_rare else 0
    points += 1 if not snores else 0
    points += 1 if not dozes else 0
    if points >= 4:
        return 2
    if points <= 2:
        return 0
    return 1


def score_smoking(status):
    """UKB field 20116: 'never' -> 0, 'previous' -> 1, 'current' -> 2."""
    return {"never": 0, "previous": 1, "current": 2}[status]


def score_alcohol(status, frequency=None):
    """UKB fields 20117 (status) + 1558 (intake frequency) -> raw score 0-4.

    `frequency` is only consulted for current drinkers, so the caller may pass
    None when the follow-up question was not asked.
    """
    if status == "never":
        return 0
    if status == "previous":
        return 1
    return {  # current drinkers, by frequency
        "special_occasions": 2, "one_to_three_monthly": 2,
        "once_or_twice_weekly": 3, "three_or_four_weekly": 3,
        "daily": 4,
    }[frequency]


def score_diet(fruit_veg_servings, fish_weekly, processed_meat_low, red_meat_low):
    """AHA food-group targets -> raw diet score 0-3 (1 point per target met).

    Targets: fruit+vegetables >= 4.5 servings/day; fish >= 2 servings/week;
    processed meat <= 2/week AND red meat <= 5/week.
    """
    score = 0
    score += 1 if fruit_veg_servings >= 4.5 else 0
    score += 1 if fish_weekly >= 2 else 0
    score += 1 if processed_meat_low and red_meat_low else 0
    return score


def score_physical(walk_days, walk_mins, mod_days, mod_mins, vig_days, vig_mins):
    """IPAQ short form (UKB fields 864/874, 884/894, 904/914) -> category 0/1/2."""
    met = (3.3 * walk_days * walk_mins
           + 4.0 * mod_days * mod_mins
           + 8.0 * vig_days * vig_mins)
    total_days = walk_days + mod_days + vig_days
    if (vig_days >= 3 and met >= 1500) or (met >= 3000 and total_days >= 7):
        return 2
    if ((vig_days >= 3 and vig_mins >= 20)
            or (mod_days >= 5 and mod_mins >= 30)
            or (walk_days >= 5 and walk_mins >= 30)
            or (total_days >= 5 and met >= 600)):
        return 1
    return 0


def score_mental(phq_low_mood, phq_disinterest, gad_nervous, gad_worry, lonely):
    """PHQ-2-style + GAD-2-style items (0-3 each) + loneliness -> raw score 0-3."""
    score = 0
    score += 1 if (phq_low_mood + phq_disinterest) >= 3 else 0
    score += 1 if (gad_nervous + gad_worry) >= 3 else 0
    score += 1 if lonely else 0
    return score


FREQUENCY_SCALE = [  # PHQ/GAD response scale, index = points
    "Not at all",
    "Several days",
    "More than half the days",
    "Nearly every day",
]

ALCOHOL_FREQ_OPTIONS = {
    "daily": "Daily or almost daily",
    "three_or_four_weekly": "Three or four times a week",
    "once_or_twice_weekly": "Once or twice a week",
    "one_to_three_monthly": "One to three times a month",
    "special_occasions": "Only on special occasions",
}


def raw_scores_to_model_row(age, sex, bmi, raw):
    """Assemble the 9-column DataFrame row the predictive pipeline expects."""
    return {
        "age": float(age),
        "bmi": float(bmi),
        "alcohol_score_0.0": float(raw["alcohol"]),
        "diet_score_0.0": float(raw["diet"]),
        "mental_score_0.0": float(raw["mental"]),
        "sleep_category_0.0": float(raw["sleep"]),
        "smoking_category_0.0": float(raw["smoking"]),
        "physical_category_0.0": float(raw["pa"]),
        "sex": float(sex),
    }


def raw_to_causal_levels(raw):
    """Collapse raw scores to the 0/1/2 levels used by the causal analyses."""
    return {domain: CAUSAL_RECODES[domain][score] for domain, score in raw.items()}
