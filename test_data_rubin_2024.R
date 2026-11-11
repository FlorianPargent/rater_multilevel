# Checks for prep_data_rubin_2024.R
#
# Sources the prep script and verifies the raw OSF files, the reshaping, the
# handling of missing ratings, and the saved output. Section 6 documents why
# the trauma type of the index event is not included in the cleaned data.
# Every check stops with an error if it fails. Blocks marked "Informational"
# only print summaries.
#
# Objects used from the prep script: data_dir, rating_files, ratings_raw,
# rater_key, crita_full (with missing ratings), crita_long (final), out_path.

source("prep_data_rubin_2024.R")

library(irr)

# 1. Raw rating files --------------------------------------------------------

raw_files <- rating_files |>
  set_names(basename) |>
  map(\(f) as_tibble(readRDS(f)))

# Each file: 100 participants, unique ids, one version label matching the
# file name.
iwalk(raw_files, \(d, f) {
  stopifnot(
    nrow(d) == 100,
    anyDuplicated(d$id) == 0,
    n_distinct(d$LEC) == 1,
    tolower(unique(d$LEC)) == str_match(f, "_(lec[13]p?)\\.rds$")[, 2]
  )
})

# Grad and phd files of the same version contain the same participants, so the
# design is fully crossed within each version.
for (v in c("lec1", "lec1p", "lec3", "lec3p")) {
  stopifnot(setequal(
    raw_files[[paste0("ratings_grad_", v, ".rds")]]$id,
    raw_files[[paste0("ratings_phd_", v, ".rds")]]$id
  ))
}

# Raw rating values are only 0, 1 or NA, before any type conversion. This also
# catches strings like "NA" or "" in the character column, which would have
# been converted to NA silently. Exactly 2 values are missing in the raw data,
# so the conversion in the prep script introduced no additional NAs.
raw_values <- raw_files |>
  map(\(d) unlist(map(select(d, starts_with("CRITA")), as.character))) |>
  unlist(use.names = FALSE)

stopifnot(
  all(raw_values %in% c("0", "1") | is.na(raw_values)),
  sum(is.na(raw_values)) == 2
)

# 2. Reshaped data with missing ratings (crita_full) -----------------------

# Version label in the files agrees with the version from the file name.
stopifnot(all(tolower(ratings_raw$LEC) == ratings_raw$version))

# Fully crossed: 400 subjects x 6 raters, 100 subjects per version x rater,
# each subject x rater combination exactly once.
stopifnot(
  nrow(crita_full) == 2400,
  n_distinct(crita_full$subject) == 400,
  n_distinct(crita_full$rater) == 6,
  all(count(crita_full, subject)$n == 6),
  all(count(crita_full, rater)$n == 400),
  all(count(crita_full, version, rater)$n == 100),
  nrow(distinct(crita_full, subject, rater)) == 2400
)

# Each version maps to exactly one combination of design factors.
design <- crita_full |>
  distinct(version, n_index, part2) |>
  mutate(across(everything(), as.character)) |>
  arrange(version)

stopifnot(
  identical(design$version, c("lec1", "lec1p", "lec3", "lec3p")),
  identical(design$n_index, c("single", "single", "up_to_3", "up_to_3")),
  identical(design$part2, c("standard", "extended", "standard", "extended"))
)

# Rater labels: three per group, mapped to the expected names; rater_group is
# constant within rater.
raters <- mutate(rater_key, across(everything(), as.character))

stopifnot(
  identical(raters$rater, c("G1", "G2", "G3", "P1", "P2", "P3")),
  identical(raters$rater_name, c("alban", "emma", "mikael", "angkaw", "powers", "smits")),
  nrow(distinct(crita_full, rater, rater_group)) == 6
)

stopifnot(all(crita_full$crit_a %in% c(0L, 1L) | is.na(crita_full$crit_a)))

# 3. Missing ratings -------------------------------------------------------

# Exactly two missing ratings, both from P3 (smits) in lec3.
missing_ratings <- filter(crita_full, is.na(crit_a))

stopifnot(
  nrow(missing_ratings) == 2,
  all(missing_ratings$rater == "P3"),
  raters$rater_name[raters$rater == "P3"] == "smits",
  all(missing_ratings$version == "lec3"),
  setequal(as.character(missing_ratings$subject), c("lec3_31494", "lec3_32098"))
)

# Informational: ratings of the other five raters for these two subjects
# (at the time of writing all were 1, i.e. no sign of an ambiguous case).
crita_full |>
  filter(subject %in% missing_ratings$subject) |>
  select(subject, rater, crit_a) |>
  pivot_wider(names_from = rater, values_from = crit_a) |>
  print()

# 4. Reproduce Fleiss' kappas reported in the paper --------------------------

# irr 0.85 needs a base data.frame or matrix: on a tibble, ratings[, i] stays a
# tibble and as.factor() fails ("cannot xtfrm data frames").
# kappam.fleiss() drops subjects with any missing rating, as in the original
# analysis, so crita_full is used here (398 subjects overall).
fleiss_kappa <- function(data) {
  data |>
    select(subject, rater, crit_a) |>
    pivot_wider(names_from = rater, values_from = crit_a) |>
    select(-subject) |>
    as.data.frame() |>
    kappam.fleiss()
}

kappa_check <- tibble(
  comparison = c(
    "all", "single index", "up to 3 indices", "standard part 2",
    "extended part 2", "graduate students", "psychologists"
  ),
  data = list(
    crita_full,
    filter(crita_full, n_index == "single"),
    filter(crita_full, n_index == "up_to_3"),
    filter(crita_full, part2 == "standard"),
    filter(crita_full, part2 == "extended"),
    filter(crita_full, rater_group == "grad_student"),
    filter(crita_full, rater_group == "psychologist")
  ),
  reported = c(0.428, 0.456, 0.397, 0.415, 0.442, 0.524, 0.380)
) |>
  mutate(
    fit = map(data, fleiss_kappa),
    subjects = map_int(fit, "subjects"),
    kappa = map_dbl(fit, "value"),
    matches = abs(round(kappa, 3) - reported) < 1e-8
  ) |>
  select(comparison, subjects, kappa, reported, matches)

print(kappa_check)

stopifnot(
  kappa_check$subjects[kappa_check$comparison == "all"] == 398,
  all(kappa_check$matches)
)

# 5. Final data (crita_long) -----------------------------------------------

stopifnot(
  nrow(crita_long) == 2398,
  n_distinct(crita_long$subject) == 400,
  all(crita_long$crit_a %in% c(0L, 1L)),
  identical(
    names(crita_long),
    c("subject", "id", "version", "n_index", "part2", "rater",
      "rater_group", "crit_a")
  ),
  # no missing values in any column
  !anyNA(crita_long),
  # crita_long equals crita_full without the two missing ratings
  identical(crita_long, filter(crita_full, !is.na(crit_a)))
)

# Saved csv has the expected header and row count, and reading it back with
# read_crita() (which recreates the factors) reproduces the object in memory,
# including column types and factor levels.
csv_lines <- readLines(out_path)
stopifnot(
  csv_lines[1] == paste0('"', paste(names(crita_long), collapse = '","'), '"'),
  length(csv_lines) == nrow(crita_long) + 1
)

crita_csv <- read_crita(out_path)

stopifnot(
  identical(crita_csv, crita_long),
  identical(levels(crita_csv$subject), levels(crita_long$subject)),
  identical(levels(crita_csv$rater), levels(crita_long$rater)),
  identical(map_chr(crita_csv, \(x) class(x)[1]), map_chr(crita_long, \(x) class(x)[1]))
)

# Informational: proportion of "Criterion A met" per rater (rater severity).
crita_long |>
  summarise(n = n(), prop_met = mean(crit_a), .by = c(rater_group, rater)) |>
  print()

# 6. Documentation: why trauma type is excluded ----------------------------

# lec_trauma_type_agreement.rds holds the LEC-5 category (LEC_code) of the
# index trauma for lec1 and lec1p (200 rows), together with the six ratings
# and the authors' agreement flags. It has no participant id. The authors'
# RR1 script built it by joining grad and phd ratings and then attaching the
# trauma codes with bind_cols(), i.e. by row position. The checks below show:
#   (a) the file's rows (ratings) align with the grad rating files by position,
#   (b) for lec1, trauma code counts and the agreement flag are consistent with
#       the paper and the ratings,
#   (c) for lec1p, neither is consistent.
# (a) does not verify that the codes belong to the right participants: a
# misordered code vector would leave the ratings intact. (b) is necessary but
# not sufficient, since codes shuffled among participants keep the counts.
# Because the link between codes and participants cannot be verified, trauma
# type is not used. If the lec1p checks below ever fail, the OSF file has
# changed and this decision should be revisited.

trauma <- readRDS(file.path(data_dir, "lec_trauma_type_agreement.rds")) |>
  as_tibble() |>
  mutate(
    version = tolower(LEC),
    row_in_version = row_number(),
    .by = LEC
  )

grad_order <- c("lec1", "lec1p") |>
  map(\(v) {
    readRDS(file.path(data_dir, paste0("ratings_grad_", v, ".rds"))) |>
      as_tibble() |>
      transmute(
        version = v,
        row_in_version = row_number(),
        subject = paste(v, id, sep = "_")
      )
  }) |>
  list_rbind()

trauma <- left_join(trauma, grad_order, by = c("version", "row_in_version"))

# (a) Row alignment: all 6 ratings of all 200 rows match crita_full for the
# subject assigned by position.
trauma_alignment <- trauma |>
  select(subject, starts_with("CRITA")) |>
  pivot_longer(
    starts_with("CRITA"),
    names_to = "rater_name",
    names_prefix = "CRITA_RATING_0_1_",
    values_to = "crit_a_trauma"
  ) |>
  left_join(select(rater_key, rater_name, rater), by = "rater_name") |>
  mutate(rater = factor(rater, levels = levels(crita_full$rater))) |>
  left_join(select(crita_full, subject, rater, crit_a), by = c("subject", "rater"))

stopifnot(
  nrow(trauma_alignment) == 1200,
  !anyNA(trauma_alignment$subject),
  all(trauma_alignment$crit_a == trauma_alignment$crit_a_trauma)
)

# Within-group agreement flags (all 3 grad / all 3 phd raters agree) are
# consistent with the group sums in both versions.
stopifnot(
  all(trauma$agreement_y_n_grad == as.numeric(trauma$sum_grad_agreement %in% c(0, 3))),
  all(trauma$agreement_y_n_phd == as.numeric(trauma$sum_phd_agreement %in% c(0, 3)))
)

# (b) + (c) Overall flag (all 6 raters agree) vs. the ratings: consistent for
# lec1, contradicted for 14 lec1p participants (non-unanimous but flagged as
# agreeing). Whether this is a row-order mismatch or a different flag
# definition cannot be determined. The authors' trauma-type brms models used
# agreement_y_n as outcome.
unanimous <- crita_full |>
  filter(version %in% c("lec1", "lec1p")) |>
  summarise(agree = as.numeric(n_distinct(crit_a) == 1), .by = subject) |>
  mutate(subject = as.character(subject))

agreement_check <- trauma |>
  select(subject, version, agreement_y_n) |>
  left_join(unanimous, by = "subject")

print(count(agreement_check, version, agree, agreement_y_n))

stopifnot(
  with(filter(agreement_check, version == "lec1"), all(agree == agreement_y_n)),
  with(filter(agreement_check, version == "lec1p"), sum(agree != agreement_y_n) == 14)
)

# (b) lec1 trauma code counts match Table 3 (Sample 1) and the RR1 script.
reported_lec1 <- tibble(
  LEC_code = c(1, 3, 4, 6, 7, 8, 10, 12, 14, 15, 17),
  n_reported = c(5, 13, 1, 9, 4, 21, 1, 18, 8, 4, 16)
)

lec1_counts <- trauma |>
  filter(version == "lec1") |>
  count(LEC_code, name = "n_file") |>
  full_join(reported_lec1, by = "LEC_code")

stopifnot(all(lec1_counts$n_file == lec1_counts$n_reported))

# (c) lec1p trauma code counts do NOT match the RR1 script (7 of 12
# categories differ). Table 3, Sample 2 agrees with the script except for
# code 15, where the paper's 15 would make the column sum to 105; the
# script's 10 is used here.
reported_lec1p <- tibble(
  LEC_code = c(1, 3, 4, 6, 7, 8, 11, 12, 13, 14, 15, 17),
  n_reported = c(9, 15, 1, 2, 2, 14, 1, 17, 1, 11, 10, 17)
)

lec1p_counts <- trauma |>
  filter(version == "lec1p") |>
  count(LEC_code, name = "n_file") |>
  full_join(reported_lec1p, by = "LEC_code") |>
  mutate(across(c(n_file, n_reported), \(x) coalesce(as.numeric(x), 0))) |>
  arrange(LEC_code)

print(lec1p_counts)

stopifnot(sum(lec1p_counts$n_file != lec1p_counts$n_reported) == 7)

message("All checks passed.")
