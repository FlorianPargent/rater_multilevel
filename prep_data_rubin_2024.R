# Prepare Criterion A rating data from Rubin et al. (2024)
#
# Source:
#   Rubin, M., Alban Foulser, A., Siegel, E., Smits, J. A. J., Powers, M. B.,
#   Angkaw, A., & Telch, M. J. (2024). Low reliability when determining
#   Criterion A for posttraumatic stress disorder from self-report descriptions
#   of traumatic events: The need for transparent methods. Psychological Trauma:
#   Theory, Research, Practice, and Policy, 16(3), 435-442.
#   https://doi.org/10.1037/tra0001477
# Data: https://osf.io/dbxpw/ (downloaded to dbxpw-osfstorage-archive/)
#
# Design: 400 participants (100 per LEC-5 version) x 6 raters (3 clinical
# psychology graduate students, 3 licensed psychologists), fully crossed. Each
# rater judged whether a participant's self-reported trauma description met
# DSM-5 PTSD Criterion A. Coding 1 = met, 0 = not met is inferred from the
# column names (CRITA_RATING_0_1_*); it is not documented in the OSF materials.
#
# Sample numbering in the paper is inconsistent (Methods text vs. Table 3), so
# design factors are derived from the file names, consistent with Table 3 and
# the authors' scripts:
#   lec1  = single index trauma, standard part 2
#   lec1p = single index trauma, extended part 2 (3 extra follow-up questions)
#   lec3  = up to 3 index traumas, standard part 2
#   lec3p = up to 3 index traumas, extended part 2
#
# Output: data/rubin_2024_crita.csv (object `crita_long`), one row per
# subject x rater. The csv does not store factor information; read_crita()
# below recreates the factor columns (subject, version, n_index, part2, rater,
# rater_group) with the same levels as in crita_long.
#   subject      participant key: version + original id (ids may repeat
#                across versions; the authors joined on id and LEC)
#   id           original participant id
#   version      lec1, lec1p, lec3, lec3p
#   n_index      number of index traumas participants could describe
#   part2        LEC-5 part 2: standard or extended
#   rater        anonymised rater: G1-G3 graduate students, P1-P3 psychologists
#   rater_group  grad_student or psychologist
#   crit_a       1 = Criterion A met, 0 = not met
#
# Not included: trauma type of the index event. The OSF file
# lec_trauma_type_agreement.rds contains LEC-5 trauma categories for lec1 and
# lec1p only (none exist for lec3/lec3p, where raters did not indicate which of
# up to 3 described events they judged). The file has no participant id, so
# codes can only be matched by row position, and their link to participants
# cannot be verified: the authors attached the codes to their ratings with
# bind_cols(), so a misordered code vector would leave the ratings in the file
# intact. For lec1p there is direct evidence of a problem: 7 of 12 category
# counts differ from the paper's Table 3 and the authors' RR1 script, and the
# file's "all 6 raters agree" flag (agreement_y_n) contradicts the ratings for
# 14 of 100 participants. lec1 passes both checks, but matching counts cannot
# rule out codes being shuffled among participants. The evidence is reproduced
# in section 6 of test_data_rubin_2024.R.
#
# Checks: test_data_rubin_2024.R

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

data_dir <- "dbxpw-osfstorage-archive/Data"
out_path <- "data/rubin_2024_crita.csv"

# Read the saved csv and restore the column types of crita_long. id stays
# character (keeps any leading zeros); factor levels are fixed explicitly.
read_crita <- function(path) {
  read.csv(
    path,
    colClasses = c(
      subject = "character", id = "character", version = "character",
      n_index = "character", part2 = "character", rater = "character",
      rater_group = "character", crit_a = "integer"
    )
  ) |>
    as_tibble() |>
    mutate(
      subject = factor(subject),
      version = factor(version, levels = c("lec1", "lec1p", "lec3", "lec3p")),
      n_index = factor(n_index, levels = c("single", "up_to_3")),
      part2 = factor(part2, levels = c("standard", "extended")),
      rater = factor(rater, levels = c(paste0("G", 1:3), paste0("P", 1:3))),
      rater_group = factor(rater_group, levels = c("grad_student", "psychologist"))
    )
}

# 1. Read rating files -----------------------------------------------------

# One file per rater group x version, each in wide format: id, three rating
# columns CRITA_RATING_0_1_<rater name>, and the version label LEC.
rating_files <- list.files(
  data_dir,
  pattern = "^ratings_(grad|phd)_lec[13]p?\\.rds$",
  full.names = TRUE
)
stopifnot(length(rating_files) == 8)

# Rating column types differ across files (integer, numeric, and character for
# smits in ratings_phd_lec3p.rds). as.character() first so that factors would
# not be converted to their level codes.
to_num <- function(x) as.numeric(as.character(x))

# Pivot each file separately because grad and phd files have different rater
# columns.
read_ratings <- function(path) {
  readRDS(path) |>
    as_tibble() |>
    mutate(
      id = as.character(id),
      LEC = as.character(LEC),
      across(starts_with("CRITA"), to_num)
    ) |>
    pivot_longer(
      starts_with("CRITA"),
      names_to = "rater_name",
      names_prefix = "CRITA_RATING_0_1_",
      values_to = "crit_a"
    ) |>
    mutate(file = basename(path))
}

# 2. Design variables from file names --------------------------------------

ratings_raw <- rating_files |>
  map(read_ratings) |>
  list_rbind() |>
  mutate(
    version = str_match(file, "_(lec[13]p?)\\.rds$")[, 2],
    rater_group = factor(
      str_match(file, "ratings_(grad|phd)_")[, 2],
      levels = c("grad", "phd"),
      labels = c("grad_student", "psychologist")
    ),
    n_index = factor(
      if_else(str_starts(version, "lec3"), "up_to_3", "single"),
      levels = c("single", "up_to_3")
    ),
    part2 = factor(
      if_else(str_ends(version, "p"), "extended", "standard"),
      levels = c("standard", "extended")
    ),
    subject = paste(version, id, sep = "_")
  )

# 3. Anonymised rater labels -----------------------------------------------

# Alphabetical within group: G1 alban, G2 emma, G3 mikael,
# P1 angkaw, P2 powers, P3 smits.
rater_key <- ratings_raw |>
  distinct(rater_group, rater_name) |>
  arrange(rater_group, rater_name) |>
  mutate(
    rater = paste0(if_else(rater_group == "grad_student", "G", "P"), row_number()),
    .by = rater_group
  )

# 4. Assemble long data (missing ratings still included) -------------------

crita_full <- ratings_raw |>
  left_join(rater_key, by = c("rater_group", "rater_name")) |>
  transmute(
    subject = factor(subject),
    id,
    version = factor(version, levels = c("lec1", "lec1p", "lec3", "lec3p")),
    n_index,
    part2,
    rater = factor(rater, levels = rater_key$rater),
    rater_group,
    crit_a = as.integer(crit_a)
  ) |>
  arrange(version, subject, rater)

# 5. Remove missing ratings (last processing step) -------------------------

# Two ratings are missing in the original data: ratings_phd_lec3.rds, rater
# smits (P3), subjects lec3_31494 and lec3_32098. Neither the paper nor the
# OSF scripts mention them; irr::kappam.fleiss() in the original analysis
# silently dropped both participants entirely (398 subjects). Here only the two
# missing rating rows are removed, so both participants keep their other five
# ratings. Result: 2398 rows, 400 subjects.
stopifnot(sum(is.na(crita_full$crit_a)) == 2)

crita_long <- crita_full |>
  filter(!is.na(crit_a))

# 6. Save --------------------------------------------------------------------

# Never overwrite: if the file exists it must be identical to the new result
# (after restoring column types).
if (!dir.exists(dirname(out_path))) {
  dir.create(dirname(out_path))
}

if (file.exists(out_path)) {
  if (!identical(read_crita(out_path), crita_long)) {
    stop(
      out_path, " exists and differs from the current result. ",
      "Rename or remove it manually before saving a new version."
    )
  }
  message(out_path, " already exists and is identical; not overwritten.")
} else {
  write.csv(crita_long, out_path, row.names = FALSE)
}