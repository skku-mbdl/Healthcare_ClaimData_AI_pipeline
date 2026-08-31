# Minimal run logger: every message goes to the console and to a
# run-specific log file under logs/, named with the same run_id used for
# the results/<run_id> output folder so a run's logs and outputs are easy
# to line up after the fact.
#
# Verbatim copy of CKM_Drug/{Japan,Korea}/scripts/utils/logging_utils.R
# (the two were already byte-identical) -- shared here so both countries'
# Stage 1 (drug screen) pipeline use exactly one copy.

new_logger <- function(run_id, logs_dir) {
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
  log_path <- file.path(logs_dir, paste0("pipeline_", run_id, ".log"))

  log_msg <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
    message(line)
    cat(line, "\n", file = log_path, append = TRUE)
  }

  list(log = log_msg, path = log_path)
}
