# Run folders (runs/<name>/) for viz/app.R: upload barcodes, start and stop
# run_pipeline.sh, read its status. The app container gets RUN_AS, DOCKER_GID and
# either DATA_CONTAINER (desktop app) or HOST_PROJECT (start.sh).

MAIN_PROJECT <- "(main project)"

STEPS <- c(
  "01" = "BLAST database", "02" = "Kraken2 database",
  "03" = "Read QC", "04" = "Filter and orient", "05" = "Kraken2",
  "06" = "Clusters and consensus (hours)", "07" = "ITSx", "08" = "OTUs",
  "09" = "BLAST and taxonomy", "10" = "Organisms", "11" = "Alignments", "12" = "Tree"
)

valid_run_name <- function(name) grepl("^[A-Za-z0-9_-]+$", name)

run_dir <- function(project, run) {
  if (identical(run, MAIN_PROJECT)) project else file.path(project, "runs", run)
}

list_runs <- function(project) {
  dirs <- list.dirs(file.path(project, "runs"), full.names = FALSE, recursive = FALSE)
  sort(dirs[valid_run_name(dirs)])
}

# Runs that have pictures to show; the main project first if it has results.
runs_with_results <- function(project) {
  has <- function(dir) file.exists(file.path(dir, "results/10_organisms/organisms.tsv"))
  runs <- Filter(function(run) has(run_dir(project, run)), list_runs(project))
  c(if (has(project)) MAIN_PROJECT, runs)
}

run_files <- function(project, run) {
  files <- list.files(file.path(run_dir(project, run), "raw_data"),
                      pattern = "\\.fastq\\.gz$", full.names = TRUE)
  if (!length(files)) return(data.frame(file = character(), MB = numeric()))
  data.frame(file = basename(files), MB = round(file.size(files) / 1e6, 1))
}

add_files <- function(project, run, uploads) {
  ok <- grepl("\\.fastq\\.gz$", uploads$name)
  target <- file.path(run_dir(project, run), "raw_data")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  file.copy(uploads$datapath[ok], file.path(target, uploads$name[ok]), overwrite = TRUE)
  list(added = uploads$name[ok], skipped = uploads$name[!ok])
}

docker <- function(...) {
  suppressWarnings(system2("docker", c(...), stdout = TRUE, stderr = TRUE))
}

runner_name <- function(run) paste0("fungiflow_runner_", run)

# Containers of this run: the runner and the step it is on.
run_containers <- function(run) {
  names <- docker("ps", "--format", "{{.Names}}")
  pattern <- paste0("^(", runner_name(run), "|fungiflow_", run, "_[0-9]{2})$")
  grep(pattern, names, value = TRUE)
}

is_running <- function(run) runner_name(run) %in% run_containers(run)

# The desktop app sets DATA_CONTAINER (a stopped container holding the /data and
# Docker socket mounts); start.sh sets HOST_PROJECT instead.
start_run <- function(project, run) {
  run_as <- Sys.getenv("RUN_AS")
  group <- Sys.getenv("DOCKER_GID")
  data_container <- Sys.getenv("DATA_CONTAINER")
  host <- Sys.getenv("HOST_PROJECT", project)
  if (!nzchar(run_as)) stop("RUN_AS is not set. Start the app with bash start.sh or the desktop app.")
  mount <- if (nzchar(data_container)) {
    c("--volumes-from", data_container, "-e", paste0("DATA_CONTAINER=", data_container))
  } else {
    c("-v", paste0(host, ":/data"), "-v", "/var/run/docker.sock:/var/run/docker.sock",
      "-e", paste0("HOST_PROJECT=", host))
  }
  log <- file.path("runs", run, "logs", "pipeline.log")
  dir.create(file.path(project, dirname(log)), recursive = TRUE, showWarnings = FALSE)
  docker("rm", "-f", runner_name(run))
  out <- docker(
    "run", "-d", "--rm", "--name", runner_name(run),
    "--user", run_as, if (nzchar(group)) c("--group-add", group),
    "-e", "HOME=/tmp", "-e", paste0("RUN_AS=", run_as), mount,
    "-w", "/data", "fungiflow-viz",
    "bash", "-c", shQuote(sprintf("bash run_pipeline.sh %s > %s 2>&1", run, log))
  )
  if (!is.null(attr(out, "status"))) stop(paste(out, collapse = "\n"))
  invisible(TRUE)
}

stop_run <- function(project, run) {
  containers <- run_containers(run)
  if (length(containers)) docker("rm", "-f", containers)
  status <- file.path(run_dir(project, run), "status.tsv")
  cat(sprintf("stopped\tstopped\t%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = status, append = TRUE)
}

# One row per step with its last state.
run_status <- function(project, run) {
  table <- data.frame(step = names(STEPS), what = unname(STEPS), state = "", time = "")
  path <- file.path(run_dir(project, run), "status.tsv")
  if (file.exists(path) && file.size(path) > 0) {
    lines <- read.delim(path, header = FALSE, colClasses = "character",
                        col.names = c("step", "state", "time"))
    for (i in seq_len(nrow(lines))) {
      row <- match(lines$step[i], table$step)
      if (!is.na(row)) table[row, c("state", "time")] <- lines[i, c("state", "time")]
    }
    if (tail(lines$state, 1) == "stopped") table$state[table$state == "running"] <- "stopped"
  }
  table[!(table$step %in% c("01", "02") & table$state == ""), ]
}

run_log <- function(project, run, lines = 40) {
  path <- file.path(run_dir(project, run), "logs", "pipeline.log")
  if (!file.exists(path)) return("No log yet.")
  text <- unlist(strsplit(readLines(path, warn = FALSE), "\r"))
  # Progress counters ("... 42%") would push everything else off the screen.
  progress <- grepl("[0-9]+%\\s*$", text)
  progress[length(progress)] <- FALSE
  paste(tail(text[!progress], lines), collapse = "\n")
}
