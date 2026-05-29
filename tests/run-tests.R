#!/usr/bin/env Rscript

find_script_path <- function() {
	args <- commandArgs(FALSE)
	file_arg <- args[grepl("^--file=", args)]
	if (length(file_arg) > 0) {
		return(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
	}
	if (file.exists(file.path("tests", "run-tests.R"))) {
		return(normalizePath(file.path("tests", "run-tests.R"), winslash = "/", mustWork = TRUE))
	}
	this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
	if (!is.na(this_file)) {
		return(this_file)
	}
	stop("Could not determine the location of tests/run-tests.R.", call. = FALSE)
}

write_lines <- function(lines, path) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	writeLines(lines, path, useBytes = TRUE)
}

format_duration <- function(start, end) {
	sprintf("%.2f seconds", as.numeric(difftime(end, start, units = "secs")))
}

as_markdown_code_block <- function(lines) {
	if (length(lines) == 0) {
		return(c("```text", "<no console output>", "```"))
	}
	c("```text", lines, "```")
}

run_logged_tests <- function() {
	if (!requireNamespace("testthat", quietly = TRUE)) {
		stop(
			"The testthat package is required to run tests. Restore or install it before running this script.",
			call. = FALSE
		)
	}
	library(testthat)

	script_path <- find_script_path()
	project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
	old_wd <- getwd()
	on.exit(setwd(old_wd), add = TRUE)
	setwd(project_root)

	run_started <- Sys.time()
	log_stamp <- format(run_started, "%Y%m%d-%H%M%S")
	log_path <- file.path("tests", "logs", paste0("test-run-", log_stamp, ".md"))
	test_files <- list.files(
		file.path("tests", "testthat"),
		pattern = "^test-.*[.]R$",
		full.names = TRUE
	)

	results <- list()
	active_file <- NA_character_

	test_that <- function(desc, code) {
		block_started <- Sys.time()
		status <- "passed"
		condition_message <- NA_character_
		console <- character()

		evaluated <- tryCatch(
			withCallingHandlers(
				capture.output(force(code), type = "output"),
				warning = function(w) {
					console <<- c(console, paste0("Warning: ", conditionMessage(w)))
					invokeRestart("muffleWarning")
				},
				message = function(m) {
					console <<- c(console, paste0("Message: ", conditionMessage(m)))
					invokeRestart("muffleMessage")
				}
			),
			skip = function(e) {
				status <<- "skipped"
				condition_message <<- conditionMessage(e)
				character()
			},
			error = function(e) {
				status <<- "failed"
				condition_message <<- conditionMessage(e)
				character()
			}
		)

		console <- c(console, evaluated)
		block_ended <- Sys.time()
		results[[length(results) + 1L]] <<- list(
			file = active_file,
			test = desc,
			status = status,
			message = condition_message,
			started = block_started,
			ended = block_ended,
			duration = format_duration(block_started, block_ended),
			console = console
		)

		invisible(NULL)
	}

	test_env <- new.env(parent = .GlobalEnv)
	test_env$test_that <- test_that

	helper_files <- list.files(
		file.path("tests", "testthat"),
		pattern = "^helper-.*[.]R$",
		full.names = TRUE
	)
	for (helper in helper_files) {
		source(helper, local = test_env)
	}

	for (test_file in test_files) {
		active_file <- normalizePath(test_file, winslash = "/", mustWork = TRUE)
		source(test_file, local = test_env)
	}

	run_ended <- Sys.time()
	statuses <- vapply(results, `[[`, character(1), "status")
	n_failed <- sum(statuses == "failed")
	n_skipped <- sum(statuses == "skipped")
	n_passed <- sum(statuses == "passed")

	lines <- c(
		"# Test Run Log",
		"",
		paste0("- Run started: ", format(run_started, "%Y-%m-%d %H:%M:%S %Z")),
		paste0("- Run ended: ", format(run_ended, "%Y-%m-%d %H:%M:%S %Z")),
		paste0("- Duration: ", format_duration(run_started, run_ended)),
		paste0("- Project root: `", project_root, "`"),
		paste0("- Result: ", if (n_failed == 0) "PASS" else "FAIL"),
		paste0("- Passed: ", n_passed),
		paste0("- Failed: ", n_failed),
		paste0("- Skipped: ", n_skipped),
		"",
		"## Test Results",
		"",
		"| Status | File | Test | Duration | Message |",
		"| --- | --- | --- | --- | --- |"
	)

	for (result in results) {
		file_label <- gsub("\\\\", "/", result$file)
		file_label <- sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/"), "", file_label)
		message <- if (is.na(result$message)) "" else gsub("[\r\n|]+", " ", result$message)
		lines <- c(
			lines,
			paste(
				result$status,
				paste0("`", file_label, "`"),
				gsub("[\r\n|]+", " ", result$test),
				result$duration,
				message,
				sep = " | "
			)
		)
	}

	lines <- c(lines, "")

	tests_with_console <- Filter(function(x) length(x$console) > 0, results)
	if (length(tests_with_console) > 0) {
		lines <- c(lines, "## Console Output", "")
		for (result in tests_with_console) {
			file_label <- gsub("\\\\", "/", result$file)
			file_label <- sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", project_root), "/"), "", file_label)
			lines <- c(
				lines,
				paste0("### ", result$status, ": ", result$test),
				"",
				paste0("File: `", file_label, "`"),
				"",
				as_markdown_code_block(result$console),
				""
			)
		}
	}

	write_lines(lines, log_path)
	message("Wrote test log: ", normalizePath(log_path, winslash = "/", mustWork = TRUE))

	if (n_failed > 0) {
		stop(n_failed, " test block(s) failed. See ", log_path, call. = FALSE)
	}

	invisible(normalizePath(log_path, winslash = "/", mustWork = TRUE))
}

tryCatch(
	run_logged_tests(),
	error = function(e) {
		message(conditionMessage(e))
		if (!interactive()) {
			quit(status = 1, save = "no")
		}
		stop(e)
	}
)
