project_root <- if (file.exists(file.path("R", "utils.R"))) {
	"."
} else {
	file.path("..", "..")
}

source(file.path(project_root, "R", "utils.R"))
source(file.path(project_root, "R", "data-processing.R"))
source(file.path(project_root, "R", "grantham-distance.R"))
source(file.path(project_root, "R", "distance-calc.R"))
source(file.path(project_root, "R", "p-epitope-calculator.R"))
source(file.path(project_root, "R", "plots-and-tables.R"))
source(file.path(project_root, "R", "subtype-contrast.R"))
