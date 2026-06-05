project_root <- if (file.exists(file.path("R", "utils.R"))) {
	"."
} else {
	file.path("..", "..")
}

source(file.path(project_root, "R", "utils.R"))
source(file.path(project_root, "R", "data-processing.R"))
source(file.path(project_root, "R", "alignment.R"))
source(file.path(project_root, "R", "provenance-audit.R"))
source(file.path(project_root, "R", "grantham-distance.R"))
source(file.path(project_root, "R", "distance-calc.R"))
source(file.path(project_root, "R", "cartography-diagnostics.R"))
source(file.path(project_root, "R", "panel-summary.R"))
source(file.path(project_root, "R", "p-epitope-calculator.R"))
source(file.path(project_root, "R", "tree-building.R"))
source(file.path(project_root, "R", "plots-and-tables.R"))
source(file.path(project_root, "R", "subtype-contrast.R"))
source(file.path(project_root, "R", "influence-analysis.R"))
source(file.path(project_root, "R", "sequence-sensitivity-audit.R"))
source(file.path(project_root, "R", "alignment-sensitivity.R"))
source(file.path(project_root, "R", "secondary-sequence-distance.R"))
