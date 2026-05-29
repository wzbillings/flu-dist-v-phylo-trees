###
# Multiple sequence alignment helpers for H1 strains
# Zane Billings
###

align_h1_sequences <- function(clean_sequences, settings = make_analysis_settings()) {
	align_subtype_sequences(clean_sequences, subtype = "h1", settings = settings)
}

#### END OF FILE ####
