###
# Multiple sequence alignment helpers for H3 strains
# Zane Billings
###

align_h3_sequences <- function(clean_sequences, settings = make_analysis_settings()) {
	align_subtype_sequences(clean_sequences, subtype = "h3", settings = settings)
}

#### END OF FILE ####
