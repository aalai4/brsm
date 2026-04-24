# Declare NSE symbols used in ggplot2 aesthetics.
utils::globalVariables(c(".data"))

.brsm_require_ggplot2 <- function() {
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		stop("package 'ggplot2' is required for plotting functions.")
	}
}