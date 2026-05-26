# .Rprofile for the UAS Gauss OSSE publication companion
# Sets project root and ensures consistent behaviour across sessions.

options(uas_gauss_osse_root = normalizePath("."))
options(stringsAsFactors = FALSE)
options(datatable.print.class = TRUE)

message("UAS Gauss OSSE project loaded. Root: ", getwd())

# Check for required packages on startup
.check_packages <- function() {
  required <- c("data.table", "ggplot2", "viridis", "cowplot")
  optional <- c("scales", "progressr")

  missing_req <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  missing_opt <- optional[!sapply(optional, requireNamespace, quietly = TRUE)]

  if (length(missing_req) > 0) {
    message("MISSING required packages: ", paste(missing_req, collapse = ", "))
    message("Install with: install.packages(c(",
            paste0('"', missing_req, '"', collapse = ', '), "))")
  }
  if (length(missing_opt) > 0) {
    message("Optional packages not found: ", paste(missing_opt, collapse = ", "))
  }
}

.check_packages()
rm(.check_packages)
