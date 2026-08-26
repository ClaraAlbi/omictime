library(stringr)
library(tibble)

extract_ldsc <- function(f) {
  
  fname <- basename(f)
  stem  <- sub("\\.log$", "", fname)
  parts <- strsplit(stem, "_vs_")[[1]]
  
  p1_name <- parts[1]
  p2_name <- parts[2]
  
  x <- readLines(f, warn = FALSE)
  x_trim <- str_trim(x)
  
  ## -----------------------------
  ## Extract heritability blocks
  ## -----------------------------
  
  # Phenotype 1 h2
  h2_1_line <- x_trim[
    grepl("^Total Observed scale h2:", x_trim) &
      cumsum(grepl("^Heritability of phenotype 1", x_trim)) == 1
  ][1]
  
  # Phenotype 2 h2
  h2_2_line <- x_trim[
    grepl("^Total Observed scale h2:", x_trim) &
      cumsum(grepl("^Heritability of phenotype 2", x_trim)) == 1
  ][1]
  
  parse_h2 <- function(line) {
    if (is.na(line)) return(c(NA, NA))
    nums <- str_match(line, "h2:\\s*([0-9.]+)\\s*\\(([0-9.]+)\\)")
    as.numeric(nums[2:3])
  }
  
  h2_1 <- parse_h2(h2_1_line)
  h2_2 <- parse_h2(h2_2_line)
  
  ## -----------------------------
  ## Extract rg line
  ## -----------------------------
  
  rg_line <- x_trim[
    grepl("/well", x_trim) &
      grepl("\\s-?[0-9]+\\.[0-9]+\\s+[0-9]+\\.[0-9]+", x_trim)
  ][1]
  
  if (is.na(rg_line)) {
    return(tibble(
      file = fname,
      p1_name = p1_name,
      p2_name = p2_name,
      h2_p1 = h2_1[1],
      h2se_p1 = h2_1[2],
      h2_p2 = h2_2[1],
      h2se_p2 = h2_2[2],
      rg = NA, se = NA, z = NA, p = NA
    ))
  }
  
  parts_line <- str_split(rg_line, "\\s+", simplify = TRUE)
  
  tibble(
    file = fname,
    p1_name = p1_name,
    p2_name = p2_name,
    h2_p1   = h2_1[1],
    h2se_p1 = h2_1[2],
    h2_p2   = h2_2[1],
    h2se_p2 = h2_2[2],
    rg      = as.numeric(parts_line[3]),
    se      = as.numeric(parts_line[4]),
    z       = as.numeric(parts_line[5]),
    p       = as.numeric(parts_line[6])
  )
}
