


### Download RAP SLEEP

cats <- c(207:216)

for (c in cats) {
  fields <- data.table::fread("data/field.tsv") %>%
    filter(main_category == c)

  # Build per-field arguments: always include eid first
  field_args <- c("-ifield_names=eid")
  if (nrow(fields) > 0) {
    field_args <- c(field_args, paste0("-ifield_names=p", fields$field_id))
  }

  # Build the full command as a single string
  cmd_parts <- c(
    "dx run table-exporter",
    "-idataset_or_cohort_or_dashboard=record-Gp1BZyjJY95YyKkqg1XqQ70f",
    paste0("-ioutput=sleep_", c),
    "-ientity=participant",
    "--tag=table_exporter",
    "--name=table_exporter",
    "--brief -y",
    # append all -ifield_names arguments (each a separate flag)
    field_args
  )

  cmd <- paste(cmd_parts, collapse = " ")

  cat(cmd, "\n\n")   # print the command
  system(cmd)
}


paste(c("dx run table-exporter -idataset_or_cohort_or_dashboard=record-Gp1BZyjJY95YyKkqg1XqQ70f -ioutput=physical  -ientity=participant --tag=table_exporter  --name=table_exporter --brief -y ",
        paste("-ifield_names", c("eid", "p4080_i0_a0", "p4080_i0_a1", "p4079_i0_a0", "p4079_i0_a1", "p46_i0", "p47_i0", "p20023_i0", "p20016_i0", "p20532", "p20534", "p20533", "p20535"), sep = "=")), collapse = " ")

