install.packages("renv")
renv::restore()  # Reads renv.lock and installs everything

renv::activate(project = "/mnt/project/renv/")
