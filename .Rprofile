# preprocessCore (DeMixSC dependency) fails to build with threading on this
# machine (pthread_create() is 22); renv::restore() then needs this flag,
# which renv.lock cannot record.
options(configure.args = list(preprocessCore = "--disable-threading"))

# Bioconductor 3.23 is already fixed in renv.lock. Giving renv the
# repository list here means activation does not load BiocManager from
# the cache before the project library is active. The CRAN URL is the
# lockfile URL; renv rewrites it to the platform binary at install time.
options(renv.bioconductor.repos = c(
  CRAN = "https://packagemanager.posit.co/cran/latest",
  BioCsoft = "https://bioconductor.org/packages/3.23/bioc",
  BioCann = "https://bioconductor.org/packages/3.23/data/annotation",
  BioCexp = "https://bioconductor.org/packages/3.23/data/experiment",
  BioCworkflows = "https://bioconductor.org/packages/3.23/workflows",
  BioCbooks = "https://bioconductor.org/packages/3.23/books"
))

# renv autoloader (created after renv::init())
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}
