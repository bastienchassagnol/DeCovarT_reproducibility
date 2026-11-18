# preprocessCore (DeMixSC dependency) fails to build with threading on this
# machine (pthread_create() is 22); renv::restore() then needs this flag,
# which renv.lock cannot record.
options(configure.args = list(preprocessCore = "--disable-threading"))

# renv autoloader (created after renv::init())
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}
