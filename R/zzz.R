###############################################
#### AUTHOR:    Arnost Komarek             ####
####            (2003)                     ####
####                                       ####
#### FILE:      zzz.R                      ####
####                                       ####
#### FUNCTIONS: .First.lib                 ####
###############################################

### =============================================
### .First.lib
### =============================================
.First.lib <- function(lib, pkg)
{
   require(survival)
   library.dynam("smoothSurv", pkg, lib)

   invisible()
}

