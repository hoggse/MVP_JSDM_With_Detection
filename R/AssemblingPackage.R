##############################
#
# package building


require(usethis)
require(devtools)
require(roxygen2)

available.packages(MVPJSDMreal)
available("MVPJSDMreal")

use_package("psych")
use_package("jagsUI")
use_package("ROCR")

use_readme_rmd()
