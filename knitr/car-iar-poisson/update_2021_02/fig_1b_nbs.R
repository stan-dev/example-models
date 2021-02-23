#spdep nb object list of integer vectors
fig_1b_nb = list(
as.vector(c(2,3), mode="integer"),       # 1
as.vector(c(1,3), mode="integer"),       # 2
as.vector(c(1,2), mode="integer"),   # 3
as.vector(c(5), mode="integer"),       # 4
as.vector(c(4), mode="integer"),       # 5
as.vector(c(0), mode="integer"))       # 6

attr(fig_1b_nb, "class") <- "nb"

attr(fig_1b_nb, "region.id") <-
c("1", "2", "3", "4", "5", "6")

