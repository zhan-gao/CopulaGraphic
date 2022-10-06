kendall_alpha_joe <- read_excel("data/kendall_alpha_joe.xlsx",
                                col_names = FALSE)
colnames(kendall_alpha_joe) <- c("alpha", "ktau")
usethis::use_data(kendall_alpha_joe)
