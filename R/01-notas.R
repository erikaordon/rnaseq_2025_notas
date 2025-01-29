## Configura tu usuario de GitHub
usethis::create_project("~/rnaseq_2025_notas")
usethis::use_r("01-notas.R")
usethis::create_github_token()
gitcreds_set()

usethis::edit_git_config()
usethis::use_github()
#Queremos usar git con nuestro proyecto

#comentario
