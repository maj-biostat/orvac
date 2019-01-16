# readme.md

This is the orvac simulation code. In a nutshell, use `run_sim_1.sh` to run simulations. You may well need to edit the script.

`main_2.R` is the entry point.

When running simulations, the following may be useful (run from the logs directory) to keep an eye on progress. You might also want to run `htop` to ensure you are using all cores.

`watch 'grep  -ie warn -ie "Cluster stopped" *.log'`

Results (from selected output, i.e. you might need to tweak) is obtained by rendering from `simulations_report.Rmd`. To get a `html` version, from `R` do `rmarkdown::render("simulation_report.Rmd", clean=TRUE)`.
