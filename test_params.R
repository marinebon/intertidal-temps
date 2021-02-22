purrr::walk(
  seq(5,25,by=5),
  function(x){
    rmarkdown::render(
      "test_params.Rmd", 
      params = list(
        n_bins = x), 
      output_file = glue::glue("test_params_{x}.html"))
  })