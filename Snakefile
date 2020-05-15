import platform

rule all:
  input: 'index/_book/dissertation.pdf'


rule install_deps:
  conda: 'envs/R.yml'
  params:
    tar = '/usr/bin/tar' if platform.system() == 'Darwin' else '/bin/tar'
  shell: """
    TAR={params.tar} R -e 'devtools::install_github("ryanpeek/aggiedown@ae99300d43bdccc16069efcc08198624c76eee0c", upgrade = "never")'
  """


rule build_thesis:
  conda: 'envs/R.yml'
  output: 'index/_book/dissertation.pdf'
  input:
    sources=expand('index/{rmd}.Rmd',
                   rmd=('index', '00-intro', '01-chap1', '02-chap2',
                        '03-chap3', '04-conclusion', '05-appendix',
                        '98-colophon', '99-references')),
    bibliography='index/bib/dissertation.bib'
  shell: """
      cd index 
      rm -f _main.Rmd
      R -e "options(tinytex.verbose = TRUE); bookdown::render_book('index.Rmd', aggiedown::thesis_pdf(latex_engine = 'xelatex'))"
      mv _book/_main.pdf _book/dissertation.pdf
  """
