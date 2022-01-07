import platform

rule all:
  input: 'index/_book/dissertation.pdf'


rule install_deps:
  conda: 'envs/R.yml'
  params:
    tar = '/usr/bin/tar' if platform.system() == 'Darwin' else '/bin/tar'
  shell: """
    TAR={params.tar} R -e 'devtools::install_github("camillescott/aggiedown@b7e13aa079afe470d1002ba18ffe6fd5c384cd4a", upgrade = "never")'
  """


rule build_thesis:
  conda: 'envs/R.yml'
  output: 'index/_book/dissertation.pdf'
  input:
    sources=expand('index/{rmd}.Rmd',
                   rmd=('index',
                        '00-intro',
                        '01-chap1',
                        '02-chap2',
                        '03-chap3',
                        '04-shmlast',
                        '05-suchtree',
                        #'06-workflows',
                        '09-conclusion',
                        '10-appendix',
                        '98-colophon',
                        '99-references')),
    bibliography='index/bib/dissertation.bib'
  shell: """
      cd index 
      rm -f _main.Rmd
      R -e "options(tinytex.verbose = TRUE); bookdown::render_book('index.Rmd', output_format='aggiedown::thesis_pdf')"
      mv _book/_main.pdf _book/dissertation.pdf
      cp _book/dissertation.pdf ../
  """
