######################
## package: ssimf ###
#####################

#inicia o pacote
#install.packages("devtools")
devtools::setup(rstudio = FALSE)

#preencher o DESCRIPTION

# crie a pasta \inst e salve o arquivo octa


# cria pasta para dados brutos
devtools::use_data_raw()

#salve este arquivo em data-raw

# Ignora essa rotina em data-raw
devtools::use_build_ignore("data-raw")

# Ignora Rproj do Rstudio
devtools::use_build_ignore("ssimf.Rproj")

#Escrevas as funcoes e salve em R

# documenta as funcoes
devtools::document()

# coloca as dependencias no pacote
devtools::use_package("RcppOctave")

# testa o pacote, provavelmente recebera um erro
devtools::check()

# corriga os possiveis erros ou adivertencias apontados 

# teste o pacote novamente apos correcao dos erros
devtools::document()
devtools::check()



## corrige os acentos
tools::showNonASCII(readLines("R/deemd.R"))

# teste o pacote novamente
devtools::check()

# instala
devtools::install()

# apenas carrega o pacote
devtools::load_all()


