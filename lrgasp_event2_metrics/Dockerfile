FROM condaforge/mambaforge-pypy3

WORKDIR /app

# Create the environment:
COPY sqanti3_lrgasp.yml .
RUN mamba env create -n sqanti_env -f sqanti3_lrgasp.yml
SHELL ["conda", "run", "-n", "sqanti_env", "/bin/bash", "-c"]

# Download and set-up cDNA_Cupcake
RUN apt-get update && wget https://github.com/Magdoll/cDNA_Cupcake/archive/refs/tags/v24.0.0.tar.gz && tar -xvf v24.0.0.tar.gz
RUN mv cDNA_Cupcake-24.0.0 cDNA_Cupcake
COPY . /app
RUN chmod +x utilities/gmst/Gibbs3 utilities/gmst/gmhmmp utilities/gmst/probuild utilities/gtfToGenePred
ENV PYTHONPATH "${PYTHONPATH}:/app/cDNA_Cupcake/sequence/"
ENV PYTHONPATH "${PYTHONPATH}:/app/cDNA_Cupcake/"

RUN apt-get install libcurl4 libcurl4-openssl-dev -y

