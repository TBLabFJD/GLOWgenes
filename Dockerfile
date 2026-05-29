FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    r-base \
    wget \
    unzip \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN pip3 install numpy pandas scipy scikit-learn networkx matplotlib

RUN R -e "install.packages(c('optparse','caret'), repos='https://cloud.r-project.org/')"

# Clonar GLOWgenes y corregir compatibilidad con sklearn moderno
RUN git clone https://github.com/TBLabFJD/GLOWgenes.git /opt/GLOWgenes \
    && sed -i 's/metrics.fbeta_score(lab, pred_label, 2)/metrics.fbeta_score(lab, pred_label, beta=2)/' /opt/GLOWgenes/GLOWgenes.py

# Descargar las redes desde figshare
RUN curl -L -o /opt/GLOWgenesNets.zip \
    https://github.com/TBLabFJD/GLOWgenes/releases/download/v1.0/GLOWgenesNets.zip \
    && unzip /opt/GLOWgenesNets.zip -d /opt/GLOWgenesNets \
    && rm /opt/GLOWgenesNets.zip \
    && sed -i 's|/home/proyectos/bioinfo/lorena/pSNOW/interactomes/networks_cfg/pSNOW_databases_new_nets/|/opt/GLOWgenesNets/GLOWgenesNets/|g' \
       /opt/GLOWgenesNets/GLOWgenesNets/networks_knowledgeCategories.cfg
       
       
WORKDIR /workspace
CMD ["/bin/bash"]
