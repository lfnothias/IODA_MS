FROM openms/executables:release-2.5.0

# Installing core dependencies for python
RUN apt-get update -y
RUN apt-get install -y python3-pip python3 build-essential

# Installing Virtual Frame Buffer
RUN apt-get install -y xvfb
RUN apt-get install -y wget

# Binder Specific Installations
RUN pip3 install --no-cache-dir notebook==5.*
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
WORKDIR /home/${NB_USER}

# Copying files
RUN echo ${HOME}
COPY . ${HOME}

RUN pwd ${HOME}
RUN ls ${HOME}
RUN pip3 install -r ${HOME}/requirements.txt

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
