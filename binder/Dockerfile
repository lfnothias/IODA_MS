FROM openms/library:release-2.5.0


# Installing core dependencies for python
RUN apt-get update -y
RUN apt-get install -y python3-pip python3 build-essential

# Building OpenMS
WORKDIR /openms-build

# configure
RUN /bin/bash -c "cmake -DCMAKE_PREFIX_PATH='/contrib-build/;/usr/;/usr/local' -DBOOST_USE_STATIC=OFF -DWITH_GUI=ON ../OpenMS"

# make OpenMS library
# RUN make OpenMS -j7

# Linking all executables
# RUN make -j7

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