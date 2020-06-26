FROM openms/library:release-2.5.0

# Installing core dependencies for python
RUN apt-get update -y
RUN apt-get install -y python3-pip python3 build-essential

RUN pip3 install ftputil
RUN pip3 install flask
RUN pip3 install gunicorn
RUN pip3 install requests
RUN pip3 install celery
RUN pip3 install joblib
RUN pip3 install vladiate
RUN pip3 install pandas
RUN pip3 install matplotlib

# Building OpenMS
WORKDIR /
RUN git clone -b Release2.4.0 --single-branch https://github.com/OpenMS/OpenMS.git && cd /OpenMS

RUN mkdir openms-build
WORKDIR /openms-build

# configure
RUN /bin/bash -c "cmake -DCMAKE_PREFIX_PATH='/contrib-build/;/usr/;/usr/local' -DBOOST_USE_STATIC=OFF -DWITH_GUI=ON ../OpenMS"

# make OpenMS library
RUN make OpenMS -j7

# Linking all executables
RUN make -j7

# Installing Virtual Frame Buffer
RUN apt-get install -y xvfb


# Binder Specific Installations
RUN pip install --no-cache-dir notebook==5.*
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}