FROM openms/executables:release-2.5.0

# Installing core dependencies for python
RUN apt-get update -y
RUN apt-get install -y python3-pip python3 build-essential

# Installing Requirements
COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN pip3 install notebook
RUN pip3 install pandas
RUN pip3 install scipy
RUN pip3 install numpy
RUN pip3 install matplotlib
RUN pip3 install logzero
RUN pip3 install xvfbwrapper
RUN pip3 install jupyterlab

# Installing Virtual Frame Buffer
RUN apt-get install -y xvfb
RUN apt-get install -y wget

# Copying files
COPY . .
CMD ["bash", "-c", "source /etc/bash.bashrc && jupyter lab --notebook-dir=/ --ip 0.0.0.0 --no-browser --allow-root --NotebookApp.token='' --NotebookApp.password=''"]
