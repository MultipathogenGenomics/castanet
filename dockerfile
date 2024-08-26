FROM continuumio/miniconda3

ENV CACHE_TTL=15
ENV CACHE_MAX_SIZE=1024

EXPOSE 80

COPY ./app /workspace/app
COPY ./install_deps.sh /workspace/
COPY ./requirements.txt /workspace/
COPY ./data/eval /workspace/data/eval

WORKDIR /workspace
RUN mkdir /logs

RUN conda create -y -n castanet python=3.10
RUN echo "source activate castanet" > ~/.bashrc
ENV PATH /opt/conda/envs/castanet/bin:$PATH

RUN pip install --upgrade pip
RUN pip install -r requirements.txt --prefer-binary
RUN bash install_deps.sh

# CMD ["python3", "-m", "app.api"]

# docker run --name castanet -p 8001:80
