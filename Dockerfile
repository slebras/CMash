FROM python:3.7-slim
COPY . /app
WORKDIR /app
RUN apt-get update -y && \
    apt-get install -y build-essential && \
    apt-get install gcc && \
    pip install --upgrade pip && \
    pip install -r requirements.txt

