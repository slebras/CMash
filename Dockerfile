FROM python:3.5.7-slim-jessie
COPY . /app
WORKDIR /app
RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get install -y gcc
RUN pip install --upgrade pip && \
    pip install -r requirements.txt && \
    pip install requests
