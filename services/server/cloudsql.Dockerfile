# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.8-slim

ENV PYTHONDONTWRITEBYTECODE 1
# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True

WORKDIR /opt
ADD ./services/server/requirements.txt ./requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt