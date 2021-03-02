# build environment
# FROM node:12-alpine as react-build
# WORKDIR /app
# COPY . ./
# ENV CONFIGFILE config/config_genbank.yaml
# RUN npm install
# RUN npm run build-only

# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.8-slim

ARG CONFIGFILE

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True

WORKDIR /opt
ADD $CONFIGFILE ./config.yaml
ADD ./static_data /static_data
ADD ./src/constants/defs.json ./defs.json

ADD ./services/server/requirements.txt ./requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt