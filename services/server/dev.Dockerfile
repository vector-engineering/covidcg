# --------------------
# FOR DEVELOPMENT ONLY
# --------------------

# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.8-slim

ENV PYTHONDONTWRITEBYTECODE 1
# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True

# Install dependencies
# We need to move the requirements.txt folder over manually.
# Although the services/server folder will be mounted at /app,
# this only happens at run-time, and we need the requirements
# file *now* during build-time.
WORKDIR /opt
ADD ./services/server/requirements.txt ./requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
