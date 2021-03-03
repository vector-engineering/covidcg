# build environment
FROM node:12-alpine as react-build
ARG CONFIGFILE
WORKDIR /app
COPY . ./
ENV CONFIGFILE $CONFIGFILE
RUN npm install
RUN npm run build-only

# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.8-slim
ARG CONFIGFILE

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True

# Copy local code to the container image.
WORKDIR /opt
COPY $CONFIGFILE ./config.yaml
COPY ./static_data /static_data
COPY ./src/constants/defs.json ./defs.json

# Copy local code to the container image.
COPY ./services/server /app
WORKDIR /app

# Install production dependencies.
RUN pip install -r ./requirements.txt

COPY --from=react-build /app/dist ./cg_server/dist
# COPY ./dist ./dist

ENV FLASK_ENV production
ENV PORT 8080
ENV CONFIGFILE /opt/config.yaml

# Run the web service on container startup. Here we use the gunicorn
# webserver, with one worker process and 8 threads.
# For environments with multiple CPU cores, increase the number of workers
# to be equal to the cores available.
CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 cg_server.main:app