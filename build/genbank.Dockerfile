# build environment
FROM node:12-alpine as react-build
WORKDIR /app
COPY . ./
ENV CONFIGFILE config/config_genbank.yaml
RUN npm install
RUN npm run build-only

# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.8-slim

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True

# Copy local code to the container image.
ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

# Install production dependencies.
# RUN pip install Flask gunicorn
RUN pip install -r flask_server/requirements.txt

COPY --from=react-build /app/dist ./flask_server/dist

ENV FLASK_ENV production
ENV PORT 8080
ENV CONFIGFILE /app/config/config_genbank.yaml

# Run the web service on container startup. Here we use the gunicorn
# webserver, with one worker process and 8 threads.
# For environments with multiple CPU cores, increase the number of workers
# to be equal to the cores available.
CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 flask_server.main:app