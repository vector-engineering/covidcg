# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.8-slim

ARG CONFIGFILE

ENV PYTHONDONTWRITEBYTECODE 1
# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED True

# Copy local code to the container image.
WORKDIR /opt
# COPY . ./

ADD $CONFIGFILE ./config.yaml
ADD ./static_data /static_data
ADD ./src/constants/defs.json ./defs.json

# Install production dependencies.
# RUN pip install Flask gunicorn
ADD ./services/server/requirements.txt ./requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt



# COPY --from=react-build /app/dist ./flask_server/dist

# Run the web service on container startup. Here we use the gunicorn
# webserver, with one worker process and 8 threads.
# For environments with multiple CPU cores, increase the number of workers
# to be equal to the cores available.
# CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 cg_server.main:app