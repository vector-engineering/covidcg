COVID-19 CG SERVER
------------------

Since the app is all front-end, this server is just a file server with HTTP authentication

The following files were written and run on Debian 10 (buster) on a Google Compute Engine VM (1 vCPU, 600MB RAM).

Setup
-----

Edit the covid.service file to point to the server/server.py script and the dist/ directory with the built app, on your machine.

All the following steps require root privileges.

 - Move the covid.service file into /etc/systemd/system/covid.service
 - Enable the service with: systemctl daemon-reload && systemctl enable covid && systemctl start covid --no-block
 - View status with: systemctl status covid
 - Stop the service: systemctl stop covid

Authentication
--------------

The HTTP user:password pair is defined in covid.service, and passed to the python server script

Attributions
------------

 - User pztrick on https://stackoverflow.com/questions/16420092/how-to-make-python-script-run-as-service
 - User tianhuil on GitHub: https://github.com/tianhuil/SimpleHTTPAuthServer